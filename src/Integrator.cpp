#include "../include/Body.h"
#include "../include/Forces.h"
#include "../include/Context.h"
#include "mpi.h"
#include <vector>


/*
  We accumulate here a list of integration schemes.  These solve 
  the system of ODEs:

  x' = v(x,t)
  v' = a(x,t)

  using various integration schemes.  The intention is to include the
  methods used by Volker Springel in Gadget to see if Adams-Bashforth
  methods may be better choices than the Leapfrog methods he uses.  
  We also verify the functionality of N-TR0PY by implementing these
  schemes and comparing to known results for the 2-body problem.  

  Right now, we intend to implement Runge-Kutta methods, the two Leapfrog
  methods Springel describes, and a series of Adams-Bashforth methods.  
  Adams-Bashforth methods are explicity s-step methods which maximize 
  order (accuracy) and stability.  For more than a 1-step method, we 
  generate the ICs via another integration scheme.   The idea is to 
  avoid making unnecessary calls to the force calculator which can 
  run in as bad as N^2 time.  It may be possible to halve the time that
  Springel's code computes and updates positions and velocities.  

  Results, if positive, will be be published in my upcoming master's
  thesis this summer (2015).  
*/

using namespace std;


/*
  The Euler Method is included for instructional purposes and for calculating 
  ICs only.  It is horrible in every possible way.  It is low order (1) and 
  is laughably unstable.  
*/

int EulerMethod(vector<bodies_t> &bodies, double evolveTime, int numberOfSteps, context_t &NBODY_CONTEXT)
{
  int nthreads, rank;

  MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  int N = bodies.size();
  vector<double> emptyForce(3,0.0);

  double h = evolveTime/(double)numberOfSteps;

  vector<vector<double> > forces(N,emptyForce);

  // Minor efficiency loss for calculating these values on each thread                                                                            

  int remainder = N % nthreads;
  int block,n0,n1;

  //double buffer[n1-n0];

  if (nthreads == 1)
    {
      n1 = N-1;
      n0 = 0;
    }

  else
    {
      block = int((N - remainder)/nthreads);


      // partitions the body list.                                                                                                                
      if (rank == (nthreads - 1))
        {
          n0 = block * rank;
          n1 = remainder + block*(rank);
        }

      else
        {
          n0 = block * rank;
          n1 = n0 + block * (rank);
        }

    }



  /*
    In principle, this process could be easily parallelized.  
    It is, however, an O(N) process and for debugging purposes,
    it will be temporarily left serial.  The force calculations,
    which are at least O(N log(N)) are left parallel.  
  */
  
  double currentTime = 0;

  
  while (currentTime < evolveTime)
    {
      N2BruteForce(bodies, forces, NBODY_CONTEXT.init3Volume, NBODY_CONTEXT.eps2);
      for (int i = n0; i<n1; i++)
	{
	  bodies[i].u1 += 0.9785 * forces[i][0] * h; //Convert u to km/s
	  bodies[i].u2 += 0.9785 * forces[i][1] * h;
	  bodies[i].u3 += 0.9785 * forces[i][2] * h;

	  forces[i][0] = 0.0;
	  forces[i][1] = 0.0;
	  forces[i][2] = 0.0;
	  
	  bodies[i].q1 += 1.022 * bodies[i].u1 * h; // convert u to kpc/Gy
	  bodies[i].q2 += 1.022 * bodies[i].u2 * h;
	  bodies[i].q3 += 1.022 * bodies[i].u3 * h;

	  //printf("[%10.5f, %10.5f, %10.5f]\n", bodies[i].q1, bodies[i].q2, bodies[i].q3);
	}

      currentTime += h;
      MPI_Barrier(MPI_COMM_WORLD);

      // Now we must broadcast and synchronize the list
      // as separate components and reconstruct them in the master
      // thread.   Then we must broadcast them back to each thread.
      // This is potentially costly, but again is an O(N) calculation
      // which should scale nicely with the number of threads.  

      // Convert changed components to C arrays whose transport is supported
      // by MPI.  The order will be preserved because the master root will 
      // take each thread in order.  Because of how the load was balanced,
      // we need not worry about particle IDs when sorting this way.  

      double q1list[n1-n0];
      double q2list[n1-n0];
      double q3list[n1-n0];
      
      double u1list[n1-n0];
      double u2list[n1-n0];
      double u3list[n1-n0];

      for (int i = 0; i < (n1-n0); i++)
	{
	  q1list[i] = bodies[i].q1;
          q2list[i] = bodies[i].q2;
          q3list[i] = bodies[i].q3;

          u1list[i] = bodies[i].u1;
          u2list[i] = bodies[i].u2;
          u3list[i] = bodies[i].u3;
	}
      
      MPI_Send(q1list,(n1-n0),MPI_DOUBLE,0,1,MPI_COMM_WORLD);
      MPI_Send(q2list,(n1-n0),MPI_DOUBLE,0,2,MPI_COMM_WORLD);
      MPI_Send(q3list,(n1-n0),MPI_DOUBLE,0,3,MPI_COMM_WORLD);
      
      MPI_Send(u1list,(n1-n0),MPI_DOUBLE,0,4,MPI_COMM_WORLD);
      MPI_Send(u2list,(n1-n0),MPI_DOUBLE,0,5,MPI_COMM_WORLD);
      MPI_Send(u3list,(n1-n0),MPI_DOUBLE,0,6,MPI_COMM_WORLD);
	
      MPI_Barrier(MPI_COMM_WORLD);

      if (rank == 0)
	{
	  int expected;
	  for (int i = 1; i<nthreads; i++)
	    {
	      
	      if(i != (nthreads - 1))
		{
		  expected = block;
		}
	      else
		{
		  expected = block + remainder;
		}
	      for (int j = 0; j < 6; j++)
		{
		  MPI_Status * status;
		  double temp[expected]; 
		  MPI_Recv(temp, expected, MPI_DOUBLE, i, j, MPI_COMM_WORLD, status);

		  

		  // Should the memory for the temporary variables be freed?
		  //allocator::deallocate(*status,sizeof(status));
		  //allocator::deallocate(*temp,sizeof(temp));
		}
	    }
	}


    }
  return 0;
}
