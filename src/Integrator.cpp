#include "../include/Body.h"
#include "../include/Forces.h"
#include "../include/Context.h"
#include "mpi.h"
#include <vector>
#include <cstdlib>
#include <cmath>
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


// NOTE: PERHAPS ONLY MAKE FORCE CALCUALTION PARALLEL?


/*
  The Euler Method is included for instructional purposes and for calculating 
  ICs only.  It is horrible in every possible way.  It is low order (1) and 
  is laughably unstable.  
*/

vector<bodies_t> EulerMethod(vector<bodies_t> bodies, double evolveTime, int numberOfSteps, context_t &NBODY_CONTEXT)
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
      n1 = N;
      n0 = 0;
    }

  else
    {
      block = int((N - remainder)/nthreads);


      // partitions the body list.                                                                                                                
      if (rank == (nthreads - 1))
        {
          n0 = block * rank;
          n1 = n0 + remainder + block*(rank);
        }

      else
        {
          n0 = block * rank;
          n1 = n0 + block * (rank + 1);
        }

    }


  //printf("I am thread %d and my (n0, n1) is (%d, %d)!\n", rank, n0, n1);
  /*
    In principle, this process could be easily parallelized.  
    It is, however, an O(N) process and for debugging purposes,
    it will be temporarily left serial.  The force calculations,
    which are at least O(N log(N)) are left parallel.  
  */
  
  double currentTime = 0;

  MPI_Datatype MPI_BODY = createMPIBody();

  //printf("MPI DATATYPE DEFINED.\n");

  //printf("MPI AND C SIZES ARE %d and %d.\n",(int)sizeof(MPI_BODY), (int)sizeof(bodies_t));

  while (currentTime < evolveTime)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      bodies_t  bodiesArr[N];
      for (int i=0; i<N; i++)
	{
	  bodies_t body;
	  body.id = -1;
	  body.type = -1;
	  body.mass = -1.0;
	  body.u1 = -1.0;
	  body.u2 = -1.0;
	  body.u3 = -1.0;
	  body.q1 = -1.0;
	  body.q2 = -1.0;
	  body.q2 = -1.0;
	  bodiesArr[i] = body;
	}

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

	  //printf("The calculated data is:\n");
	  //printf("r = [%5.10f, %5.10f, %5.10f]\n", bodies[i].q1,bodies[i].q2,bodies[i].q3);
	  //printf("u = [%5.10f, %5.10f, %5.10f]\n", bodies[i].u1,bodies[i].u2,bodies[i].u3);

	  bodiesArr[i] = bodies[i];

	  //printf("[%10.10f, %10.10f, %10.10f]\n", bodies[i].q1, bodies[i].q2, bodies[i].q3);
	}

      // printf("FORCES CALCULATED ON THREAD %d.\n",rank);
      

      //bodies_t * bodiesArr[n1 - n0]
      // MPI_Datatype MPI_BODY = createMPIBody();

      //printf("MPI DATATYPE DEFINED.\n");

      if (rank != 0)
	{
	  //printf("Preparing to send data of size %d from thread %d...\n",N,rank);
	  MPI_Send(&bodiesArr, N, MPI_BODY, 0, rank, MPI_COMM_WORLD);
	  //printf("Sent data from thread %d.\n",rank);
	}

      if (rank == 0)
	{
	  int Recvn0, Recvn1, expected;
	  //printf("PERFORMING SERIAL OPERATIONS.\n");
     
	  for(int j = 1; j<nthreads; j++)
	    {
	      if (rank == (nthreads - 1))
		{
		  Recvn0 = block * j;
		  Recvn1 =  remainder + (block + 1)*j;
		  // printf("[Recvn1, Recvn2] = [%d,%d]",Recvn0,Recvn1);
		}

	      else
		{
		  Recvn0 = block * j;
		  Recvn1 = Recvn0 + block * (j );
		  // printf("[Recvn1, Recvn2] = [%d,%d]", Recvn0,Recvn1);

		}

	      expected = Recvn1 - Recvn0;
	      bodies_t tempBody[N];
	      MPI_Status status;
	      MPI_Request  request;
	      int  flag;
	      
	      //printf("Preparing to receive data of size %d from thread %d...\n",N,j);
	      int ret = MPI_Recv(&tempBody, N, MPI_BODY, j,j,MPI_COMM_WORLD,&status);
	      //printf("Received value staus is %d\n",ret);
	      //printf("Data received.");
	      
	      for(int k = Recvn0; k < Recvn1; k++ )
		{
		  //printf("K = %d!\n",k);
		  //printf("TempBody %d has RECEIVED values:\n",k);
		  //printf("r = [%5.10f, %5.10f, %5.10f]\n", tempBody[k].q1,tempBody[k].q2,tempBody[k].q3);
		  //printf("u = [%5.10f, %5.10f, %5.10f]\n", tempBody[k].u1,tempBody[k].u2,tempBody[k].u3);
		  bodiesArr[k] = tempBody[k];
		}

	      //MPI_Bcast(bodies,)
	      //MPI_Test(request,flag,status);
	    }
	  //for (int m = 0; m < N; m++)
	  // {
	  //  bodies[m] = bodiesArr[m];
	  //}

	  // MPI_Test(status,request);
	  //printf("Beginning serial code.\n");
	}
      MPI_Barrier(MPI_COMM_WORLD);

      MPI_Bcast(bodiesArr, N, MPI_BODY,0, MPI_COMM_WORLD );
      for (int m = 0; m < N; m++)
	{
	  bodies[m] = bodiesArr[m];
	}

      //MPI_Bcast(bodiesArr, N, MPI_BODY,0, MPI_COMM_WORLD );

      currentTime += h;
      //      printf("Current time is %10.10f.\n",currentTime);

      if(rank == 0)
	{
	  for (int i = 0; i < N; i++)
	    {
	      //printf("Current time is %10.10f.\n",currentTime);
	      //printf("The data for body %d is:\n",i);
	      //printf("r = [%5.10f, %5.10f, %5.10f]\n", bodies[i].q1,bodies[i].q2,bodies[i].q3);
	      //printf("u = [%5.10f, %5.10f, %5.10f]\n", bodies[i].u1,bodies[i].u2,bodies[i].u3);
	      printf("%5.10f\n", pow(pow(bodies[i].q1,2) + pow(bodies[i].q2, 2) + pow(bodies[i].q3,2),0.5));
	    }
	}

      
      MPI_Barrier(MPI_COMM_WORLD);
      

      //MPI_Scatter(Q1, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

     //MPI_Bcast(Q1,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
  return bodies;
    
}
