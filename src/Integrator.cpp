#include "../include/Body.h"
#include "../include/Forces.h"
#include "../include/Context.h"
#include "mpi.h"
#include <vector>
#include <cstdlib>

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
      MPI_Barrier(MPI_COMM_WORLD);
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

	  //printf("[%10.10f, %10.10f, %10.10f]\n", bodies[i].q1, bodies[i].q2, bodies[i].q3);
	}

      currentTime += h;
      printf("Current time is %10.10f.\n",currentTime);
      MPI_Barrier(MPI_COMM_WORLD);
      printf("Thread %d passed barrier 1.\n", rank);
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

      MPI_Status send_status;
      MPI_Request send_request;
      printf("Thead %d preparing to send %d elements...\n",rank, (n1-n0));
      
      if (rank != 0)
	{
	  MPI_Send(q1list, (n1-n0), MPI_DOUBLE, 0, 1 + 6 * rank, MPI_COMM_WORLD);
          MPI_Send(q2list, (n1-n0), MPI_DOUBLE, 0, 2 + 6 * rank, MPI_COMM_WORLD);
          MPI_Send(q3list, (n1-n0), MPI_DOUBLE, 0, 3 + 6 * rank, MPI_COMM_WORLD);
          MPI_Send(u1list, (n1-n0), MPI_DOUBLE, 0, 4 + 6 * rank, MPI_COMM_WORLD);
          MPI_Send(u2list, (n1-n0), MPI_DOUBLE, 0, 5 + 6 * rank, MPI_COMM_WORLD);
          MPI_Send(u3list, (n1-n0), MPI_DOUBLE, 0, 6 + 6 * rank, MPI_COMM_WORLD);

	}

      //MPI_Isend(q1list,(n1-n0),MPI_DOUBLE,0,1 + 6 * rank,MPI_COMM_WORLD,&send_request);
      //MPI_Isend(q2list,(n1-n0),MPI_DOUBLE,0,2 + 6 * rank,MPI_COMM_WORLD,&send_request);
      //MPI_Isend(q3list,(n1-n0),MPI_DOUBLE,0,3 + 6 * rank,MPI_COMM_WORLD,&send_request);      
      //MPI_Isend(u1list,(n1-n0),MPI_DOUBLE,0,4 + 6 * rank,MPI_COMM_WORLD,&send_request);
      //MPI_Isend(u2list,(n1-n0),MPI_DOUBLE,0,5 + 6 * rank,MPI_COMM_WORLD,&send_request);
      //MPI_Isend(u3list,(n1-n0),MPI_DOUBLE,0,6 + 6 * rank,MPI_COMM_WORLD,&send_request);
	
      //MPI_Wait(&send_request, &send_status);
      //MPI_Barrier(MPI_COMM_WORLD);
      
      printf("Data sent by thread %d.\n",rank);
      
      // Synchronize before serial operations

      double  Q1[N];
      double  Q2[N];
      double  Q3[N];

      double  U1[N];
      double  U2[N];
      double  U3[N];

      MPI_Barrier(MPI_COMM_WORLD);

      if (rank == 0)
	{
	  int block = int((N - remainder)/nthreads);

	  for (int k = 0; k < block; k++)
	    {
	      Q1[k] = q1list[k];
	      Q2[k] = q2list[k];
	      Q3[k] = q3list[k];
	      U1[k] = u1list[k];
	      U2[k] = u2list[k];
	      U3[k] = u3list[k];
	    }
	
	  printf("Performing serial operations.\n");
	  int expected;

	  // Master list of positions and velocities
	  //double  Q1[N];
	  //double  Q2[N];
	  //double  Q3[N];
	  
	  ////double  U1[N];
	  //double  U2[N];
	  //double  U3[N];

	  // required since block was defined out of scope
	  //int block = int((N - remainder)/nthreads);
	  int flag;
	  MPI_Status recv_status;
	  MPI_Request recv_request;

	  /*
	  for (int k = 0; k<block; k++)
	    {
	      Q1[k] = temp[k];
	      
	    }
	  */
	  for (int i = 1; i<nthreads; i++)
	    {
	      
	      if(i == (nthreads - 1))
		{
		  n0 = block * i;
		  n1 = remainder + n0 + block;
		  expected = n1 - n0;
		}
	      else
		{
		  n0 = block * i;
		  n1 = n0 + (block);
		  expected = n1 - n0;
                  printf("Expecting %d elements for thread %d.\n", expected,i);
		}

	      //printf("Memory assigned.\n");

	      double tempQ1[expected];
	      double tempQ2[expected];
	      double tempQ3[expected];
	      //double * Q_buff[3 * expected];

	      double tempU1[expected];
	      double tempU2[expected];
	      double tempU3[expected];
	      //double * U_buff[3 * expected];
	      
	      //MPI_Status recv_status;
	      //MPI_Request recv_request;

	      
	      MPI_Recv(tempQ1, expected, MPI_DOUBLE, i, 1 + 6 * i, MPI_COMM_WORLD, &recv_status);
	      MPI_Recv(tempQ2, expected, MPI_DOUBLE, i, 2 + 6 * i, MPI_COMM_WORLD, &recv_status);
	      MPI_Recv(tempQ3, expected, MPI_DOUBLE, i, 3 + 6 * i, MPI_COMM_WORLD, &recv_status);
	      
	      MPI_Recv(tempU1, expected, MPI_DOUBLE, i, 4 + 6 * i, MPI_COMM_WORLD, &recv_status);
	      MPI_Recv(tempU2, expected, MPI_DOUBLE, i, 5 + 6 * i, MPI_COMM_WORLD, &recv_status);
	      MPI_Recv(tempU3, expected, MPI_DOUBLE, i, 6 + 6 * i, MPI_COMM_WORLD, &recv_status);
	      
	      /*
		MPI_Irecv(tempQ1, expected, MPI_DOUBLE, i, 1 + 6 * i, MPI_COMM_WORLD, &recv_request);
		MPI_Irecv(tempQ2, expected, MPI_DOUBLE, i, 2 + 6 * i, MPI_COMM_WORLD, &recv_request);
		MPI_Irecv(tempQ3, expected, MPI_DOUBLE, i, 3 + 6 * i, MPI_COMM_WORLD, &recv_request);
		MPI_Irecv(tempU1, expected, MPI_DOUBLE, i, 4 + 6 * i, MPI_COMM_WORLD, &recv_request);
		MPI_Irecv(tempU2, expected, MPI_DOUBLE, i, 5 + 6 * i, MPI_COMM_WORLD, &recv_request);
		MPI_Irecv(tempU3, expected, MPI_DOUBLE, i, 6 + 6 * i, MPI_COMM_WORLD, &recv_request);*/
	      
	      /*
		else 
		{
		for (int k = 0; i < block; i++)
		{
		Q1[k] = q1list[k];
		Q2[k] = q2list[k];
		Q3[k] = q3list[k];
		U1[k] = u1list[k];
		U2[k] = u2list[k];
		U3[k] = u3list[k];
		}*/
	      //printf("List Q1 = [%10.10f, %10.10f]\n",tempQ1[0],tempQ1[1]);
	      //printf("List Q2 = [%10.10f, %10.10f]\n",tempQ2[0],tempQ2[1]);
	      //printf("List Q3 = [%10.10f, %10.10f]\n",tempQ3[0],tempQ3[1]);
	      //printf("List U1 = [%10.10f, %10.10f]\n",tempU1[0],tempU1[1]);
	      //printf("List U2 = [%10.10f, %10.10f]\n",tempU2[0],tempU2[1]);
	      //printf("List U3 = [%10.10f, %10.10f]\n",tempU3[0],tempU3[1]);
	      
	      
	    

	      printf("Q1 = [%10.10f, %10.10f]\n",Q1[0],Q1[1]);
	      printf("Q2 = [%10.10f, %10.10f]\n",Q2[0],Q2[1]);
	      printf("Q3 = [%10.10f, %10.10f]\n",Q3[0],Q3[1]);
	      printf("U1 = [%10.10f, %10.10f]\n",U1[0],U1[1]);
	      printf("U2 = [%10.10f, %10.10f]\n",U2[0],U2[1]);
	      printf("U3 = [%10.10f, %10.10f]\n",U3[0],U3[1]);
	      
	      
	      
	      //printf("The size of the temporary arrays is %d and their elements have size %d.\n", );
	      
	      //	      MPI_Wait(&send_request, &send_status);
	      //	      MPI_Wait(&recv_request, &recv_status);


	      //double *U_buff = (double *)malloc(3 * expected * sizeof(double));
	      //double *Q_buff = (double *)malloc(3 * expected * sizeof(double));
	      printf("The block size is %d\n",block);
	     
	      for (int j = block * i; j < (block)*i + expected; j++)
		{
		  Q1[j] = tempQ1[j];
		  Q2[j] = tempQ2[j];
		  Q3[j] = tempQ3[j];
		  
		  U1[j] = tempU1[j];
		  U2[j] = tempU2[j];
		  U3[j] = tempU3[j];
		  printf("(i,j) = (%d,%d)\n",i,j);
		  
		  printf("Temp body %d has position [%10.10f,%10.10f, %10.10f], velocity [%10.10f,%10.10f,%10.10f]\n",j, tempQ1[j], tempQ2[j], tempQ3[j], tempU1[j], tempU2[j], tempU3[j]);
		}

	    }
	
	      
  
	      
		  // Should the memory for the temporary variables be freed?
		  //allocator::deallocate(*status,sizeof(status));
		  //allocator::deallocate(*temp,sizeof(temp));
      
	      //int flag;
	      //MPI_Test(&send_request,&flag, &send_status);
	      //MPI_Test(&recv_request,&flag, &recv_status);
      
    }
	  //MPI_Test(&send_request,&flag, &send_status);
	  //MPI_Test(&recv_request,&flag, &recv_status);

	  // Broadcast back to other threads

      printf("Preparing to broadcast %d value arrays from %d sized arrays...\n",N, (int)sizeof(Q1)/(int)sizeof(double));
      //printf("The position of body 2 is: [%10.10f, %10.10f, %10.10f]\n",Q1[1],Q2[1],Q3[1]);
      printf("Q1 = [%10.10f, %10.10f]\n",Q1[0],Q1[1]);
      printf("Q2 = [%10.10f, %10.10f]\n",Q2[0],Q2[1]);
      printf("Q3 = [%10.10f, %10.10f]\n",Q3[0],Q3[1]);
      printf("U1 = [%10.10f, %10.10f]\n",U1[0],U1[1]);
      printf("U2 = [%10.10f, %10.10f]\n",U2[0],U2[1]);
      printf("U3 = [%10.10f, %10.10f]\n",U3[0],U3[1]);
      
      
      //MPI_Scatter(Q1, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      //MPI_Scatter(Q2, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      //MPI_Scatter(Q3, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
      //MPI_Scatter(U1, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      //MPI_Scatter(U2, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      //MPI_Scatter(U3, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      


      MPI_Bcast(Q1,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(Q2,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(Q3,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(U1,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(U2,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(U3,N,MPI_DOUBLE,0,MPI_COMM_WORLD);


      /*
	double updatedQ1[N];
	double updatedQ2[N];
	double updatedQ3[N];
	
	double updatedU1[N];
	double updatedU2[N];
	double updatedU3[N];
	
	MPI_Scatter(Q1, N, MPI_DOUBLE, updatedQ1, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(Q2, N, MPI_DOUBLE, updatedQ2, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(Q3, N, MPI_DOUBLE, updatedQ3, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Scatter(U1, N, MPI_DOUBLE, updatedU1, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(U2, N, MPI_DOUBLE, updatedU2, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(U3, N, MPI_DOUBLE, updatedU3, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      */
      printf("Received Q1 = [%10.10f, %10.10f]\n",Q1[0],Q1[1]);
      printf("Received Q2 = [%10.10f, %10.10f]\n",Q2[0],Q2[1]);
      printf("Received Q3 = [%10.10f, %10.10f]\n",Q3[0],Q3[1]);
      printf("Received U1 = [%10.10f, %10.10f]\n",U1[0],U1[1]);
      printf("Received U2 = [%10.10f, %10.10f]\n",U2[0],U2[1]);
      printf("Received U3 = [%10.10f, %10.10f]\n",U3[0],U3[1]);
      
      for (int j = 0; j<N; j++)
	{
	  bodies[j].q1 = Q1[j];
	  bodies[j].q2 = Q2[j];
	  bodies[j].q3 = Q3[j];
	  
	  bodies[j].u1 = U1[j];
	  bodies[j].u2 = U2[j];
	  bodies[j].u3 = U3[j];
	  
	}
      
      // Synchronize
      //MPI_Barrier(MPI_COMM_WORLD);
  
      
      
      
    }  
  return 0;
    
}
