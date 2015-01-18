#include "../include/Body.h"
#include "../include/Domain.h"
#include <vector>
#include "mpi.h"
#include <cmath>

using namespace std;

/*
  This module contains force calculation methods.
*/


// Compute the force between body p and q

vector<double> PairForceCalculation(bodies_t &p, bodies_t &q, double eps2)
{  
  double dx = p.q1 - q.q1;
  double dy = p.q2 - q.q2;
  double dz = p.q3 - q.q3;

  double dr2 = pow(dx,2) + pow(dy,2) + pow(dz,2);

  double forceMag =  (p.mass * q.mass) / (dr2 + eps2);
  double forceQuotient = forceMag/dr2;
  
  double forceArr[] = {forceQuotient * dx, forceQuotient * dy, forceQuotient * dz};
  vector<double> force(forceArr, forceArr + sizeof(forceArr));
  return force;
}

int N2BruteForce(vector<bodies_t> &bodies, vector<vector<double> > &forces, double spaceVolume, double eps2)
{ 
  // check if all particles are in the domain.  This should be
  // an O(n) operation. Throw exception if a particle has left 
  // the simulation.


  // Compute forces. Since each operation has the same cost, there
  // is no need to consider loading factors between the threads.  
  
  int nthreads, rank;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int N = bodies.size();
  
  double L = 0.5 * pow(spaceVolume, 1.0/3.0);
  domain_t masterDomain;
  
  masterDomain.x0 = -L;
  masterDomain.x1 = L;
  
  masterDomain.y0 = -L;
  masterDomain.y1 = L;
  
  masterDomain.z0 = -L;
  masterDomain.z1 = L;
    

  if (nthreads > N && rank == 0)
    {
      printf("There are %d processes and only %d bodies: runtime error expected.\n", nthreads, N);
    }
  
  // Equally split load
  // Minor efficiency loss for calculating these values on each thread
  
  int remainder = N % nthreads;
  int block,n0,n1;

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

  
  printf("My rank is %d, and my search segment is [%d,%d]\n",rank,n0,n1);
  /*
    We now directly compute the forces.  
    It should be noted that this is a 
    stand-in for a double sum.  Since MPI
    cannot know if F_ij = -F_ji, we cannot
    use MPI to exploit the symmetry.  N(N-1)
    operations are performed, so MPI is not
    cost effective for fewer than two threads.
  */
  
  MPI_Barrier(MPI_COMM_WORLD);

  for (int i = n0; i < (n1 + 1);i++)
    {
      //printf("Computing force on Body %d...\n", bodies[i].id);
      vector<double> forceOnI(3,0.0);
      for (int j = 0; j < N; j++)
	{
	  if (i != j) // This is necessary to avoid computing "self" force.
	    {
	      vector<double> force = PairForceCalculation(bodies[i],bodies[j], 0.0); // temporarily set eps2 = 0     	      
	      //printf("The force is [%10.5f, %10.5f, %10.5f]\n", force[0], force[1], force[2]);
	      forceOnI[0] += force[0];
	      forceOnI[1] += force[1];
	      forceOnI[2] += force[2];
	    }	  
	  forces[i] = forceOnI;
	}
      
    }

  MPI_Barrier(MPI_COMM_WORLD); // Synchronize before integration.
  

  /*
    We shouldn't do the integration here.  If a method requires multiple calls
    to the force calculation, having segments of force operations is messy.  Integration
    should be called separately, and this function should be called from the
    integrator.  There are potential issues with not doing this, namely that we will
    mix up the order of the forces.  Each MPI thread has its OWN copy of bodies[] and
    will thus need to be synchronized.  
  */  
}


/*
  Will compute O(n log(n)) time forces for the ODE scheme. 
  Algorithm extracted from Barnes and Hut 1986.
*/


int BarnesHut(vector<bodies_t> &bodies, vector<vector<double> > &forces,double spaceVolume, double eps2)
{
  int N = bodies.size(); // number of bodies

  // build the master node, check to make sure that all bodies are contained within the domain. 
  // Exceptions should be thrown when bodies are outside the simulation domain.

  // Assume a cube
  
  double L = 0.5 * pow(spaceVolume, 1.0/3.0); 
  domain_t masterDomain;
  
  masterDomain.x0 = -L;
  masterDomain.x1 = L;
 
  masterDomain.y0 = -L;
  masterDomain.y1 = L;
  
  masterDomain.z0 = -L;
  masterDomain.z1 = L;

  // Parallel tree construction is tricky.  


  return 0;
}
