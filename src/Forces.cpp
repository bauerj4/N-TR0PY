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

double PairForceCalculation(bodies_t &p, bodies_t &q, double eps2)
{  
  double dx = p.q1 - q.q1;
  double dy = p.q2 - q.q2;
  double dz = p.q3 - q.q3;

  double dr2 = pow(dx,2) + pow(dy,2) + pow(dz,2);

  double force =  (p.mass * q.mass) / (dr2 + eps2);
  return force;
}

int N2BruteForce(vector<bodies_t> &bodies, vector<vector<double> > &forces, double spaceVolume, double eps2)
{
  int N = bodies.size();
  
  double L = 0.5 * pow(spaceVolume, 1.0/3.0);
  domain_t masterDomain;

  masterDomain.x0 = -L;
  masterDomain.x1 = L;

  masterDomain.y0 = -L;
  masterDomain.y1 = L;
  
  masterDomain.z0 = -L;
  masterDomain.z1 = L;
 
  // check if all particles are in the domain.  This should be
  // an O(n) operation. Throw exception if a particle has left 
  // the simulation.



  // Compute forces. Since each operation has the same cost, there
  // is no need to consider loading factors between the threads.  

  MPI::Init();
  
  int nthreads, rank;
  
  MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (nthreads > N && rank ==0)
    {
      printf("There are %d processes and only %d bodies: runtime error expected.\n", nthreads, N);
    }
  
  // Equally split load
  // Minor efficiency loss for calculating these values on each thread
  
  int remainder = N % nthreads;

  int block = int((N - remainder)/nthreads);

  int n0,n1;

  // partitions the body list.
  if (rank == (nthreads - 1))
    {
      n0 = block * rank;
      n1 = remainder + block*(rank + 1);
    }

  else
    {
      n0 = block * rank;
      n1 = n0 + block * (rank + 1);
    }
  

  /*
    We now directly compute the forces.  
    It should be noted that this is a 
    stand-in for a double sum.  Since MPI
    cannot know if F_ij = -F_ji, we cannot
    use MPI to exploit the symmetry.  N(N-1)
    operations are performed, so MPI is not
    cost effective for fewer than two threads.
  */

  for (int i = n0; i < (n1 + 1),i++)
    {
      for (int j = 0; i < N; j++)
	{
	  if (i != j) // This is necessary to avoid computing "self" force.
	    {
	      vector<vector<double> > force = PairForceCalculation(bodies[i],bodies[j], 0.0); // temporarily set eps2 = 0     
	    }
	}
      
    }

  MPI::Finalize();
  
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