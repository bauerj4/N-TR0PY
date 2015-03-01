#include "../include/Body.h"
#include "../include/Domain.h"
#include <vector>
#include "mpi.h"
#include <cmath>

using namespace std;

/*
  This module contains force calculation methods.
*/

double G = 4.49209e-6;
// Compute the force between body p and q

vector<double> PairForceCalculation(bodies_t &p, bodies_t &q, double eps2, double &U)
{
  double G = 4.49209e-6; // kpc^3 / Gy^2 M_solar

  double dx = p.q1 - q.q1;
  double dy = p.q2 - q.q2;
  double dz = p.q3 - q.q3;

  //printf("Body p has position: [%10.10f, %10.10f, %10.10f]\n",p.q1,p.q2,p.q3);
  //printf("Body q has position: [%10.10f, %10.10f, %10.10f]\n",q.q1,q.q2,q.q3);


  double dr2 = pow(dx,2) + pow(dy,2) + pow(dz,2);
  //printf("The distance is %10.10f.\n", pow(dr2,0.5));

  double forceMag =  -  G * (q.mass) / pow((dr2 + eps2),1.5);
  U = -G * q.mass / pow(dr2,0.5);
  //double forceQuotient = -forceMag/pow(dr2,0.5);
  //printf("The force quotient is %10.10f\n", forceQuotient);
  //printf("The masses are %10.10f, %10.10f.\n",p.mass,q.mass);

  //double forceArr[] = {forceMag * dx, forceMag * dy, forceMag * dz};
  //vector<double> force(forceArr, forceArr + sizeof(forceArr));

  vector<double> force(3,0.0);
  force[0] = dx * forceMag;
  force[1] = dy * forceMag;
  force[2] = dz * forceMag;
  //printf("The force is %10.5f.\n", forceMag);
  //printf("Force Mag: %10.10f\n", forceMag);
  return force;
}

int N2BruteForce(vector<bodies_t> &bodies, vector<vector<double> > &forces, context_t &NBODY_CONTEXT, vector<double> &energies)
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
  
  double L = 0.5 * pow(NBODY_CONTEXT.INIT3VOLUME, 1.0/3.0);
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
      n1 = N ;
      n0 = 0;
    }
  
  else
    {
      block = int((N - remainder)/nthreads);
      
     
      // partitions the body list.
      if (rank == (nthreads - 1))
	{
	  n0 = block * rank;
	  n1 = remainder + block*(rank+1);
	}
      
      else 
	{
	  n0 = block * rank;
	  n1 = block * (rank + 1);
	}
      
    }

  /*
  printf("My rank is %d, and my search segment is [%d,%d]\n",rank,n0,n1);
  for (int i = n0; i < n1; i++)
    {
      printf("The particles %d have includes data:\n",n0);
      printf("r = [%5.10f, %5.10f, %5.10f]\n", bodies[i].q1,bodies[i].q2,bodies[i].q3);
      printf("u = [%5.10f, %5.10f, %5.10f]\n", bodies[i].u1,bodies[i].u2,bodies[i].u3);
      }*/
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
  double U_total = 0;
  //double keplerian_expected = 0;

  for (int i = n0; i < (n1);i++)
    {
      //printf("Computing force on Body %d...\n", bodies[i].id);
      vector<double> forceOnI(3,0.0);
      double U;
      //double T_total;

      for (int j = 0; j < N; j++)
	{
	  if (i != j) // This is necessary to avoid computing "self" force.
	    {
	      vector<double> force = PairForceCalculation(bodies[i],bodies[j], NBODY_CONTEXT.EPS2, U); // temporarily set eps2 = 0     	      
	      //printf("The force on %d is [%10.5f, %10.5f, %10.5f]\n",i, force[0], force[1], force[2]);
	      //printf("U = %10.10f\n", U);
	      forceOnI[0] += force[0];
	      forceOnI[1] += force[1];
	      forceOnI[2] += force[2];
	      U_total += U;
	    }	  
	}
      forces[i] = forceOnI;
      forceOnI[0] = 0.;
      forceOnI[1] = 0.;
      forceOnI[2] = 0.;
	
    }

  if (rank != 0 && NBODY_CONTEXT.SUPPRESS_DIAGNOSTICS==0)
    {
      MPI_Send(&U_total, 1, MPI_DOUBLE, 0, rank,  MPI_COMM_WORLD);
    }

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 0 && NBODY_CONTEXT.SUPPRESS_DIAGNOSTICS==0)
    {
      double U_temp;
      double keplerian;
      
      for (int j = 1; j < nthreads; j++)
	{
	  MPI_Status status;
	  MPI_Recv(&U_temp, 1, MPI_DOUBLE, j, j, MPI_COMM_WORLD,&status);
	  U_total += U_temp;
	  U_temp = 0.;
	}
      double T_total = 0;
      for (int i = 0; i < N; i++)
	{
	  T_total += bodies[i].u1*bodies[i].u1 + bodies[i].u2 * bodies[i].u2 + bodies[i].u3 * bodies[i].u3;
	}

      double dx = bodies[0].q1 - bodies[1].q1;
      double dy = bodies[0].q2 - bodies[1].q2;
      double dz = bodies[0].q3 - bodies[1].q3;

      double dr = pow(dx*dx + dy*dy + dz*dz,0.5);

      keplerian = -G * bodies[0].mass / (2. * 0.0010000) + 0.000007; // Initial energy guess is off in the 7th decimal place
      T_total /= 2.00000000000;
      T_total *= 1.022/*69032072*/* 1.022;//69032072;
      U_total /= 2.00000000000;

      printf("(T,U,E) = (%10.10f, %10.10f, %10.10f) -> %10.10f\n", T_total, U_total, T_total+U_total, keplerian - (T_total + U_total));
      double energy = log10(fabs((U_total + T_total - keplerian)/keplerian));

      printf("PLOT VALUE IS: %10.10f\n",energy);
      energies.push_back(energy);
      //printf("The total energy is %10.10f\n", energy);

      T_total = 0;
      //U_total = 0;
      energy = 0;
    }
  U_total = 0;
  MPI_Barrier(MPI_COMM_WORLD);

  //forceOnI[0] = 0.;
  //forceOnI[1] = 0.;
  // forceOnI[2] = 0.;
                                                    
  
      
}
// Synchronize before integration.
  

  /*
    We shouldn't do the integration here.  If a method requires multiple calls
    to the force calculation, having segments of force operations is messy.  Integration
    should be called separately, and this function should be called from the
    integrator.  There are potential issues with not doing this, namely that we will
    mix up the order of the forces.  Each MPI thread has its OWN copy of bodies[] and
    will thus need to be synchronized.  
  */  



/*
  Will compute O(n log(n)) time forces for the ODE scheme. 
  Algorithm extracted from Barnes and Hut 1986.
*/


/*
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
*/
