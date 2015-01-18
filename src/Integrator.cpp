#include "../include/Body.h"
#include "../include/Forces.h"
#include "../include/Context.h"
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

  int N = bodies.size();
  vector<double> emptyForce(3,0.0);

  double h = evolveTime/(double)numberOfSteps;

  vector<vector<double> > forces(N,emptyForce);

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
      for (int i = 0; i<N; i++)
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
      
    }
  return 0;
}
