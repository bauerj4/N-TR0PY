#include "../include/Body.h"
#include "../include/Forces.h"
#include "../include/Context.h"
#include "../include/Integrator.h"
#include "mpi.h"
#include <vector>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  string path = "./test/TwoBody.dat";
  vector<bodies_t> bodies;
  readASCII(bodies, path);

  context_t NBODY_CONTEXT;
  NBODY_CONTEXT.eps2 = 0.0;
  NBODY_CONTEXT.init3Volume = 10.0;
  NBODY_CONTEXT.cosmology = "FLAT LAMBDA-CDM";

  //vector<vector<double> > forces;

  double evolveTime = 1.0; // Gy
  int NStep = 1000;
  
  EulerMethod(bodies, evolveTime, NStep, NBODY_CONTEXT); // This function should take an integration scheme specifier
  // maybe just pass &NBODY_CONTEXT?

  /*
    Get simulation context and pass it to the integrator.
    Choose the integrator.  
   */
  MPI_Finalize();
  return 0;
}
