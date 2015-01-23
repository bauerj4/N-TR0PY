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
  /*
  long int body_address, id_address, type_address, mass_address, u1_address, u2_address, u3_address, q1_address, q2_address, q3_address;

  int id_offset, type_offset, mass_offset, u1_offset, u2_offset, u3_offset, q1_offset, q2_offset, q3_offset;

  bodies_t test_body = bodies[0];

  MPI_Get_address(&test_body, &body_address);
  MPI_Get_address(&test_body.id, &id_address);
  MPI_Get_address(&test_body.type, &type_address);
  MPI_Get_address(&test_body.q1, &q1_address);
  MPI_Get_address(&test_body.q2, &q2_address);
  MPI_Get_address(&test_body.q3, &q3_address);
  MPI_Get_address(&test_body.u1, &u1_address);
  MPI_Get_address(&test_body.u2, &u2_address);
  MPI_Get_address(&test_body.u3, &u3_address);
  MPI_Get_address(&test_body.mass, &mass_address);

  
  id_offset = (int)(id_address - body_address);
  type_offset = (int)(type_address - body_address);
  q1_offset = (int)(q1_address - body_address);
  q2_offset = (int)(q2_address - body_address);
  q3_offset = (int)(q3_address - body_address);
  u1_offset = (int)(u1_address - body_address);
  u2_offset = (int)(u2_address - body_address);
  u3_offset = (int)(u3_address - body_address);
  mass_offset = (int)(mass_address - body_address);

  printf("The offset array is: {%d,%d,%d,%d,%d,%d,%d,%d,%d}\n", id_offset, type_offset, mass_offset, q1_offset, q2_offset,
	 q3_offset, u1_offset, u2_offset, u3_offset);
  */


  double evolveTime = 1.0; // Gy
  int NStep = 500;
  
  bodies = EulerMethod(bodies, evolveTime, NStep, NBODY_CONTEXT); // This function should take an integration scheme specifier
  // maybe just pass &NBODY_CONTEXT?

  /*
    Get simulation context and pass it to the integrator.
    Choose the integrator.  
   */
  MPI_Finalize();
  return 0;
}
