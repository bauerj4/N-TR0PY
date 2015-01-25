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

  char * ctx_path = argv[1];
  //char * ctx_path = (char*)ctx_path_string.c_str();

  int current_snapshot = 0;
  
  context_t NBODY_CONTEXT;

  readContext(ctx_path, NBODY_CONTEXT);
  printContext(NBODY_CONTEXT);

  vector<bodies_t> bodies;
  readASCII(bodies, NBODY_CONTEXT.BODIES_FILE_PATH);
  //printf("Bodies loaded.\n");

  writeSnapshot(bodies, NBODY_CONTEXT, current_snapshot);
  current_snapshot += 1;
  
  bodies = EulerMethod(bodies, NBODY_CONTEXT, current_snapshot); // This function should take an integration scheme specifier
  // maybe just pass &NBODY_CONTEXT?

  /*
    Get simulation context and pass it to the integrator.
    Choose the integrator.  
   */

  int rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    {
      printf("Writing output...\n");
      writeSnapshot(bodies, NBODY_CONTEXT, current_snapshot);
    }
  
  MPI_Finalize();
  return 0;
}
