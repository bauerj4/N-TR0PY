#include "../include/Body.h"
#include "../include/Forces.h"
#include <vector>
#include <string>

using namespace std;

int main(void)
{
  string path = "./test/TwoBody.dat";
  vector<bodies_t> bodies;
  readASCII(bodies, path);

  vector<vector<double> > forces;
  
  N2BruteForce(bodies, forces, 10.0, 0.0); // This function should take an integration scheme specifier

  /*
    Get simulation context and pass it to the integrator.
    Choose the integrator.  
   */
  
  return 0;
}
