#include "../include/Body.h"
#include "../include/Forces.h"
#include "../include/Context.h"
#include <vector>

using namespace std;

//double computeTimestep();

vector<bodies_t> EulerMethod(vector<bodies_t> &bodies, context_t &NBODY_CONTEXT, int &snapshot_number);
