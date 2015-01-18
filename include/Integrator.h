#include "../include/Body.h"
#include "../include/Forces.h"
#include "../include/Context.h"
#include <vector>

using namespace std;

//double computeTimestep();

int EulerMethod(vector<bodies_t> &bodies, double evolveTime, int numberOfSteps, context_t &NBODY_CONTEXT);
