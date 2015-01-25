#include <vector>
#include "../include/Body.h"
#include "../include/Context.h"


vector<double> PairForceCalculation(bodies_t &p, bodies_t &q, double eps2);

int N2BruteForce(vector<bodies_t> &bodies, vector<vector<double> > &forces, context_t &NBODY_CONTEXT);

//int BarnesHut(vector<bodies_t> &bodies, vector<vector<double> > &forces, context_t &NBODY_CONTEXT);
