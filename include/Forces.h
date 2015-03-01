#include <vector>
#include "../include/Body.h"
#include "../include/Context.h"


vector<double> PairForceCalculation(bodies_t &p, bodies_t &q, double eps2, double &U);

int N2BruteForce(vector<bodies_t> &bodies, vector<vector<double> > &forces, context_t &NBODY_CONTEXT, vector<double> &energies);

//int BarnesHut(vector<bodies_t> &bodies, vector<vector<double> > &forces, context_t &NBODY_CONTEXT);
