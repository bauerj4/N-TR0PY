#include <vector>
#include "../include/Body.h"

vector<double> PairForceCalculation(bodies_t &p, bodies_t &q, double eps2);

int N2BruteForce(vector<bodies_t> &bodies, vector<vector<double> > &forces,double spaceVolume,  double eps2);

int BarnesHut(vector<bodies_t> &bodies, vector<vector<double> > &forces, double spaceVolume, double eps2);
