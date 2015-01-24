                                                                                                         
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

#pragma once
struct context_t
{
  int SUPPRESS_DIAGNOSTICS;
  int SNAPSHOTS_EVERY_N;
  int NSTEP;
  double EPS2;
  double INIT3VOLUME;
  double BARNES_OPENING_ANGLE;
  double SIMULATION_TIME;
  string FORCE_CALCULATOR;
  string EXTERNAL_POTENTIAL; // Each will have its own method of taking params
  string COSMOLOGY; // In the future, non-flat geometries and non-Lambda CDM will be supported.
  string GEOMETRY;
  string BODIES_FILE_PATH;
  string OUTPUT_PATH;
  string DOMAIN_TYPE;
};

int readContext(char * PATH, context_t &context);
int printContext(context_t &context);
