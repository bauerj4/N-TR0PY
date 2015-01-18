                                                                                                         
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

#pragma once
struct context_t
{
  double eps2;
  double init3Volume;
  string cosmology; // In the future, non-flat geometries and non-Lambda CDM will be supported.
};
