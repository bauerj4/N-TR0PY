#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "mpi.h"
#include "Context.h"

using namespace std;

#pragma once
struct bodies_t  
{
  int id;
  int type;
  
  double mass;
  
  double q1;
  double q2;
  double q3;
  
  double u1;
  double u2;
  double u3;
};

/*
  For now, I intend to only support an ASCII output format,
  but I intend to implement a read method which will support 
  IFrIT's binary format in the future.
*/

MPI_Datatype createMPIBody();
int delimitString(vector<string>& tokens, string line, string delimiter);
int readASCII(vector<bodies_t>& bodies,string path);
int writeSnapshot(vector<bodies_t> &bodies, context_t &context, vector<double> &times, vector<double> &e, int &num);
