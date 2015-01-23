#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "mpi.h"

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

//The offset array is: {0,4,8,16,24,32,40,48,56}
/*MPI_Datatype MPI_BODY;
int  counts[] = {1,1,1,1,1,1,1,1};
int  offsets[] = {0,4,8,16,24,32,40,48,56};
MPI_Datatype types[] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

MPI_Type_struct(9, counts, offsets,types, &MPI_BODY);
MPI_Type_commmit(MPI_BODY);*/

/*
  For now, I intend to only support an ASCII output format,
  but I intend to implement a read method which will support 
  IFrIT's binary format in the future.
*/

MPI_Datatype createMPIBody();
int delimitString(vector<string>& tokens, string line, string delimiter);
int readASCII(vector<bodies_t>& bodies,string path);
