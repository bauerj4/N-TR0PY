#include <iostream>
#include <string.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "../include/Body.h"


using namespace std;

/*

  This class includes methods for defining, reading, and writing
  information about simulation bodies.

*/


/*
  For now, I intend to only support an ASCII output format,
  but I intend to implement a read method which will support 
  IFrIT's binary format in the future.
*/

MPI_Datatype createMPIBody()
{

  MPI_Datatype  MPI_BODY;
  int  counts[] = {1,1,1,1,1,1,1,1,1}; 
  long int  offsets[] = {0,4,8,16,24,32,40,48,56};
  MPI_Datatype types[] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  
  MPI_Type_struct(9, counts, offsets,types, &MPI_BODY);
  MPI_Type_commit(&MPI_BODY);
  int size;
  MPI_Type_size(MPI_BODY,&size);
  //printf("MPI_BODY IS DECLARED WITH SIZE: %d\n",(int) size);
  return MPI_BODY;
}

int delimitString(vector<string> &tokens, string line, string delimiter)
{
  int pos = 0;
  int newpos = 0;
  string token;
  while(1)
    {
      newpos = line.find(' ');
      token = line.substr(0,newpos);
      line.erase(0,newpos+1);
      //cout << "(pos, newpos) is (" + pos + ',' + newpos +')' << endl; 
      //cout << "Token is " + token << endl;
      
      tokens.push_back(token);
      if (newpos == string::npos)
	{
	  break;
	}
    }
}

int readASCII(vector<bodies_t>& bodies,string path)
{
  ifstream file;
  //cout << "THE PATH IS: " << path << '\n';
  file.open(path.c_str(), ifstream::in);
  string line;
  
  if (file.is_open())
    {
      int i=0;
      while(getline(file,line))
	{
	  //int type,id;
	  bodies_t Body;
	  //double body_q1, body_q2, body_q3;
	  //double body_u1, body_u2, body_u3;
	  
	  //cout << line << "\n";
	  vector<string> values;
	  string delim = " ";
	  delimitString(values, line, delim);
	 
	  Body.id = atoi(values[0].c_str());
	  
	  Body.type = atoi(values[1].c_str());
	
	  Body.q1 = atof(values[2].c_str());
	  
	  Body.q2 = atof(values[3].c_str());
	 
	  Body.q3 = atof(values[4].c_str());	  
      
	  Body.u1 = atof(values[5].c_str());
	  
	  Body.u2 = atof(values[6].c_str());
  
	  Body.u3 = atof(values[7].c_str());
	 
	  Body.mass = atof(values[8].c_str());
	  
	  bodies.push_back(Body);
     
	  i++;
	  //printf("Iteration %d in reading completeted.\n",i);
	}
   
      file.close();
      return 0;
    }
  
  else
    {
      cout << "N-TR0PY ERROR: Read file " + path + " does not exist.\n";
      return 1;
    }
}

int writeSnapshot(vector<bodies_t> &bodies, context_t &context, int snapshot_number)
{
  ofstream output;
  stringstream ss;
  string PATH;


  ss << context.OUTPUT_PATH << "_" << snapshot_number;
  PATH = ss.str();
  output.open(PATH.c_str());

  for (int i = 0; i < bodies.size(); i++)
    {
      output << bodies[i].id << " ";
      output << bodies[i].type << " ";
      output << bodies[i].mass << " ";
      output << bodies[i].q1 << " ";
      output << bodies[i].q2 << " ";
      output << bodies[i].q3 << " ";
      output << bodies[i].u1 << " ";
      output << bodies[i].u2 << " ";
      output << bodies[i].u3 << "\n";
    }
  return 0;
}
