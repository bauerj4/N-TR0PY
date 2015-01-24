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
  cout << "THE PATH IS: " << path << '\n';
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
	  

	  vector<string> values;
	  string delim = " ";
	  delimitString(values, line, delim);
	 
	  Body.id = atoi(values[0].c_str());
	  //cout << "ID allocated. " << endl;
	  Body.type = atoi(values[1].c_str());
	  //cout << "Type allocated." << endl;
	  Body.q1 = atof(values[2].c_str());
	  //cout << "Body position 1 allocated. " << endl;
	  Body.q2 = atof(values[3].c_str());
	  //cout << "Body position 2 allocated. " << endl;
	  Body.q3 = atof(values[4].c_str());	  
	  //cout << "Body position 3 allocated. " << endl;
	  Body.u1 = atof(values[5].c_str());
	  //cout << "Body velocity 1 allocated. " << endl;	  
	  Body.u2 = atof(values[6].c_str());
	  //cout << "Body velocity 2 allocated. " << endl;
	  Body.u3 = atof(values[7].c_str());
	  //cout << "Body velocity 3 allocated. " << endl;
	 
	  Body.mass = atof(values[8].c_str());
	  //cout << "Body properties allocated." << endl;
	  bodies.push_back(Body);
	  //cout << "Body pushed." << endl;
	  //	  delete[] (&Body);
	  
	  //cout << "Iteration " << i << " completed." <<endl;
	  /* cout << "Position of body " << i << " is [" 
	    + values[2] + ", " + values[3] + ", " + values[4]
	    + "]" << endl;*/
	  i++;
	}
      //cout << bodies[0].type << endl << bodies[1].type << endl;
      file.close();
      return 0;
    }
  
  else
    {
      cout << "N-TR0PY ERROR: Read file " + path + " does not exist.\n";
      return 1;
    }
}

