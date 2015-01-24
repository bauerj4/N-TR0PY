#include "../include/Context.h"
#include "../include/Body.h"
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

/*
  We should define methods for reading context from a file that should
  be given as an argument to the simulator.  Define similar methods as 
  the ones defined in Body.cpp, but we need not use a binary format since
  the number of entries SHOULD be small.
*/

using namespace std;

int readContext(char * PATH, context_t &context)
{
  ifstream file; 
  file.open(PATH,ifstream::in);
  vector<string> tokens;

  if (file.is_open())
    {
      int i=0;
      string line;
      while(getline(file,line))
        {
	  delimitString(tokens,line, " = ");
	  i++;
	}

      for (int j = 0; j < tokens.size(); j++)
	{
	  if(tokens[j] == "SUPPRESS_DIAGNOSTICS")
            {
              context.SNAPSHOTS_EVERY_N =  atoi(tokens[j + 2].c_str());
            }

          if(tokens[j] == "SNAPSHOTS_EVERY_N")
            {
              context.SNAPSHOTS_EVERY_N =  atoi(tokens[j + 2].c_str());
            }

          if(tokens[j] == "NSTEP")
            {
	      context.NSTEP =  atoi(tokens[j + 2].c_str());
            }


	  if(tokens[j] == "BARNES_OPENING_ANGLE")
	    {
	      context.BARNES_OPENING_ANGLE = atof(tokens[j + 2].c_str());
	    }
         
	  if(tokens[j] == "SIMULATION_TIME")
            {
              context.SIMULATION_TIME = atof(tokens[j + 2].c_str());
            }

	  if (tokens[j] == "EPS2")
	    {
	      context.EPS2 = atof(tokens[j + 2].c_str());
	    }

          if (tokens[j] == "DOMAIN_VOLUME")
            {
              context.INIT3VOLUME = atof(tokens[j + 2].c_str());
            }

          if(tokens[j] == "FORCE_CALCULATOR")
            {
              context.FORCE_CALCULATOR = tokens[j + 2];
            }


          if(tokens[j] == "COSMOLOGY")
            {
	      context.COSMOLOGY= tokens[j + 2];
            }

          if(tokens[j] == "GEOMETRY")
            {
              context.GEOMETRY = tokens[j + 2];
            }

          if(tokens[j] == "DOMAIN_TYPE")
            {
              context.DOMAIN_TYPE = tokens[j + 2];
            }

          if(tokens[j] == "EXTERNAL_POTENTIAL")
            {
              context.EXTERNAL_POTENTIAL = tokens[j + 2];
            }

          if(tokens[j] == "BODIES_FILE_PATH")
            {
              context.BODIES_FILE_PATH = tokens[j + 2];
            }

          if(tokens[j] == "SNAPSHOT_PATH")
            {
              context.OUTPUT_PATH = tokens[j + 2];
            }

	  //cout << "Read token is " + tokens[j] << endl;
	}
    }
  file.close();

  return 0;
}


int printContext(context_t &context)
{
  printf("THE CONTEXT IS: \n");
  printf("SUPPRESS_DIAGNOSTICS = %d\n",context.SUPPRESS_DIAGNOSTICS);
  printf("SNAPSHOTS_EVERY_N = %d\n", context.SNAPSHOTS_EVERY_N);
  printf("NSTEP = %d\n", context.NSTEP);
  printf("EPS2 = %10.10f\n", context.EPS2);
  printf("DOMAIN_VOLUME = %10.10f\n", context.INIT3VOLUME);
  printf("BARNES_OPENING_ANGLE = %10.10f\n", context.BARNES_OPENING_ANGLE);
  printf("SIMULATION_TIME = %10.10f\n", context.SIMULATION_TIME);
  cout << "FORCE_CALCULATOR = "  << context.FORCE_CALCULATOR  << '\n';
  cout << "COSMOLOGY = " << context.COSMOLOGY << '\n';
  cout << "GEOMETRY = " << context.GEOMETRY << '\n';
  cout << "EXTERNAL_POTENTIAL = " << context.EXTERNAL_POTENTIAL << '\n';
  cout << "BODIES_FILE_PATH = " << context.BODIES_FILE_PATH << '\n';
  cout << "SNAPSHOT_PATH = " << context.OUTPUT_PATH << '\n';
  cout << "DOMAIN_TYPE = " << context.DOMAIN_TYPE << '\n';
}
//delimitString(vector<string> &tokens, string line, string delimiter);
