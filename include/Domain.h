#include <vector>
#include "Body.h"

/*
  Used for constructing simulation domain as well as subdomains
  for tree codes.
*/

using namespace std;


typedef struct domain_t
{
  double x0;
  double y0;
  double z0;

  double x1;
  double y1;
  double z1;

  domain_t* isInDomain;
  domain_t* containsDomains[];
  
  vector<bodies_t> /* should be pointer? */ bodiesInDomain;

} Domain;

double computeVolume(domain_t &domain);
