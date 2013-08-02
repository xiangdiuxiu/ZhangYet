#include "Individual.h"

void Family::copy(const Family &rhs)
{
  include = rhs.include;
  pat = rhs.pat;
  mat = rhs.mat;
  kid.clear();
  for (unsigned int c=0; c<rhs.kid.size(); c++)
    kid.push_back(rhs.kid[c]);
}

