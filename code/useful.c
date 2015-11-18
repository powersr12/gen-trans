#include "useful.h"

int max(int num1, int num2)
{
  if(num1>= num2)
    return num1;
  else
    return num2;
}


double myRandNum (double min, double max)
{
  return min + (max-min)*( (double) rand() / (double) RAND_MAX);
}


//is a point (atom) inside a polygon (al antidot)?
int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) 
  {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
         (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}
