#include "useful.h"

int max(int num1, int num2)
{
  if(num1>= num2)
    return num1;
  else
    return num2;
}


int mymin(int num1, int num2)
{
  if(num1<= num2)
    return num1;
  else
    return num2;
}
double myRandNum (double min, double max)
{
  return min + (max-min)*( (double) rand() / (double) RAND_MAX);
}

//generates a random INTEGER between min and max INCLUSIVE
int myRandInt (int min, int max)
{
  return min + ((rand()) % (max-min+1)) ;
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

void swapxy (double **posarray, int dim)
{
  int i, j;
  double **temp = createNonSquareDoubleMatrix(dim, 2);
  
  for(i=0; i<dim; i++)
  {
    temp[i][0] = posarray[i][1];
    temp[i][1] = posarray[i][0];
    
    posarray[i][0] = temp[i][0];
    posarray[i][1] = temp[i][1];
  }
  free(temp[0]); free(temp);

  
  
}
