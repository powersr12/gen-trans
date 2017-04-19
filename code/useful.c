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

double dmax(double num1, double num2)
{
  if(num1>= num2)
    return num1;
  else
    return num2;
}


double dmin(double num1, double num2)
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

//is a point (atom) inside (or on the edge of) a polygon ?
int pnoronpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) 
  {
    if ( ((verty[i]>=testy) != (verty[j]>=testy)) &&
         (testx <= (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}


//two-dim cross product (returns z direct scalar)
double twodimCP(double *v1, double *v2)
{
  return v1[0]*v2[1] - v2[0]*v1[1];
}

//are two points on the same side of a line?
int SameSide (double *p1, double *p2, double *a, double *b)
{
	int i;
	double a1[2], a2[2], a3[2];
	
	for(i=0; i<2; i++)
	{
		a1[i] = b[i]-a[i];
		a2[i] = p1[i] -a[i];
		a3[i] = p2[i] -a[i];
	}
	
	if (twodimCP(a1, a2) * twodimCP(a1, a3) >=0)
		return 1;
	else return 0;
}

//is a point (atom) inside or on a triangle?
int pntriangle(double *vertx, double *verty, double testx, double testy)
{
	double testpoint[2] ={testx, testy};
	double pointA[2] = {vertx[0], verty[0]};
	double pointB[2] = {vertx[1], verty[1]};
	double pointC[2] = {vertx[2], verty[2]};
	
	if(SameSide(testpoint, pointA, pointB, pointC) && SameSide(testpoint, pointB, pointA, pointC) && SameSide(testpoint, pointC, pointA, pointB))
		return 1;
	else return 0;
}


// function PointInTriangle(p, a,b,c)
//     if SameSide(p,a, b,c) and SameSide(p,b, a,c)
//         and SameSide(p,c, a,b) then return true
//     else return false

double areaTriangle(double *vertx, double *verty)
{
  return 0.5*fabs( (vertx[0] -vertx[2])*(verty[1]-verty[0]) - (vertx[0] -vertx[1])*(verty[2] -verty[0])) ;
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
