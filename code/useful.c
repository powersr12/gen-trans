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
