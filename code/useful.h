#pragma once
#include <stdlib.h>
#include "../libs/matrices.h"


int max(int num1, int num2);
double myRandNum (double min, double max);
int mymin(int num1, int num2);
int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy);
int pnoronpoly(int nvert, double *vertx, double *verty, double testx, double testy);

void swapxy (double **posarray, int dim);
int myRandInt (int min, int max);
double dmax(double num1, double num2);
double dmin(double num1, double num2);
double areaTriangle(double *vertx, double *verty);

double twodimCP(double *v1, double *v2);
int SameSide (double *p1, double *p2, double *a, double *b);
int pntriangle(double *vertx, double *verty, double testx, double testy);
