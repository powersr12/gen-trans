

#include "matrices.h"
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "graphene_gf.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

double I1int (double intvar, void  *p);
double I2int (double intvar, void  *p);
double I3int (double intvar, void  *p);
void HFsusMatrix(double _Complex **HFsus, double omega, double fermi, double eta, double *zeeman, int N, double *m, double *shift, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams);
void HFsusMatrix_comp(double _Complex **HFsus, double omega, double fermi, double eta, double *zeeman, int N, double *m, double *shift, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, double _Complex **I1, double _Complex **I2, double _Complex **I3 );

double _Complex HFsusElement(double omega, double fermi, double eta, double *zeeman, int N, double *m, double *shift, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, int i1, int j1, int k1, int l1);


double I1int_pristine (double intvar, void  *p);
double I2int_pristine (double intvar, void  *p);
double I3int_pristine (double intvar, void  *p);
void HFsusMatrix_pristine(double _Complex **HFsus, double omega, double fermi, double eta, double *zeeman, int N, double *m, double *shift, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams);

void HFsusMatrix_pristine_comp(double _Complex **HFsus, double omega, double fermi, double eta, double *zeeman, int N, double *m, double *shift, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, double _Complex **I1, double _Complex **I2, double _Complex **I3);

typedef struct	{
	double omega;
	double omegaf;
	double eta; 
	double *zeeman;
	int i1;
	int j1;
	int k1;
	int l1;
	int reim;
	int N;
	double *m;
	double *shift;
	void (*gfroutine)(double _Complex , double _Complex **, void *);
	void *gfparams;
	double hubU;
	}Iint_params;

