
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include "matrices.h"

typedef struct{
	double realenergy;
	double *m ;
	double *shift;
	double *zeeman;
	double hubU;
	void (*gfroutine)(double _Complex , double _Complex **, void * );
	void *gfparams;
}intcoupling_params;

typedef struct{
	double realenergy;
	double *m ;
	double *shift;
	double hubU;
	int Ntotal;
	int Nblack;
	int Nwhite;
	void (*gfroutine)(double _Complex , double _Complex **, void * );
	void *gfparams;
}intEnsemble_params;

	double GetCoupling(double *m, double *shift, double *zeeman, double _Complex En, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams);
	double intcoupling (double impart, void *p);
	double GetEnsembleDiff(double *m, double *shift, double _Complex En, double hubU, int Ntotal, int Nblack, int Nwhite, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams);
	double intEnsembleDiff (double impart, void *p);
