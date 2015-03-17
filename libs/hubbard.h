#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../libs/matrices.h"


#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>


#define t -1.0


typedef struct {
	int i;
	int N;
	void (*gfroutine2)(double _Complex , double _Complex **, void * );
	void *gfparams;
	double fermi;
	double imagEn;
	double hubU;
	double *m ;
	double *mnew;
	double *shift;
	double *n;
	double *zeeman;
}nroots_params2;


typedef struct {
	int N;
	double fermi;
	double imagEn;
	double init_m;
	double hubU;
	double msumdiff_thresh;
	double max_shift;
	double min_shift;
	int sym;
	double *zeeman;
}newHub_params;

	

typedef struct {
	int i;
	int N;
	double fermi;
	double imagEn;
	double _Complex **V;
	void (*gfroutine2)(double _Complex , double _Complex **, void * );
	void *gfparams;
	}intdos_params2;


double nroots_new (double shifti, void *p);
void newHub (double *m, double *n, double *shift, void *hub_params, void (*gfroutine2) (double _Complex, double _Complex **, void *), void *gfparams);
double intdos_new(double impart, void *p);



