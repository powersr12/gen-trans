#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>



double graphenegf(double _Complex En, int a1, int a2, int reim, int diagt);
double intg(double kx, void *p);
void getS1S2 (double _Complex *S1, double _Complex *S2, double x, double y, double kx, double _Complex cosky1, double _Complex cosky2, double _Complex En);
void graph_NNS (double _Complex **g, double _Complex En, double a1, double a2, int diagt);
int getdiagt (int d1, int d2);
void MagGreen(double _Complex upenergy, double _Complex downenergy, double _Complex **Gup, double _Complex **Gdown, double *m, double *shift, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, double hubU, int N, double *zeemanterm);
	double graphene_SPA_AC (double _Complex En, int a1, int a2, int reim, int diagt);
	void grapheneSPAmatrix( double _Complex En, double _Complex **G, void *p);
void grapheneGFmatrix( double _Complex En, double _Complex **G, void *p);
void graphenebridgematrix( double _Complex En, double _Complex **G, void *p);
void graphenetopmatrix( double _Complex En, double _Complex **G, void *p);
void graphenecentrematrix( double _Complex En, double _Complex **G, void *p);
double _Complex graphenecentresum ( double _Complex En, double _Complex **G, void *p);

void NonMagGreen(double _Complex upenergy, double _Complex downenergy, double _Complex **Gup, double _Complex **Gdown, double *m, double *shift, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, double hubU, int N, double *zeemanterm);
double getlocalm (double *m, double *shift, double *zeemanterm, double hubU, void (*gfroutine2) (double _Complex, double _Complex **, void *), void *gfparams, int i, int numsites, double fermi, double imagEn);
double intdos_other(double impart, void *p);

double graphenegf_strain(double _Complex En, int a1, int a2, int reim, int diagt, double t1, double t2);
double intg_strain(double kx, void *p);
void getS1S2_strain (double _Complex *S1, double _Complex *S2, double x, double y, double kx, double _Complex cosky1, double _Complex cosky2, double _Complex En, double t1, double t2);
void get_t1t2 (double epsilon, double sigma, int direction, double *t1, double *t2);
	void grapheneGFmatrix_strain( double _Complex En, double _Complex **G, void *p);



typedef struct {
	int i;
	int N;
	double fermi;
	double imagEn;
	double _Complex **V;
	void (*gfroutine2)(double _Complex , double _Complex **, void * );
	void *gfparams;
	}intdosother_params;

typedef struct {
	double x;
	double y;
	double _Complex En;
	int reim;
	int diagt; // 0 for diagonal, 1 for 12-type off diagonal, 2 for 21-type offdiagonal
}intg_params;

typedef struct {
	int x;
	int y;
	double _Complex En;
	int reim;
	int diagt; // 0 for diagonal, 1 for 12-type off diagonal, 2 for 21-type offdiagonal
	double t1;
	double t2;
}intgstrain_params;

	

typedef struct {
	int **pos;
	int N;
	int psym;
}gfmat_p ;
		

typedef struct {
	int **pos;
	int N;
	double tprime;
}gfmatb_p ;
			


typedef struct {
	int **pos;
	int N;
	int psym;
	double t1;
	double t2;
}gfmats_p ;
			

//tight binding parameters
#define dis 1.0
#define t -1.0
