
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include "../libs/matrices.h"

typedef struct {
	int size;
	int numsites;
	int *pos;
	void (*gfroutine)(double _Complex **, double _Complex **, double _Complex **, int , double, double _Complex );
	double rubio_error;
	double hopping;
	}GFsites_params;
	
	typedef struct {
	int size;
	int numsites;
	int *pos;
	void (*gfroutine)(double _Complex **, double _Complex **, double _Complex **, double , double, double _Complex );
	double rubio_error;
	double hopping;
	double kpar;
	}GFsites2_params;

	typedef struct {
	int size;
	int numsites;
	int *pos;
	void (*gfroutine)(double _Complex **, double _Complex **, double _Complex **, double , double, double _Complex );
	double rubio_error;
	double hopping;
	double _Complex En;
	int cell1_edge_dist;
	int cell2_edge_dist;
	int cell_para_dist;
	int cell1site ;
	int cell2site ; 
	int reim;
	}GFsites2int_params;

typedef struct {
	void (*gfroutine)(double _Complex **, double _Complex **, double _Complex **, int , double, double _Complex );
	int size;
	double hopping;
	double rubio_error;
	int numsites;
	int *pos;
	int blockfunc;
	int blocksize;
	int block2size;
	double _Complex *DA;
	int numij;
	int *posij;
	}RibbonGF_params;
	
	typedef struct {
	void (*gfroutine)(double _Complex **, double _Complex **, double _Complex **, int , double, double _Complex );
	int size;
	double hopping;
	double rubio_error;
	int numcells;
	int numsites;
	int *pos;
	double _Complex *DA;
	}RibbonGFs_params;

typedef struct {
	void (*gfroutine)(double _Complex **, double _Complex **, double _Complex **, int , double, double _Complex );
	int size;
	double hopping;
	double rubio_error;
	int numsites1;
	int numcells1;
	int *pos1;
	int numsites2;
	int numcells2;
	int *pos2;
	int blockfunc;
	int blocksize;
	int block2size;
	double _Complex *DA;
	int ads;
	double ads_hopping;
	}RibbonGF_coupparams;


typedef struct {
	int i;
	int N;
	double fermi;
	double imagEn;
	double _Complex **V;
	void (*gfroutine2)(double _Complex , double _Complex **, void * );
	void *gfparams;
	}intdosother_params;


// 	void MagGreen(double _Complex upenergy, double _Complex downenergy, double _Complex **Gup, double _Complex **Gdown, double *m, double *shift, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, double hubU, int N, double *zeemanterm);
// 
// 	double getlocalm (double *m, double *shift, double *zeemanterm, double hubU, void (*gfroutine2) (double _Complex, double _Complex **, void *), void *gfparams, int i, int numsites, double fermi, double imagEn);
// 
// 	double intdos_other(double impart, void *p);

	void ACinf (double _Complex **g11inv, double _Complex **V01, double _Complex **V10, double kpar, double hopping, double _Complex Eng);

	
	void RubioSGF (double _Complex **S, double _Complex **g00, double _Complex **V01, double _Complex **V10, int N, int *count, double error);

	void ZGNR (double _Complex **g11inv, double _Complex **V01, double _Complex **V10, int N, double hopping, double _Complex Eng);

	void AGNR (double _Complex **g11inv, double _Complex **V01, double _Complex **V10, int N, double hopping, double _Complex Eng);

	void AGNR_E (double _Complex **g11inv, double _Complex **V01, double _Complex **V10, int N, double hopping, double _Complex Eng);

	void GFsites (double _Complex En, double _Complex **G, void *p);
	
	void GFsites2 (double _Complex En, double _Complex **G, void *p);

	
	double GFsites2int (double kpar, void *p);


	void SimpleRubioBlock (double _Complex **G11, double _Complex **GLL, double _Complex **G1L, double _Complex **GL1, double _Complex **g00, double _Complex **V01, double _Complex **V10, int numits, int dim);

	void SimpleBlock (double _Complex **G11, double _Complex **GLL, double _Complex **G1L, double _Complex **GL1, double _Complex **g00, double _Complex **V01, double _Complex **V10, int numits, int dim, double _Complex *DisorderArray);
	
	void VSimpleBlock (double _Complex **G11, double _Complex **GLL, double _Complex **G1L, double _Complex **GL1, double _Complex **g00, double _Complex **V01, double _Complex **V10, int numits, int dim);



	void ConnectBlock (double _Complex **g00, double _Complex **g11, double _Complex **gLL, double _Complex **g1L, double _Complex **gL1, double _Complex **newSF, double _Complex **gnn, double _Complex **gn0, double _Complex **g0n, double _Complex **V01, double _Complex **V10, int dim1, int dim2);
	
	void ConnectSides (double _Complex **gLL, double _Complex **gRR, double _Complex **GLR, double _Complex **GRL, double _Complex **gnn, double _Complex **gnL, double _Complex **gLn, double _Complex **GnR, double _Complex **GRn, double _Complex **VLR, double _Complex **VRL, int dim1, int dim2);

	void RibbonGF (double _Complex En, double _Complex **G, void *p);

	double getConduc (double _Complex En, void *p, double _Complex **V01, double _Complex **V10);

	void zzedgedisorder (double _Complex *DisorderArray, double potential, double concentration, int length, int size, int depth, int seed);

	void acedgedisorder (double _Complex *DisorderArray, double potential, double concentration, int length, int size, int depth, int seed);

	int **zzconnections (int size, int *numconnections);

	int **acconnections (int size, int *numconnections);
