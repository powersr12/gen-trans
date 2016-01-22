#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_cblas.h>




void MatrixCopy(double _Complex **Orig, double _Complex **Copy, int dim);
void MatrixCopyPart(double _Complex **Orig, double _Complex **Copy, int istartdim1, int jstartdim1,  int istartdim2, int jstartdim2, int dim1, int dim2);
void IntMatrixCopyPart(int **Orig, int **Copy, int istartdim1, int jstartdim1,  int istartdim2, int jstartdim2, int dim1, int dim2);

double _Complex **createSquareMatrix(int dim)	;
int **createNonSquareIntMatrix(int dim, int dim2);	
double _Complex **createNonSquareMatrix(int dim, int dim2);	
double **createNonSquareDoubleMatrix(int dim, int dim2)	;

void MatrixMult(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAB, int dim);
void MatrixMultNS (double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAB, int r1, int c1r2, int c2);
void MatrixAdd(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAplusB, int dim);
void MatrixAddNS(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAplusB, int dim1, int dim2);
void MatrixSubtract(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAminusB, int dim);
void MatrixSubtractNS(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAminusB, int dim1, int dim2);

void EmptyMatrix(double _Complex **MatrixA, int dim);
	void EmptyDoubleMatrix(double  **MatrixA, int dim1, int dim2);

void InvertMatrixGSL(double _Complex **origMatrix, double _Complex **InvMatrix, int dim);
double _Complex GetDeterminantGSL(double _Complex **origMatrix, int dim);
void EigenvaluesGSL(double _Complex **origMatrix, double _Complex *evals, int dim);
void EigenvectorsGSL(double _Complex **origMatrix, double _Complex *evals, double _Complex **evecs, int dim);


void printEMatrix(double _Complex **aMatrix, int dim);
void printIntMatrix(int **aMatrix, int dim);

double _Complex MatrixTrace(double _Complex **Matrix, int dim);
void MatrixReal(double _Complex **Orig, double _Complex **MatReal, int dim);
int *createIntArray(int dimension);
double _Complex *createCompArray(int dimension);
double *createDoubleArray(int dimension);
void FreeMatrix(double _Complex **MatrixA);
void dyson(double _Complex **g, double _Complex **V, double _Complex **G, int N);
	void DysonConnect(double _Complex **g00, double _Complex **g11, double _Complex **V01, double _Complex **V10, double _Complex **G00, double _Complex **G11, double _Complex **G01, double _Complex **G10, int m, int edge, int n);
		void DysonConnect2(double _Complex **g00, double _Complex **g11, double _Complex **V01, double _Complex **V10, double _Complex **G00, double _Complex **G11, double _Complex **G01, double _Complex **G10, int m, int edge, int n);
	void DysonConnect3(double _Complex **g00, double _Complex **g11, double _Complex **V01, double _Complex **V10, double _Complex **G00, double _Complex **G11, double _Complex **G01, double _Complex **G10, int m, int edge, int edge2, int n);

	
	
	void compareMatrices(double _Complex **MatrixA, double _Complex **MatrixB, int dim1, int dim2);
	void listNonZero(double _Complex **Matrix, int dim1, int dim2);

	void LinEqnDouble (double **MatrixA, double *VectorB, double *VectorX, int dim);
