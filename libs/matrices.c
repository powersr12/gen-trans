#include "matrices.h"

/////////////////

//MATRIX ROUTINES

/////////////////

		
	//Copy Matrix
	void MatrixCopy(double _Complex **Orig, double _Complex **Copy, int dim)
	{
		int i, j;
		
		for(i=0; i<dim;i++)
		{
			for(j=0; j<dim; j++)
			{
				Copy[i][j] = Orig[i][j];
				
			}
		}
		
	}

// //Function to copy parts of a complex matrix into parts of another complex matrix
	void MatrixCopyPart(double _Complex **Orig, double _Complex **Copy, int istartdim1, int jstartdim1,  int istartdim2, int jstartdim2, int dim1, int dim2)
	{
		int i, j;
		
		for(i=0; i<dim1;i++)
		{
			for(j=0; j<dim2; j++)
			{
				Copy[istartdim2 + i][jstartdim2 + j] = Orig[istartdim1 + i][jstartdim1 + j];
				
			}
		}
		
	}
	
	//Function to copy parts of a complex matrix into parts of another complex matrix
	void IntMatrixCopyPart(int **Orig, int **Copy, int istartdim1, int jstartdim1,  int istartdim2, int jstartdim2, int dim1, int dim2)
	{
		int i, j;
		
		for(i=0; i<dim1;i++)
		{
			for(j=0; j<dim2; j++)
			{
				Copy[istartdim2 + i][jstartdim2 + j] = Orig[istartdim1 + i][jstartdim1 + j];
				
			}
		}
		
	}

	//Routine to create a matrix of run-time specified dimension*/
	double _Complex **createSquareMatrix(int dim)	
	{
		int i,j;
		double _Complex **a;	/*pointer to a vector of vectors (ie matrix) */

		a = (double _Complex **)malloc(dim * sizeof(double _Complex *));	
					/*Asks for memory space for the array of 
					  pointers to arrays.*/

		a[0] = (double _Complex *)malloc(dim * dim * sizeof(double _Complex)); 
					/*allocates memory for whole matrix, 
					returned pointe r= pointer to 1st array*/

		for(i=1; i<dim; i++)
		{
			a[i] = a[i-1] + dim;
		}
					/*Calculates the rest of the needed pointers. 
					The pointer to the second array must be placed 
					N elements ahead of the first array etc.*/

		for(i=0; i<dim; i++)
		{
			for(j=0; j<dim; j++)
			{
				a[i][j] = 0.0 + 0.0*I;
			}
		}	

					/*Sets all matrix elements to zero*/

		return a;		/*Returns address of array of arrays*/

	
	}


//Routine to create a non-square matrix of run-time specified dimension
	double _Complex **createNonSquareMatrix(int dim, int dim2)	
	{
	int i,j;
	double _Complex **a;	/*pointer to a vector of vectors (ie matrix) */

	a = (double _Complex **)malloc(dim * sizeof(double _Complex *));	
				/*Asks for memory space for the array of 
				  pointers to arrays.*/

	a[0] = (double _Complex *)malloc(dim * dim2 * sizeof(double _Complex)); 
				/*allocates memory for whole matrix, 
				returned pointe r= pointer to 1st array*/

	for(i=1; i<dim; i++)
	{
		a[i] = a[i-1] + dim2;
	}
				/*Calculates the rest of the needed pointers. 
				The pointer to the second array must be placed 
				N elements ahead of the first array etc.*/

	for(i=0; i<dim; i++)
	{
		for(j=0; j<dim2; j++)
		{
			a[i][j] = 0.0 + 0.0*I;
		}
	}	

				/*Sets all matrix elements to zero*/

	return a;		/*Returns address of array of arrays*/

	
}

//Routine to create a non-square double matrix of run-time specified dimension
	double **createNonSquareDoubleMatrix(int dim, int dim2)	
	{
	int i,j;
	double  **a;	/*pointer to a vector of vectors (ie matrix) */

	a = (double **)malloc(dim * sizeof(double *));	
				/*Asks for memory space for the array of 
				  pointers to arrays.*/

	a[0] = (double*)malloc(dim * dim2 * sizeof(double)); 
				/*allocates memory for whole matrix, 
				returned pointe r= pointer to 1st array*/

	for(i=1; i<dim; i++)
	{
		a[i] = a[i-1] + dim2;
	}
				/*Calculates the rest of the needed pointers. 
				The pointer to the second array must be placed 
				N elements ahead of the first array etc.*/

	for(i=0; i<dim; i++)
	{
		for(j=0; j<dim2; j++)
		{
			a[i][j] = 0.0 + 0.0*I;
		}
	}	

				/*Sets all matrix elements to zero*/

	return a;		/*Returns address of array of arrays*/

	
}



	//Routine to create a matrix of run-time specified dimension*/
	//Routine to create a non-square INTEGER matrix of run-time specified dimension
	int **createNonSquareIntMatrix(int dim, int dim2)	
	{
		int i,j;
		int **a;	/*pointer to a vector of vectors (ie matrix) */
	
		a = (int **)malloc(dim * sizeof(int *));	
					/*Asks for memory space for the array of 
					pointers to arrays.*/
	
		a[0] = (int *)malloc(dim * dim2 * sizeof(int)); 
					/*allocates memory for whole matrix, 
					returned pointe r= pointer to 1st array*/
	
		for(i=1; i<dim; i++)
		{
			a[i] = a[i-1] + dim2;
		}
					/*Calculates the rest of the needed pointers. 
					The pointer to the second array must be placed 
					N elements ahead of the first array etc.*/
	
		for(i=0; i<dim; i++)
		{
			for(j=0; j<dim2; j++)
			{
				a[i][j] = 0;
			}
		}	
	
					/*Sets all matrix elements to zero*/
	
		return a;		/*Returns address of array of arrays*/
	
		
	}


	//Function to Multiply two SQUARE Matrices together
// 	void MatrixMult(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAB, int dim)
// 	{
// 		int i, j, k;
// 		EmptyMatrix(MatrixAB, dim);
// 		for(i=0; i<dim;i++)
// 		{
// 			for(j=0; j<dim; j++)
// 			{
// 				for(k=0; k<dim; k++)
// 				{
// 					MatrixAB[i][j] += (MatrixA[i][k] * MatrixB[k][j]);
// 				}
// 			}
// 		}
// 		
// 	}
	
	     
	void MatrixMult(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAB, int dim)
	{
	
	 double _Complex **C = createSquareMatrix(dim);
	 double _Complex alpha, beta;
	 alpha = 1.0;
	 beta = 0.0;
	 
	 cblas_zgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &alpha, MatrixA[0], dim, MatrixB[0], dim, &beta, C[0], dim);
	  
	  MatrixCopy(C, MatrixAB, dim);
	  FreeMatrix(C);
	}
// 	
//      int
//      main (void)
//      {
//        int lda = 3;
//      
//        float A[] = { 0.11, 0.12, 0.13,
//                      0.21, 0.22, 0.23 };
//      
//        int ldb = 2;
//        
//        float B[] = { 1011, 1012,
//                      1021, 1022,
//                      1031, 1032 };
//      
//        int ldc = 2;
//      
//        float C[] = { 0.00, 0.00,
//                      0.00, 0.00 };
//      
//        /* Compute C = A B */
//      
//        cblas_sgemm (CblasRowMajor, 
//                     CblasNoTrans, CblasNoTrans, 2, 2, 3,
//                     1.0, A, lda, B, ldb, 0.0, C, ldc);
//      
//        printf ("[ %g, %g\n", C[0], C[1]);
//        printf ("  %g, %g ]\n", C[2], C[3]);
//      
//        return 0;  
//      }

	

//Function to Multiply two NONSQUARE Matrices together
// 	void MatrixMultNS (double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAB, int r1, int c1r2, int c2)
// 	{
// 		int i, j, k;
// 		//MatrixAB = createNonSquareMatrix(r1, c2);
// 		
// 		for(i=0; i<r1;i++)
// 		{
// 			for(j=0; j<c2; j++)
// 			{
// 				MatrixAB[i][j] = 0.0;
// 				for(k=0; k<c1r2; k++)
// 				{
// 					MatrixAB[i][j] += (MatrixA[i][k] * MatrixB[k][j]);
// 				}
// 			}
// 		}
// 		
// 	}
// 	
	void MatrixMultNS(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAB, int r1, int c1r2, int c2)
	{
	
	 double _Complex **C = createNonSquareMatrix(r1, c2);
	 double _Complex alpha, beta;
	 alpha = 1.0;
	 beta = 0.0;
	 
	 cblas_zgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, r1, c2, c1r2, &alpha, MatrixA[0], c1r2, MatrixB[0], c2, &beta, C[0], c2);
	  
	  MatrixCopyPart(C, MatrixAB, 0, 0, 0, 0, r1, c2);
	  FreeMatrix(C);
	}

	//Function to Add two SQUARE Matrices together
	void MatrixAdd(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAplusB, int dim)
	{
		int i, j, k;
		EmptyMatrix(MatrixAplusB, dim);
		for(i=0; i<dim;i++)
		{
			for(j=0; j<dim; j++)
			{
				
				MatrixAplusB[i][j] = (MatrixA[i][j] + MatrixB[i][j]);
				
			}
		}

	}

		//Function to Add two NONSQUARE Matrices together
	void MatrixAddNS(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAplusB, int dim1, int dim2)
	{
		int i, j, k;
		//EmptyMatrix(MatrixAplusB, dim);
		for(i=0; i<dim1;i++)
		{
			for(j=0; j<dim2; j++)
			{
				
				MatrixAplusB[i][j] = (MatrixA[i][j] + MatrixB[i][j]);
				
			}
		}

	}
	
	//Function to Subtract one SQUARE Matrix from another
	void MatrixSubtract(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAminusB, int dim)
	{
		int i, j, k;
		EmptyMatrix(MatrixAminusB, dim);
		for(i=0; i<dim;i++)
		{
			for(j=0; j<dim; j++)
			{
				MatrixAminusB[i][j] = (MatrixA[i][j] - MatrixB[i][j]);
			}
		}

	}

	//Function to Subtract one SQUARE Matrix from another
	void MatrixSubtractNS(double _Complex **MatrixA, double _Complex **MatrixB, double _Complex **MatrixAminusB, int dim1, int dim2)
	{
		int i, j, k;
		//EmptyMatrix(MatrixAminusB, dim);
		for(i=0; i<dim1;i++)
		{
			for(j=0; j<dim2; j++)
			{
				MatrixAminusB[i][j] = (MatrixA[i][j] - MatrixB[i][j]);
			}
		}

	}
	//Function to set a matrix fully to zero
	void EmptyMatrix(double _Complex **MatrixA, int dim)
	{
		int i, j;
		for(i=0; i<dim;i++)
		{
			for(j=0; j<dim; j++)
			{
				MatrixA[i][j] = 0.0 + 0.0*I;
			}
		}
	}
	
		//Function to set a matrix fully to zero
	void EmptyDoubleMatrix(double  **MatrixA, int dim1, int dim2)
	{
		int i, j;
		for(i=0; i<dim1;i++)
		{
			for(j=0; j<dim2; j++)
			{
				MatrixA[i][j] = 0.0 ;
			}
		}
	}



	//Inverts a square complex matrix using GSL
	void InvertMatrixGSL(double _Complex **origMatrix, double _Complex **InvMatrix, int dim)
	{
		gsl_matrix_complex *orig;
		gsl_matrix_complex *inverse;
		orig = gsl_matrix_complex_calloc (dim, dim);
		inverse = gsl_matrix_complex_calloc (dim, dim);
		gsl_complex elem;
		gsl_permutation * p = gsl_permutation_calloc (dim);
		double elem1, elem2;
		double repart, impart;
		int signum;
		int i, j;

		//convert matrix to gsl format
		for(i = 0; i<dim; i++)
		{
			for(j=0; j<dim; j++)
			{	
				elem1 = creal(origMatrix[i][j]);
				elem2 = cimag(origMatrix[i][j]);
				GSL_SET_COMPLEX (&elem, elem1, elem2);
				gsl_matrix_complex_set (orig, i,  j, elem);
			}
		}

		gsl_linalg_complex_LU_decomp (orig,p, &signum);
		gsl_linalg_complex_LU_invert (orig,p, inverse);

		//convert gsl matrix to standard format
		for(i = 0; i<dim; i++)
		{
			for(j=0; j<dim; j++)
			{	
				elem = gsl_matrix_complex_get (inverse, i, j);
				repart = GSL_REAL(elem);
				impart = GSL_IMAG(elem);
				InvMatrix[i][j] = repart + impart*I;
			}
		}

		gsl_matrix_complex_free (orig);
		gsl_matrix_complex_free (inverse);
		gsl_permutation_free (p);
	}


//Inverts a square complex matrix using GSL
	double _Complex GetDeterminantGSL(double _Complex **origMatrix, int dim)
	{
		gsl_matrix_complex *orig;
		orig = gsl_matrix_complex_calloc (dim, dim);
		gsl_complex elem;
		gsl_complex determinant;
		gsl_permutation * p = gsl_permutation_calloc (dim);
		double elem1, elem2;
		double _Complex answer;
		double repart, impart;
		int signum;
		int i, j;

		//convert matrix to gsl format
		for(i = 0; i<dim; i++)
		{
			for(j=0; j<dim; j++)
			{	
				elem1 = creal(origMatrix[i][j]);
				elem2 = cimag(origMatrix[i][j]);
				GSL_SET_COMPLEX (&elem, elem1, elem2);
				gsl_matrix_complex_set (orig, i,  j, elem);
			}
		}

		gsl_linalg_complex_LU_decomp (orig,p, &signum);
		determinant = gsl_linalg_complex_LU_det (orig, signum);

		//convert gsl determinant to standard format
		repart = GSL_REAL(determinant);
		impart = GSL_IMAG(determinant);
		answer = repart + impart*I;
		
		gsl_matrix_complex_free (orig);
		gsl_permutation_free (p);
		
		return answer;
	}

	
//Returns eigenvalues of a Complex Hermitian matrix using GSL
	void EigenvaluesGSL(double _Complex **origMatrix, double _Complex *evals, int dim)
	{
		gsl_matrix_complex *orig;
		orig = gsl_matrix_complex_calloc (dim, dim);
		gsl_vector *eval = gsl_vector_alloc(dim);
		
		
		gsl_complex elem;
		
		double elem1, elem2;
		double repart, impart;
		int signum;
		int i, j;

		//convert matrix to gsl format
		for(i = 0; i<dim; i++)
		{
			for(j=0; j<dim; j++)
			{	
				elem1 = creal(origMatrix[i][j]);
				elem2 = cimag(origMatrix[i][j]);
				GSL_SET_COMPLEX (&elem, elem1, elem2);
				gsl_matrix_complex_set (orig, i,  j, elem);
			}
		}

		
		
		gsl_eigen_herm_workspace *w = gsl_eigen_herm_alloc(3*dim);
		gsl_eigen_herm(orig, eval, w);
		gsl_eigen_herm_free(w);
		gsl_sort_vector(eval);
		
		//convert eigenvalues matrix to usual format (output to complex matrix, but of course answers are real)
		for(i=0; i<dim; i++)
		{
		  elem1 = gsl_vector_get(eval, i);
		  evals[i] = elem1;
		}
	
		gsl_matrix_complex_free (orig);
		gsl_vector_free (eval);
	}
	
//Returns eigenvalues AND eigenvectors of a Complex Hermitian matrix using GSL
	void EigenvectorsGSL(double _Complex **origMatrix, double _Complex *evals, double _Complex **evecs, int dim)
	{
		gsl_matrix_complex *orig;
		orig = gsl_matrix_complex_calloc (dim, dim);
		gsl_vector *eval = gsl_vector_alloc(dim);
		gsl_matrix_complex *evec = gsl_matrix_complex_calloc (dim, dim);
		
		
		gsl_complex elem;
		
		double elem1, elem2;
		double repart, impart;
		int signum;
		int i, j;

		//convert matrix to gsl format
		for(i = 0; i<dim; i++)
		{
			for(j=0; j<dim; j++)
			{	
				elem1 = creal(origMatrix[i][j]);
				elem2 = cimag(origMatrix[i][j]);
				GSL_SET_COMPLEX (&elem, elem1, elem2);
				gsl_matrix_complex_set (orig, i,  j, elem);
			}
		}

		
		
		gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(5*dim);
		gsl_eigen_hermv(orig, eval, evec, w);
		gsl_eigen_hermv_free(w);
		gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
		
		//convert eigenvalues matrix to usual format (output to complex matrix, but of course answers are real)
		for(i=0; i<dim; i++)
		{
		  elem1 = gsl_vector_get(eval, i);
		  evals[i] = elem1;
		}
		
		//convert eigenvectors matrix to required output format
		//GSL stores in columns, we will store in rows
		for(i=0; i<dim; i++)
		{
		  for(j=0; j<dim; j++)
		  {
		    elem = gsl_matrix_complex_get (evec, j, i);
		    repart = GSL_REAL(elem);
		    impart = GSL_IMAG(elem);
		    evecs[i][j] = repart + impart*I;
		  }
		}
	
		gsl_matrix_complex_free (orig);
		gsl_matrix_complex_free (evec);
		gsl_vector_free (eval);
	}

//Function to print a Matrix (in exp form) prettily on screen*/
	void printEMatrix(double _Complex **aMatrix, int dim)

	{
		int k, l;
		printf("\n");
		for(k=0; k<dim; k++) {
			for(l=0; l<dim; l++) {
				//printf("%.1e %+.1e i  \t", creal(aMatrix[k][l]), cimag(aMatrix[k][l]));
				printf("%+.1e \t", creal(aMatrix[k][l]));
			}
			printf("\n");
		}
	}
	
	void printDMatrix(double  **aMatrix, int dim)

	{
		int k, l;
		printf("\n");
		for(k=0; k<dim; k++) {
			for(l=0; l<dim; l++) {
				//printf("%.1e %+.1e i  \t", creal(aMatrix[k][l]), cimag(aMatrix[k][l]));
				printf("%+.1e \t", aMatrix[k][l]);
			}
			printf("\n");
		}
	}
	
	//Function to print a int Matrix prettily on screen*/
	void printIntMatrix(int **aMatrix, int dim)

	{
		int k, l;
		printf("\n");
		for(k=0; k<dim; k++) {
			for(l=0; l<dim; l++) {
				//printf("% .5lf %+.5lf i  \t", creal(aMatrix[k][l]), cimag(aMatrix[k][l]));
				printf("%d\t", aMatrix[k][l]);
			}
			printf("\n");
		}
	}

//Function to take the real part of a complex matrix
	void MatrixReal(double _Complex **Orig, double _Complex **MatReal, int dim)
	{
		int i, j;
		
		for(i=0; i<dim;i++)
		{
			for(j=0; j<dim; j++)
			{
				MatReal[i][j] = creal(Orig[i][j]);
				
			}
		}
		
	}

//Function to copy a complex matrix
	double _Complex MatrixTrace(double _Complex **Matrix, int dim)
	{
		int i;
		double _Complex ans = 0.0 + 0.0*I;
		
		for(i=0; i<dim;i++)
		{
			ans += Matrix[i][i];			
			
		}

		return ans;
		
	}

int *createIntArray(int dimension)
	{
  		int i; 	/* used in a for-loop */
  		int *a; 	/* declare a pointer to a vector */

  		/* Allocate memory.*/
  		a = (int *)malloc( dimension * sizeof(int) );

  		/* Set all elements to 0.0 */
  		for ( i = 0; i < dimension; i++) 
		{
    			a[i] = 0;
  		}	
  
  		/* Return Pointer */
  		return a;

	}

	double _Complex *createCompArray(int dimension)
	{
  		int i; 	/* used in a for-loop */
  		double _Complex *a; 	/* declare a pointer to a vector */

  		/* Allocate memory.*/
  		a = (double _Complex *)malloc( dimension * sizeof(double _Complex) );

  		/* Set all elements to 0.0 */
  		for ( i = 0; i < dimension; i++) {
    			a[i] = 0.0 + 0.0*I;
  		}
  
  		/* Return Pointer */
  		return a;

	}

	double *createDoubleArray(int dimension)
	{
  		int i; 	/* used in a for-loop */
  		double  *a; 	/* declare a pointer to a vector */

  		/* Allocate memory.*/
  		a = (double *)malloc( dimension * sizeof(double) );

  		/* Set all elements to 0.0 */
  		for ( i = 0; i < dimension; i++) {
    			a[i] = 0.0 + 0.0*I;
  		}
  
  		/* Return Pointer */
  		return a;

	}
//Function to free the memory alloacted to a matrix
	void FreeMatrix(double _Complex **MatrixA)
	{
		free(MatrixA[0]);
		free(MatrixA);
	}


//returns the Green's function matrix G of a system (g) perturbed by a potential V.
	void dyson(double _Complex **g, double _Complex **V, double _Complex **G, int N)
	{
		int i, j, k;
		double _Complex **T = createSquareMatrix(N);
		double _Complex **temp1 = createSquareMatrix(N);
		double _Complex **temp2 = createSquareMatrix(N);
		double _Complex **unit = createSquareMatrix(N);

		for(i=0; i<N; i++)
			unit[i][i] = 1.0 + 0.0*I;

		MatrixMult(g, V, temp1, N);
		MatrixSubtract(unit, temp1, temp2, N);
		InvertMatrixGSL(temp2, T, N);


		MatrixMult(temp1, T, temp2, N);
		MatrixMult(temp2, g, temp1, N); 

		MatrixAdd(g, temp1, G, N);

		free(T[0]);
		free(temp1[0]);
		free(temp2[0]);
		free(unit[0]);
		free(T);
		free(temp1);
		free(temp2);
		free(unit);		
		

	}

	
	//this routine compares two complex matrices to see if they are the same
	//prints to screen the sum of the absolute values of individual matrix element differences, and the average value of this quantity
	void compareMatrices(double _Complex **MatrixA, double _Complex **MatrixB, int dim1, int dim2)
	{
	      int i, j;
	      double sum=0.0;
	      
	      for(i=0; i<dim1; i++)
	      {
		for(j=0; j<dim2; j++)
		{
		  sum += cabs(MatrixA[i][j] - MatrixB[i][j]);
		}
	      }
	      
	      printf("\n\n MATRIX DIFFERENCE: %e \t\t per element %e \n\n", sum, sum/(dim1*dim2) );
	      
	}
	
	//lists the row and column indices of nonzero elements of a complex double matrix
	//useful for checking if sparse connections matrices are correct
	void listNonZero(double _Complex **Matrix, int dim1, int dim2)
	{
	  printf("\n # NON ZERO MATRIX ELEMENTS\n");
	  int i, j;
	  for(i=0; i<dim1; i++)
	  {
	    for(j=0; j<dim2; j++)
	    {
	     // printf("%d	%d	%e\n", i, j, creal(Matrix[i][j]));
	      if(Matrix[i][j] != 0.0)
	      {
		printf("#	%d	%d	%e	%e\n", i, j, creal(Matrix[i][j]), cimag(Matrix[i][j]));
	      }
	    }
	  }
	  
	}
	

//Function to connect two disconnected systems
	void DysonConnect(double _Complex **g00, double _Complex **g11, double _Complex **V01, double _Complex **V10, double _Complex **G00, double _Complex **G11, double _Complex **G01, double _Complex **G10, int m, int edge, int n)
	{
		//g00 is the big system, g11 is a small addition, connecting to the last "edge" atoms of it
		int i;
		double _Complex **g00_small = createSquareMatrix(edge);
		double _Complex **temp1 = createNonSquareMatrix(n, edge);
		double _Complex **temp2 = createNonSquareMatrix(n, edge);
		double _Complex **temp3 = createSquareMatrix(n);
		double _Complex **temp4 = createSquareMatrix(n);
		double _Complex **temp5 = createNonSquareMatrix(n, m);
		double _Complex **temp6 = createNonSquareMatrix(m, n);
		double _Complex **temp7 = createSquareMatrix(m);

		

		double _Complex **unit = createSquareMatrix(n);
			for(i=0; i<n; i++)
			{
				unit[i][i] = 1.0 + 0.0*I;
			}

		//printf("connectok1\n");
		//calc G11 - connected new edge
		MatrixCopyPart(g00, g00_small, m - edge, m - edge, 0, 0, edge, edge);
		MatrixMultNS(g11, V10, temp1, n, n, edge);
		MatrixMultNS(temp1, g00_small, temp2, n, edge, edge);
		MatrixMultNS(temp2, V01, temp3, n, edge, n);
		MatrixSubtract(unit, temp3, temp4, n);
		
		//listNonZero(temp4, n, n);
		InvertMatrixGSL(temp4, temp3, n);
		MatrixMultNS(temp3, g11, G11, n, n, n);
		//printf("connectok2\n");


		double _Complex **V10_big = createNonSquareMatrix(n, m);
		double _Complex **V01_big = createNonSquareMatrix(m, n);

		MatrixCopyPart(V10, V10_big, 0, 0, 0, m-edge, n, edge);
		MatrixCopyPart(V01, V01_big, 0, 0, m-edge, 0, edge, n);

		int  j;
		
		
		
		//calc G10
		MatrixMultNS(G11, V10_big, temp5, n, n, m);
		MatrixMultNS(temp5, g00, G10, n, m , m);
		
	
		
		
		//calc G01
		MatrixMultNS(g00, V01_big, temp6, m, m, n);
		MatrixMultNS(temp6, G11, G01, m, n, n);

		//calc G00
		MatrixMultNS(temp6, G10, temp7, m, n, m);
		MatrixAdd(g00, temp7, G00, m);

		FreeMatrix(g00_small);
		FreeMatrix(temp1);
		FreeMatrix(temp2);
		FreeMatrix(temp3);
		FreeMatrix(temp4);
		FreeMatrix(temp5);
		FreeMatrix(temp6);
		FreeMatrix(temp7);
		FreeMatrix(unit);
		FreeMatrix(V10_big);
		FreeMatrix(V01_big);


	}

	void DysonConnect2(double _Complex **g00, double _Complex **g11, double _Complex **V01, double _Complex **V10, double _Complex **G00, double _Complex **G11, double _Complex **G01, double _Complex **G10, int m, int edge, int n)
	{
		//g00 is the big system, g11 is a small addition, connecting to the last "edge" atoms of it
		
		double _Complex **gbig = createSquareMatrix(m+n);
		double _Complex **Gbig = createSquareMatrix(m+n);
		double _Complex **Vbig = createSquareMatrix(m+n);

		MatrixCopyPart(g00, gbig, 0, 0, 0, 0, m, m);
		MatrixCopyPart(g11, gbig, 0, 0, m, m, n, n);
		MatrixCopyPart(V01, Vbig, 0, 0, m-edge, m, edge, n);
		MatrixCopyPart(V10, Vbig, 0, 0, m, m-edge, n, edge);
		
		dyson(gbig, Vbig, Gbig, m+n);
		
		MatrixCopyPart(Gbig, G00, 0, 0, 0, 0, m, m);
		MatrixCopyPart(Gbig, G11, m, m, 0, 0, n, n);
		MatrixCopyPart(Gbig, G01, 0, m, 0, 0, m, n);
		MatrixCopyPart(Gbig, G10, m, 0, 0, 0, n, m);
		


		FreeMatrix(gbig);
		FreeMatrix(Gbig);
		FreeMatrix(Vbig);
		

	}

	
	//different number of 'edge' sites in both structures
	void DysonConnect3(double _Complex **g00, double _Complex **g11, double _Complex **V01, double _Complex **V10, double _Complex **G00, double _Complex **G11, double _Complex **G01, double _Complex **G10, int m, int edge, int edge2, int n)
	{
		//g00 is the big system, g11 is a small addition, connecting to the last "edge" atoms of it
		
		double _Complex **gbig = createSquareMatrix(m+n);
		double _Complex **Gbig = createSquareMatrix(m+n);
		double _Complex **Vbig = createSquareMatrix(m+n);

		MatrixCopyPart(g00, gbig, 0, 0, 0, 0, m, m);
		MatrixCopyPart(g11, gbig, 0, 0, m, m, n, n);
		MatrixCopyPart(V01, Vbig, 0, 0, m-edge, m, edge, edge2);
		MatrixCopyPart(V10, Vbig, 0, 0, m, m-edge, edge2, edge);
				//printEMatrix(gbig, m+n);

		dyson(gbig, Vbig, Gbig, m+n);
		
		MatrixCopyPart(Gbig, G00, 0, 0, 0, 0, m, m);
		MatrixCopyPart(Gbig, G11, m, m, 0, 0, n, n);
		MatrixCopyPart(Gbig, G01, 0, m, 0, 0, m, n);
		MatrixCopyPart(Gbig, G10, m, 0, 0, 0, n, m);
		


		FreeMatrix(gbig);
		FreeMatrix(Gbig);
		FreeMatrix(Vbig);
		

	}
	
	void LinEqnDouble (double **MatrixA, double *VectorB, double *VectorX, int dim)
	{
	  int i, j, k;
	  gsl_matrix *orig;
	  orig = gsl_matrix_calloc(dim, dim);
	  gsl_vector *b, *x;
	  b = gsl_vector_calloc(dim);
	  x = gsl_vector_calloc(dim);
	  
	  for(i=0; i<dim; i++)
	  {
	    for(j=0; j<dim; j++)
	    {
	      gsl_matrix_set(orig, i, j, MatrixA[i][j]);
	    }
	    gsl_vector_set(b, i, VectorB[i]);
	  }
	  
	  int s;
	  gsl_permutation * p = gsl_permutation_calloc (dim);
	  
	  gsl_linalg_LU_decomp(orig, p, &s);
	  gsl_linalg_LU_solve(orig, p, b, x);
	  
	  
	  for(i=0; i<dim; i++)
	  {
	    VectorX[i] = gsl_vector_get(x, i);
	  }
	  
	  gsl_matrix_free(orig);
	  gsl_vector_free(b);
	  gsl_vector_free(x);
	  gsl_permutation_free(p);
	  
	}
	  
	  
	    
	  
	  
	  
	  
	  
	  
// 	  	gsl_matrix_complex *orig;
// 		gsl_matrix_complex *inverse;
// 		orig = gsl_matrix_complex_calloc (dim, dim);
// 		inverse = gsl_matrix_complex_calloc (dim, dim);
// 		gsl_complex elem;
// 		gsl_permutation * p = gsl_permutation_calloc (dim);
// 		double elem1, elem2;
// 		double repart, impart;
// 		int signum;
// 		int i, j;
// 
// 		//convert matrix to gsl format
// 		for(i = 0; i<dim; i++)
// 		{
// 			for(j=0; j<dim; j++)
// 			{	
// 				elem1 = creal(origMatrix[i][j]);
// 				elem2 = cimag(origMatrix[i][j]);
// 				GSL_SET_COMPLEX (&elem, elem1, elem2);
// 				gsl_matrix_complex_set (orig, i,  j, elem);
// 			}
// 		}
// 
// 		gsl_linalg_complex_LU_decomp (orig,p, &signum);
// 		gsl_linalg_complex_LU_invert (orig,p, inverse);
// 
// 		//convert gsl matrix to standard format
// 		for(i = 0; i<dim; i++)
// 		{
// 			for(j=0; j<dim; j++)
// 			{	
// 				elem = gsl_matrix_complex_get (inverse, i, j);
// 				repart = GSL_REAL(elem);
// 				impart = GSL_IMAG(elem);
// 				InvMatrix[i][j] = repart + impart*I;
// 			}
// 		}

	  
	  
	