
#include "greenfns.h"





/*returns the (up AND down spin) Green Fn matrices required for the integrands of I1, I2 and I3.
		-gfroutine is a function of the form gfroutine (En, **G, *gfparams) 
		-gfparams is a parameter structure containing anything else needed to calculate the GFs (positions, size, disorder locations etc)
		-N is the dimension of the Green function matrix.*/

// 	void MagGreen(double _Complex upenergy, double _Complex downenergy, double _Complex **Gup, double _Complex **Gdown, double *m, double *shift, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, double hubU, int N, double *zeemanterm)
// 	{
// 			//Memory declarations
// 			double _Complex **Gun= createSquareMatrix(N);
// 			double _Complex **Vup= createSquareMatrix(N);
// 			double _Complex **Vdown= createSquareMatrix(N);
// 			int count, i, j;
// 
// 		//define Vup and Vdown
// 			for(i=0; i<N; i++)
// 			{
// 				Vup[i][i] = -(hubU/2)*(m[i]) + shift[i] -zeemanterm[i]/2; 
// 				Vdown[i][i] = (hubU/2)*(m[i]) + shift[i] +zeemanterm[i]/2;
// 			}
// 		//printf("ok\n");
// 
// 		//Gup
// 			(*gfroutine)(upenergy, Gun, gfparams);			
// 			dyson(Gun, Vup, Gup, N);
// 
// 		//Gdown
// 			(*gfroutine)(downenergy, Gun, gfparams);			
// 			dyson(Gun, Vdown, Gdown, N);
// 
// 		//Free up memory
// 			FreeMatrix(Gun);
// 			FreeMatrix(Vup);
// 			FreeMatrix(Vdown);
// 
// 
// 	}


/*Used to calculate the induced magnetisation at a site when a spin-dependent potential is applied (usually elsewhere)
		-integrates using an integral calculated using intdos_other.
		-up and down spin potential matrices are calculated from the arrays *m, *shift and *zeemanterm
		-gfroutine is a function of the form gfroutine (En, **G, *gfparams) 
		-gfparams is a parameter structure containing anything else needed to calculate the GFs (positions, size, disorder locations etc)
		-uses the parameter structure "intdosother_params", defined in greenfns.h and containing the elements
			int i; int N; double fermi; double imagEn; double _Complex **V; void (*gfroutine2)(double _Complex , double _Complex **, void * ); void *gfparams;
		-i is the site index where the local moment is being calculated
		-N is the dimension of the Green function matrix.	*/

// 	double getlocalm (double *m, double *shift, double *zeemanterm, double hubU, void (*gfroutine2) (double _Complex, double _Complex **, void *), void *gfparams, int i, int numsites, double fermi, double imagEn)
// 	{
// 		intdosother_params para = {};
// 	
// 		double nup, ndown, error1, error2;
// 		int j;
// 	
// 		//DOS integrator - initialisation
// 		gsl_function integ;
// 		integ.function = &intdos_other;
// 		integ.params = &para;
// 	
// 		//perturbation matrices for up and down spin electrons
// 		double _Complex **Vup = createSquareMatrix(numsites);
// 		double _Complex **Vdown = createSquareMatrix(numsites);
// 	
// 		//parameters passed to DOS integrator
// 		para.i = i;
// 		para.N = numsites;
// 		para.fermi = fermi;
// 		para.imagEn = imagEn;
// 		para.gfroutine2 = gfroutine2;
// 		para.gfparams = gfparams;
// 		
// 	
// 		//set-up perturbation potential
// 		for(j=0; j<numsites; j++)
// 		{
// 			Vup[j][j] = 	-(hubU/2)*(m[j]) + shift[j] -zeemanterm[j]/2;
// 			Vdown[j][j] = (hubU/2)*(m[j]) + shift[j] -zeemanterm[j]/2;
// 		}
// 	
// 	
// 		gsl_integration_workspace *w;
// 		w = gsl_integration_workspace_alloc (100000);
// 		para.V = Vup;
// 		gsl_integration_qags (&integ, 0.0, 1.0, 0.00000, 0.0001, 100000,  w, &nup, &error1);
// 		nup += 0.5;
// 	
// 		para.V = Vdown;
// 		gsl_integration_qags (&integ, 0.0, 1.0, 0.00000, 0.0001, 100000,  w, &ndown, &error2);
// 		ndown += 0.5;
// 	
// 		gsl_integration_workspace_free (w);
// 	
// 	
// 	
// 		return nup - ndown;
// 	}


/*calculates the integrand used in getlocalm
		-uses the parameter structure "intdosother_params", defined in greenfns.h and containing the elements
			int i; int N; double fermi; double imagEn; double _Complex **V; void (*gfroutine2)(double _Complex , double _Complex **, void * ); void *gfparams;
		-i is the site index where the local moment is being calculated */

// 	double intdos_other(double impart, void *p)
// 	{
// 		intdosother_params *params = (intdosother_params *)p;
// 
// 		double ans;
// 
// 		int i = (params->i);	//which element to return
// 		int N = (params -> N);	//matrix size
// 		double imagEn = (params->imagEn);
// 		double fermi = (params->fermi);
// 		double _Complex **V = (params->V);
// 		void (*gfroutine2)(double _Complex, double _Complex **, void *) = (params->gfroutine2);
// 		void *gfparams = (params->gfparams);
// 
// 
// 		double _Complex **G = createSquareMatrix(N);
// 		double _Complex **Gpert = createSquareMatrix(N);
// 
// 		double _Complex En = fermi + ((1 + imagEn - impart)/impart)*I;
// 
// 		(*gfroutine2) (En, G, gfparams);
// 		
// 		dyson(G, V, Gpert, N);
// 		ans = ((1 + imagEn)/M_PI) * (creal(Gpert[i][i])) / (impart*impart) ;
// 
// 
// 		FreeMatrix(Gpert);
// 		FreeMatrix(G);
// 		return ans;
// 	}


//creates RubioSGF inputs for a semi-infinite graphene sheet with armchair edge 
// indices 0 and 1 are on the very edge, 2 and 3 next to edge
//al matrices are 4x4 and should be declared externally beforehand
//kpar is the parallel k-vector that should be integrated over eventually to get the full GF

void ACinf (double _Complex **g11inv, double _Complex **V01, double _Complex **V10, double kpar, double hopping, double _Complex Eng)
	{
		int i, j;
		EmptyMatrix(g11inv, 4);
		EmptyMatrix(V01,4);
		EmptyMatrix(V10, 4);

		g11inv[0][1] = -hopping;
		g11inv[0][3] = -hopping * cexp(I*kpar);
		g11inv[1][0] = -hopping;
		g11inv[1][2] = -hopping;
		g11inv[2][1] = -hopping;
		g11inv[2][3] = -hopping;
		g11inv[3][0] = -hopping * cexp(-I*kpar);
		g11inv[3][2] = -hopping;

		g11inv[0][0] = Eng;
		g11inv[1][1] = Eng;
		g11inv[2][2] = Eng;
		g11inv[3][3] = Eng;

		V01[2][1] = hopping;
		V01[3][0] = hopping * cexp(-I*kpar);
		V10[0][3] = hopping * cexp(I*kpar);
		V10[1][2] = hopping;
		

	}





/*Uses the Rubio-Sancho (transfer-matrix) approach to calculate the surface Green function
		-inputs required describe the unit cell and are calculated from another routine (AGNR, ZGNR etc)
		-g00 - green function of unit cell
		-V01 and V10 - connection matrices between neighbouring cells
		-note that here "00" denotes the surface, so V01 and V10 have opposite syntax than usual*/
	
	void RubioSGF (double _Complex **S, double _Complex **g00, double _Complex **V01, double _Complex **V10, int N, int *count, double error)
	{

		int i, j;
		double calc_error = error + 1.0;

		double _Complex **unit = createSquareMatrix(N);

		for(i=0; i<N; i++)
		{
			unit[i][i] = 1.0;
		}

		double _Complex **temp1 = createSquareMatrix(N);
		double _Complex **temp2 = createSquareMatrix(N);
		double _Complex **temp3 = createSquareMatrix(N);

		double _Complex **t_old = createSquareMatrix(N);
		double _Complex **t_new = createSquareMatrix(N);
		double _Complex **td_old = createSquareMatrix(N);
		double _Complex **td_new = createSquareMatrix(N);

		double _Complex **T_old = createSquareMatrix(N);
		double _Complex **T_new = createSquareMatrix(N);
		double _Complex **B_old = createSquareMatrix(N);
		double _Complex **B_new = createSquareMatrix(N);


		MatrixMult(g00, V10, t_old, N);
		MatrixMult(g00, V01, td_old, N);
		MatrixCopy(t_old, T_old, N);
		MatrixCopy(td_old, B_old, N);

		*count = 0;
		while(calc_error > error)
		{
			*count = *count + 1;
			calc_error = 0.0;
			//t_new
			MatrixMult(t_old, td_old, temp1, N);
			MatrixSubtract(unit, temp1, temp2, N);
			MatrixMult(td_old, t_old, temp1, N);
			MatrixSubtract(temp2, temp1, temp3, N);
			InvertMatrixGSL(temp3, temp1, N);

			MatrixMult(temp1, t_old, temp2, N);
			MatrixMult(temp2, t_old, t_new, N);


			//td_new
			MatrixMult(temp1, td_old, temp2, N);
			MatrixMult(temp2, td_old, td_new, N);			
			//printEMatrix(td_new, 4);
			//printEMatrix(td_new, N);

			
			//T_new
			MatrixMult(B_old, t_new, temp1, N);
			MatrixAdd(T_old, temp1, T_new, N);

			//B_new
			MatrixMult(B_old, td_new, B_new, N);

			for(i=0; i < N; i++)	
			{
				for(j=0; j<N; j++)
				{
					calc_error += cabs(t_new[i][j] - t_old[i][j]);
					calc_error += cabs(td_new[i][j] - td_old[i][j]);
				}
			}

			//printf("# calc_error: %e	allowed error: %e \n", calc_error, error);
			MatrixCopy(t_new, t_old, N);
			MatrixCopy(td_new, td_old, N);
			MatrixCopy(T_new, T_old, N);
			MatrixCopy(B_new, B_old, N);

		}
		
		MatrixMult(g00, V01, temp1, N);
		MatrixMult(temp1, T_new, temp2, N);
		MatrixSubtract(unit, temp2, temp1, N);
		InvertMatrixGSL(temp1, temp2, N);
		MatrixMult(temp2, g00, S, N);

		FreeMatrix(unit);
		FreeMatrix(temp1);
		FreeMatrix(temp2);
		FreeMatrix(temp3);
		FreeMatrix(t_old);
		FreeMatrix(t_new);
		FreeMatrix(td_old);
		FreeMatrix(td_new);
		FreeMatrix(T_old);
		FreeMatrix(T_new);
		FreeMatrix(B_old);
		FreeMatrix(B_new);
	}	

/*calculates unit cell INVERSE gf, and connection matrices for a zigzag GNR, for use in Rubio SGF routine and building disordered regions
		-g11inv must be inverted before being sent to RubioSGF
		-this form allows easier inclusion of disorder when building disordered regions*/

	void ZGNR (double _Complex **g11inv, double _Complex **V01, double _Complex **V10, int N, double hopping, double _Complex Eng)
	{
		int i, j;
		EmptyMatrix(g11inv, 2*N);
		EmptyMatrix(V01, 2*N);
		EmptyMatrix(V10, 2*N);
		g11inv[0][1] = -hopping;
		g11inv[2*N -1][2*N-2] = -hopping;				
		g11inv[0][0] = Eng;
		g11inv[2*N -1][2*N -1] = Eng;
		V01[0][1] = hopping;

		for(i=1; i<2*N-1; i++)
		{
			g11inv[i][i] = Eng;
			g11inv[i][i+1] = -hopping;
			g11inv[i][i-1] = -hopping;
		}

		for(i=1; i<2*N-3; i+=4)
		{
			V10[i][i-1] = hopping;
			V10[i+1][i+2] = hopping;
			V01[i+2][i+1] = hopping;
			V01[i+3][i+4] = hopping;
		}

		if(N % 2 == 0)
		{
			V10[2*N-3][2*N-4] = hopping;
			V10[2*N-2][2*N-1] = hopping;
			V01[2*N-1][2*N-2] = hopping;
		}

		if(N % 2 == 1)
		{
			V10[2*N-1][2*N-2] = hopping;
		}	

	}

/*calculates unit cell INVERSE gf, and connection matrices for a armchair GNR, for use in Rubio SGF routine and building disordered regions
		-g11inv must be inverted before being sent to RubioSGF
		-this form allows easier inclusion of disorder when building disordered regions
		-unit cell constructed so that 	sites 0 and N are not connected - i.e. constitute the open part of an edge
						sites N and 2N ARE connected - i.e. the closed part of an edge 


					0		N
				_______			  _________
					\		/
					  \	      /	
					    \_______/

					     1	    N+1
	
*/

	void AGNR (double _Complex **g11inv, double _Complex **V01, double _Complex **V10, int N, double hopping, double _Complex Eng)
	{
		int i, j;
		EmptyMatrix(g11inv, 2*N);
		EmptyMatrix(V01, 2*N);
		EmptyMatrix(V10, 2*N);

		g11inv[0][1] = -hopping;
		g11inv[N][N+1] = -hopping;
		g11inv[N -1][N-2] = -hopping;
		g11inv[2*N -1][2*N-2] = -hopping;

		g11inv[0][0] = Eng;
		g11inv[N-1][N-1] = Eng;

		g11inv[N][N] = Eng;
		g11inv[2*N -1][2*N -1] = Eng;

		V10[0][N] = hopping;
		V01[N][0] = hopping;

		for(i=1; i<N-1; i++)
		{
			g11inv[i][i] = Eng;
			g11inv[N+i][N+i] = Eng;

			g11inv[i][i+1] = -hopping;
			g11inv[N+i][N+i+1] = -hopping;

			g11inv[i][i-1] = -hopping;
			g11inv[N+i][N+i-1] = -hopping;

		}
		
		for(i=1; i<N-1; i+=2)
		{
			
			g11inv[i][N+i] = -hopping;
			g11inv[N+i][i] = -hopping;

			V10[i+1][N+i+1] = hopping;
			V01[N+i+1][i+1] = hopping;

		}

		if(N % 2 == 0)
		{
			g11inv[N-1][2*N-1] = -hopping;
			g11inv[2*N-1][N-1] = -hopping;
		}
		

	}

/*calculates unit cell INVERSE gf, and connection matrices for a armchair GNR, for use in Rubio SGF routine and building disordered regions
		-g11inv must be inverted before being sent to RubioSGF
		-unit cell construction as per AGNR, EXCEPT:
			12% increase in t parameter for appropriate edges, corresponding to edge reconfiguration
			should lead to AGNRs always being semiconducting */

	void AGNR_E (double _Complex **g11inv, double _Complex **V01, double _Complex **V10, int N, double hopping, double _Complex Eng)
	{
		int i, j;
		EmptyMatrix(g11inv, 2*N);
		EmptyMatrix(V01, 2*N);
		EmptyMatrix(V10, 2*N);

		g11inv[0][1] = -hopping;
		g11inv[N][N+1] = -hopping;
		g11inv[N -1][N-2] = -hopping;
		g11inv[2*N -1][2*N-2] = -hopping;


				
		g11inv[0][0] = Eng;
		g11inv[N-1][N-1] = Eng;

		g11inv[N][N] = Eng;
		g11inv[2*N -1][2*N -1] = Eng;

		V10[0][N] = 1.12*hopping;
		V01[N][0] = 1.12*hopping;

		for(i=1; i<N-1; i++)
		{
			g11inv[i][i] = Eng;
			g11inv[N+i][N+i] = Eng;

			g11inv[i][i+1] = -hopping;
			g11inv[N+i][N+i+1] = -hopping;

			g11inv[i][i-1] = -hopping;
			g11inv[N+i][N+i-1] = -hopping;

		}

		for(i=1; i<N-2; i+=2)
		{
			
			g11inv[i][N+i] = -hopping;
			g11inv[N+i][i] = -hopping;

			V10[i+1][N+i+1] = hopping;
			V01[N+i+1][i+1] = hopping;

		}

		if(N % 2 == 1)
		{
			V10[N-1][2*N-1] = 1.12*hopping;
			V01[2*N-1][N-1] = 1.12*hopping;
		}

		if(N % 2 == 0)
		{
			g11inv[N-1][2*N-1] = -1.12*hopping;
			g11inv[2*N-1][N-1] = -1.12*hopping;
		}
		

	}


//returns GF for selected sites in a GNR or CNT unit cell
	//acts as a wrapper for Rubio or other GF routine
	void GFsites (double _Complex En, double _Complex **G, void *p)
	{
		GFsites_params *params = (GFsites_params *)p;
		
		int N = (params->size);
		int numsites = (params->numsites);
		int *pos = (params->pos);
		double error = (params->rubio_error);
		double hopping = (params->hopping);
		int count;

		int i, j;
		double _Complex **g11i = createSquareMatrix(2*N);
		double _Complex **g11 = createSquareMatrix(2*N);
		double _Complex **SL = createSquareMatrix(2*N);
		double _Complex **SR = createSquareMatrix(2*N);
		double _Complex **V01 = createSquareMatrix(2*N);
		double _Complex **V10 = createSquareMatrix(2*N);
		double _Complex **G00 = createSquareMatrix(2*N);
		double _Complex **temp1 = createSquareMatrix(2*N);
		double _Complex **temp2 = createSquareMatrix(2*N);
		double _Complex **unit = createSquareMatrix(2*N);
			for(j=0; j<2*N; j++)
			{
				unit[j][j] = 1.0;
			}

			(params->gfroutine) (g11i, V01, V10, N, hopping, En);
			InvertMatrixGSL(g11i, g11, 2*N);
			RubioSGF(SL, g11, V10, V01, 2*N, &count, error);
			RubioSGF(SR, g11, V01, V10, 2*N, &count, error);

			MatrixMult(SL, V01, temp1, 2*N);
			MatrixMult(temp1, SR, temp2, 2*N);
			MatrixMult(temp2, V10, temp1, 2*N);
			MatrixSubtract(unit, temp1, temp2, 2*N);
			InvertMatrixGSL(temp2, temp1, 2*N);
			MatrixMult(temp1, SL, G00, 2*N);


		for(i=0; i< numsites; i++)
		{
			for(j=0; j< numsites; j++)
			{
				G[i][j] = G00[pos[i]][pos[j]] ;
			}
		}


		FreeMatrix(g11i);
		FreeMatrix(g11);
		FreeMatrix(SL);
		FreeMatrix(SR);
		FreeMatrix(V01);
		FreeMatrix(V10);
		FreeMatrix(G00);
		FreeMatrix(temp1);
		FreeMatrix(temp2);
		FreeMatrix(unit);


	}
	
	//returns GF for selected sites in a semi-inf graphene sheet
	//acts as a wrapper for Rubio or other GF routine
	void GFsites2 (double _Complex En, double _Complex **G, void *p)
	{
		GFsites2_params *params = (GFsites2_params *)p;
		
		int N = (params->size);
		int numsites = (params->numsites);
		int *pos = (params->pos);
		double error = (params->rubio_error);
		double hopping = (params->hopping);
		double kpar = (params->kpar);
		int count;

		int i, j;
		double _Complex **g11i = createSquareMatrix(2*N);
		double _Complex **g11 = createSquareMatrix(2*N);
		double _Complex **SL = createSquareMatrix(2*N);
		double _Complex **SR = createSquareMatrix(2*N);
		double _Complex **V01 = createSquareMatrix(2*N);
		double _Complex **V10 = createSquareMatrix(2*N);
		double _Complex **G00 = createSquareMatrix(2*N);
		double _Complex **temp1 = createSquareMatrix(2*N);
		double _Complex **temp2 = createSquareMatrix(2*N);
		double _Complex **unit = createSquareMatrix(2*N);
			for(j=0; j<2*N; j++)
			{
				unit[j][j] = 1.0;
			}

			(params->gfroutine) (g11i, V01, V10, kpar, hopping, En);
			//printEMatrix(g11i, 4);
			//printEMatrix(V01, 4);
			//printEMatrix(V10, 4);


			InvertMatrixGSL(g11i, g11, 2*N);

			RubioSGF(SL, g11, V10, V01, 2*N, &count, error);
			//comment out next couple of lines for a semi infinite system
			RubioSGF(SR, g11, V01, V10, 2*N, &count, error);


			MatrixMult(SL, V01, temp1, 2*N);
			MatrixMult(temp1, SR, temp2, 2*N);
			MatrixMult(temp2, V10, temp1, 2*N);
			MatrixSubtract(unit, temp1, temp2, 2*N);
			InvertMatrixGSL(temp2, temp1, 2*N);
			MatrixMult(temp1, SL, G00, 2*N);


		for(i=0; i< numsites; i++)
		{
			for(j=0; j< numsites; j++)
			{
				G[i][j] = G00[pos[i]][pos[j]] ;
			}
		}


		FreeMatrix(g11i);
		FreeMatrix(g11);
		FreeMatrix(SL);
		FreeMatrix(SR);
		FreeMatrix(V01);
		FreeMatrix(V10);
		FreeMatrix(G00);
		FreeMatrix(temp1);
		FreeMatrix(temp2);
		FreeMatrix(unit);


	}
	
	//returns kpar-dependent integrand for selected sites (on the edge of...) a semi-inf graphene sheet
	//acts as a wrapper for Rubio or other GF routine
	double GFsites2int (double kpar, void *p)
	{
		GFsites2int_params *params = (GFsites2int_params *)p;
		
		int N = (params->size);
		int numsites = (params->numsites);
		int *pos = (params->pos);
		double error = (params->rubio_error);
		double hopping = (params->hopping);
		double _Complex En = (params->En);
		int cell1_edge_dist = (params->cell1_edge_dist);	
		int cell2_edge_dist = (params->cell2_edge_dist);
		int cell_para_dist = (params->cell_para_dist);
		int cell1site = (params->cell1site);
		int cell2site = (params->cell2site);
		int reim = (params->reim);

		int count, swap, numsep;
		
		//make cell 2 the furthest from the edge
		if(cell2_edge_dist < cell1_edge_dist)
		{
		  swap =cell2_edge_dist;
		  cell2_edge_dist = cell1_edge_dist;
		  cell1_edge_dist = swap;
		  swap = cell2site;
		  cell2site = cell1site;
		  cell1site = cell2site;
		}

		int i, j, size_of_GNN;
		double _Complex **G = createSquareMatrix(2*N);
		double _Complex **GLR = createSquareMatrix(2*N);
		double _Complex **GRL = createSquareMatrix(2*N);
		double _Complex **G00 = createSquareMatrix(2*N);
		double _Complex **G11 = createSquareMatrix(2*N);
		double _Complex **Gnn = createSquareMatrix(2*N);
		double _Complex **Gnn_old = createSquareMatrix(2*N);
		double _Complex **GnL_old = createSquareMatrix(2*N);
		double _Complex **GLn_old = createSquareMatrix(2*N);
		double _Complex **GnR = createSquareMatrix(2*N);
		double _Complex **GRn = createSquareMatrix(2*N);
		
		double _Complex **bigGnn_old = createSquareMatrix(4*N);
		double _Complex **bigGnL_old = createNonSquareMatrix(4*N, 2*N);
		double _Complex **bigGLn_old = createNonSquareMatrix(2*N, 4*N);
		double _Complex **bigGnR = createNonSquareMatrix(4*N,2*N);
		double _Complex **bigGRn = createNonSquareMatrix(2*N, 4*N);
		
		
		double _Complex **g11i = createSquareMatrix(2*N);
		double _Complex **g11 = createSquareMatrix(2*N);
		double _Complex **SL = createSquareMatrix(2*N);
		double _Complex **SR = createSquareMatrix(2*N);
		double _Complex **V01 = createSquareMatrix(2*N);
		double _Complex **V10 = createSquareMatrix(2*N);
		double _Complex **temp1 = createSquareMatrix(2*N);
		double _Complex **temp2 = createSquareMatrix(2*N);
		double _Complex **unit = createSquareMatrix(2*N);
		
		
			for(j=0; j<2*N; j++)
			{
				unit[j][j] = 1.0;
			}

			(params->gfroutine) (g11i, V01, V10, kpar, hopping, En);
			//printEMatrix(g11i, 4);
			//printEMatrix(V01, 4);
			//printEMatrix(V10, 4);


			InvertMatrixGSL(g11i, g11, 2*N);

			
			//SGF with outer cell becoming cell 2 (furthest from edge in final calculation)
			RubioSGF(SL, g11, V10, V01, 2*N, &count, error);
			
			MatrixCopy(SL, G00, 2*N);
			MatrixCopy(g11, G11, 2*N);

			numsep = cell2_edge_dist - cell1_edge_dist;
			
			//add blocks between the two cells of interest, and then the second cell  
			if(numsep > 0)
			{
			  //connect first cell
			  ConnectSides(G00, G11, GLR, GRL, Gnn_old, GnL_old, GLn_old, GnR, GRn, V01, V10, 2*N, 0);

			  MatrixCopy(G00, Gnn_old, 2*N);
			  MatrixCopy(GLR, GnL_old, 2*N);
			  MatrixCopy(GRL, GLn_old, 2*N);
			  MatrixCopy(G11, G00, 2*N);
			  MatrixCopy(g11, G11, 2*N);
			  
			  //connect remaining cells
			  for(i=1; i<numsep; i++)
			  {
				ConnectSides(G00, G11, GLR, GRL, Gnn_old, GnL_old, GLn_old, GnR, GRn, V01, V10, 2*N, 2*N);
				MatrixCopy(GnR, GnL_old, 2*N);
				MatrixCopy(GRn, GLn_old, 2*N);
				MatrixCopy(G11, G00, 2*N);
				MatrixCopy(g11, G11, 2*N);
			  }
			  
			}
			
			
			if(numsep = 0)
			{
			 
			  MatrixCopy(G00, Gnn_old, 2*N);
			  MatrixCopy(GLR, GnL_old, 2*N);
			  MatrixCopy(GRL, GLn_old, 2*N);
			  MatrixCopy(G11, G00, 2*N);
			  MatrixCopy(g11, G11, 2*N);
			  
			 			  
			}
			

			//from here on assumes 2 different cells
			//work out nice way later on of generalising this...
			
			//if cell1_edge_dist=0, we are at edge, else add more cells
			if(cell1_edge_dist==0)
			{
			    MatrixCopyPart(Gnn_old, bigGnn_old, 0, 0, 0, 0, 2*N, 2*N);
			    MatrixCopyPart(G00, bigGnn_old, 0, 0, 2*N, 2*N, 2*N, 2*N);
			    MatrixCopyPart(GnL_old, bigGnn_old, 0, 0, 0, 2*N, 2*N, 2*N);
			    MatrixCopyPart(GLn_old, bigGnn_old, 0, 0, 2*N, 0, 2*N, 2*N);
			}
			
			if(cell1_edge_dist > 0)
			{
	    
			    //first cell to be added
			    ConnectSides(G00, G11, GLR, GRL, Gnn_old, GnL_old, GLn_old, GnR, GRn, V01, V10, 2*N, 2*N);
			    
			    //this relabelling is different because the n matrix has 2 cells from now on
			    MatrixCopyPart(Gnn_old, bigGnn_old, 0, 0, 0, 0, 2*N, 2*N);
			    MatrixCopyPart(G00, bigGnn_old, 0, 0, 2*N, 2*N, 2*N, 2*N);
			    MatrixCopyPart(GnL_old, bigGnn_old, 0, 0, 0, 2*N, 2*N, 2*N);
			    MatrixCopyPart(GLn_old, bigGnn_old, 0, 0, 2*N, 0, 2*N, 2*N);

			  
			  
			}
			
			
			


		      if (reim ==0)
				return  creal(Gnn_old[cell1site][cell2site]);
		      
		      if (reim ==1)
				return  cimag(Gnn_old[cell1site][cell2site]);

		FreeMatrix(G);
		FreeMatrix(GLR);
		FreeMatrix(GRL);
		FreeMatrix(G00);
		FreeMatrix(G11);
		FreeMatrix(Gnn);
		FreeMatrix(Gnn_old);
		FreeMatrix(GnL_old);
		FreeMatrix(GLn_old);

		FreeMatrix(GnR);
		FreeMatrix(GRn);
		FreeMatrix(bigGnn_old);
		FreeMatrix(bigGnL_old);
		FreeMatrix(bigGLn_old);

		FreeMatrix(bigGnR);
		FreeMatrix(bigGRn);
		FreeMatrix(g11i);
		FreeMatrix(g11);
		FreeMatrix(SL);
		FreeMatrix(SR);
		FreeMatrix(V01);
		FreeMatrix(V10);
		FreeMatrix(temp1);
		FreeMatrix(temp2);
		FreeMatrix(unit);


	}



//returns the start, end and cross term GFs for a large period system of 2^numits unit cells
	void SimpleRubioBlock (double _Complex **G11, double _Complex **GLL, double _Complex **G1L, double _Complex **GL1, double _Complex **g00, double _Complex **V01, double _Complex **V10, int numits, int dim)
	{
		double _Complex **A = createSquareMatrix(dim);
		double _Complex **B = createSquareMatrix(dim);
		double _Complex **t1 = createSquareMatrix(dim);
		double _Complex **t2 = createSquareMatrix(dim);
		double _Complex **t3 = createSquareMatrix(dim);
		double _Complex **t4 = createSquareMatrix(dim);

		double _Complex **g11old = createSquareMatrix(dim);
		double _Complex **gLLold = createSquareMatrix(dim);
		double _Complex **g1Lold = createSquareMatrix(dim);
		double _Complex **gL1old = createSquareMatrix(dim);



		double _Complex **unit= createSquareMatrix(dim);


		int i, j, k;
		for(i=0; i<dim; i++)
			unit[i][i] = 1.0;


		//1st iteration - 2 unit chains 
			MatrixMult (g00, V01, t1, dim);		//t1 = g00 V01
			MatrixMult (g00, V10, t2, dim);		//t2 = g00 V10
	
			MatrixMult (t1, t2, t3, dim);
			MatrixSubtract (unit, t3, t4, dim);
			InvertMatrixGSL(t4, t3, dim);
			MatrixMult(t3, g00, G11, dim);
	
			MatrixMult(t2, t1, t3, dim);
			MatrixSubtract (unit, t3, t4, dim);
			InvertMatrixGSL(t4, t3, dim);
			MatrixMult(t3, g00, GLL, dim);
	
			MatrixMult(t1, GLL, G1L, dim);
			MatrixMult(t2, G11, GL1, dim);

			MatrixCopy(G11, g11old, dim);
			MatrixCopy(GLL, gLLold, dim);
			MatrixCopy(G1L, g1Lold, dim);
			MatrixCopy(GL1, gL1old, dim);		


		//all remaining iterations: 2^(j+1) unit cells

		for(j=1; j<numits; j++)
		{
			
			MatrixMult(g11old, V10, t1, dim);	//t1 = g11 V_RL
			MatrixMult(gLLold, V01, t2, dim);	//t2 = gLL V_LR

			MatrixMult(t2, t1, t3, dim);
			MatrixSubtract(unit, t3, t4, dim);
			InvertMatrixGSL(t4, t3, dim);
			MatrixMult(t3, t2, t4, dim);
			MatrixMult(t4, g1Lold, A, dim);

			MatrixMult(t1, t2, t3, dim);
			MatrixSubtract(unit, t3, t4, dim);
			InvertMatrixGSL(t4, t3, dim);
			MatrixMult(t3, t1, t4, dim);
			MatrixMult(t4, gL1old, B, dim);
		

			//new GLL
			MatrixMult (gL1old, V10, t3, dim);
			MatrixMult (t3, A, t4, dim);
			MatrixAdd (gLLold, t4, GLL, dim);

			//new G1L
			MatrixMult (t1, A, t3, dim);
			MatrixAdd (g1Lold, t3, t4, dim);
			MatrixMult(V01, t4, t3, dim);
			MatrixMult(g1Lold, t3, G1L, dim);

			//new G11
			MatrixMult (g1Lold, V01, t3, dim);
			MatrixMult (t3, B, t4, dim);
			MatrixAdd (g11old, t4, G11, dim);

			//new GL1
			MatrixMult (t2, B, t3, dim);
			MatrixAdd (gL1old, t3, t4, dim);
			MatrixMult(V10, t4, t3, dim);
			MatrixMult(gL1old, t3, GL1, dim);

			MatrixCopy(G11, g11old, dim);
			MatrixCopy(GLL, gLLold, dim);
			MatrixCopy(G1L, g1Lold, dim);
			MatrixCopy(GL1, gL1old, dim);

		}


		FreeMatrix(A);
		FreeMatrix(B);
		FreeMatrix(t1);
		FreeMatrix(t2);
		FreeMatrix(t3);
		FreeMatrix(t4);

		FreeMatrix(g11old);
		FreeMatrix(gLLold);
		FreeMatrix(g1Lold);
		FreeMatrix(gL1old);
		FreeMatrix(unit);

	}





//returns the start, end and cross term GFs for a system of numits unit cells, with onsite disorder given by DisorderArray
	void SimpleBlock (double _Complex **G11, double _Complex **GLL, double _Complex **G1L, double _Complex **GL1, double _Complex **g00, double _Complex **V01, double _Complex **V10, int numits, int dim, double _Complex *DisorderArray)
	{
		
		double _Complex **t1 = createSquareMatrix(dim);
		double _Complex **t2 = createSquareMatrix(dim);
		double _Complex **t3 = createSquareMatrix(dim);
		double _Complex **t4 = createSquareMatrix(dim);
		double _Complex **Vtemp = createSquareMatrix(dim);
		double _Complex **G00 = createSquareMatrix(dim);

		double _Complex **g11old = createSquareMatrix(dim);
		double _Complex **gLLold = createSquareMatrix(dim);
		double _Complex **g1Lold = createSquareMatrix(dim);
		double _Complex **gL1old = createSquareMatrix(dim);


		int i, j, k;

		double _Complex **unit= createSquareMatrix(dim);

		for(i=0; i<dim; i++)
			unit[i][i] = 1.0;



			for(i=0; i<dim; i++)
				Vtemp[i][i] = DisorderArray[i] ;

			dyson(g00, Vtemp, G00, dim);

			MatrixCopy(G00, g11old, dim);
			MatrixCopy(G00, gLLold, dim);
			MatrixCopy(G00, g1Lold, dim);
			MatrixCopy(G00, gL1old, dim);		


		



		for(j=1; j<numits; j++)
		{
			
			for(i=0; i<dim; i++)
				Vtemp[i][i] = DisorderArray[j*dim + i] ;
			
			dyson(g00, Vtemp, G00, dim);

			MatrixMult(G00, V10, t1, dim);
			MatrixMult(t1, gLLold, t2, dim);
			MatrixMult(t2, V01, t1, dim);
			MatrixSubtract(unit, t1, t2, dim);
			InvertMatrixGSL(t2, t1, dim);
			MatrixMult(t1, G00, GLL, dim);		//GLL

			MatrixMult(GLL, V10, t1, dim);
			MatrixMult(t1, gL1old, GL1, dim);	//GL1

			MatrixMult(g1Lold, V01, t1, dim);	
			MatrixMult(t1, GL1, t2, dim);		
			MatrixAdd(g11old, t2, G11, dim);	//G11

			MatrixMult(t1, GLL, G1L, dim);		//G1L


			MatrixCopy(G11, g11old, dim);
			MatrixCopy(GLL, gLLold, dim);
			MatrixCopy(G1L, g1Lold, dim);
			MatrixCopy(GL1, gL1old, dim);	
			
						

		}


		

		FreeMatrix(t1);
		FreeMatrix(t2);
		FreeMatrix(t3);
		FreeMatrix(t4);
		FreeMatrix(Vtemp);
		FreeMatrix(G00);

		FreeMatrix(g11old);
		FreeMatrix(gLLold);
		FreeMatrix(g1Lold);
		FreeMatrix(gL1old);
		FreeMatrix(unit);

	}

//version without disorder
void VSimpleBlock (double _Complex **G11, double _Complex **GLL, double _Complex **G1L, double _Complex **GL1, double _Complex **g00, double _Complex **V01, double _Complex **V10, int numits, int dim)
	{
		
		double _Complex **t1 = createSquareMatrix(dim);
		double _Complex **t2 = createSquareMatrix(dim);
		double _Complex **t3 = createSquareMatrix(dim);
		double _Complex **t4 = createSquareMatrix(dim);

		double _Complex **g11old = createSquareMatrix(dim);
		double _Complex **gLLold = createSquareMatrix(dim);
		double _Complex **g1Lold = createSquareMatrix(dim);
		double _Complex **gL1old = createSquareMatrix(dim);


		int i, j, k;

		double _Complex **unit= createSquareMatrix(dim);

		for(i=0; i<dim; i++)
			unit[i][i] = 1.0;



			

			MatrixCopy(g00, g11old, dim);
			MatrixCopy(g00, gLLold, dim);
			MatrixCopy(g00, g1Lold, dim);
			MatrixCopy(g00, gL1old, dim);		

			MatrixCopy(g11old, G11, dim);
			MatrixCopy(g11old, GLL, dim);
			MatrixCopy(g11old, GL1, dim);
			MatrixCopy(g11old, G1L, dim);




		for(j=1; j<numits; j++)
		{
			
			
			MatrixMult(g00, V10, t1, dim);
			MatrixMult(t1, gLLold, t2, dim);
			MatrixMult(t2, V01, t1, dim);
			MatrixSubtract(unit, t1, t2, dim);
			InvertMatrixGSL(t2, t1, dim);
			MatrixMult(t1, g00, GLL, dim);		//GLL

			MatrixMult(GLL, V10, t1, dim);
			MatrixMult(t1, gL1old, GL1, dim);	//GL1

			MatrixMult(g1Lold, V01, t1, dim);	
			MatrixMult(t1, GL1, t2, dim);		
			MatrixAdd(g11old, t2, G11, dim);	//G11

			MatrixMult(t1, GLL, G1L, dim);		//G1L


			MatrixCopy(G11, g11old, dim);
			MatrixCopy(GLL, gLLold, dim);
			MatrixCopy(G1L, g1Lold, dim);
			MatrixCopy(GL1, gL1old, dim);	
			
						

		}

		
		

		FreeMatrix(t1);
		FreeMatrix(t2);
		FreeMatrix(t3);
		FreeMatrix(t4);

		FreeMatrix(g11old);
		FreeMatrix(gLLold);
		FreeMatrix(g1Lold);
		FreeMatrix(gL1old);
		FreeMatrix(unit);

	}


/*routine for connecting a block (GFs known to the system)
		- g00 is the left side SGF before the block is added 
		- g11, gLL, g1L and gL1 are the relevant GFs of the block
		- newSF is the surface GF after the block is added (i.e. GLL)
		- gnn, gn0, g0n are elements inside the old system that need updating. gnn has dimension dim2 - if dim2=0, these are not updated or used (all other matrices have dimension dim1)
			otherwise they are replaced with Gnn, and GnL, GLn - ie connections to new surface
*/

	void ConnectBlock (double _Complex **g00, double _Complex **g11, double _Complex **gLL, double _Complex **g1L, double _Complex **gL1, double _Complex **newSF, double _Complex **gnn, double _Complex **gn0, double _Complex **g0n, double _Complex **V01, double _Complex **V10, int dim1, int dim2)
	{
		
		double _Complex **A = createSquareMatrix(dim1);
		double _Complex **t1 = createSquareMatrix(dim1);
		double _Complex **t2 = createSquareMatrix(dim1);
		double _Complex **t3 = createSquareMatrix(dim1);
		double _Complex **t4 = createSquareMatrix(dim1);


		double _Complex **unit= createSquareMatrix(dim1);


		int i, j, k;
		for(i=0; i<dim1; i++)
			unit[i][i] = 1.0;


		MatrixMult(g00, V01, t1, dim1);		//t1 = g00 V01
		MatrixMult(g11, V10, t2, dim1);  	//t2 = g11 V10

		MatrixMult(t1, t2, t3, dim1);
		MatrixSubtract(unit, t3, t4, dim1);
		InvertMatrixGSL(t4, t3, dim1);
		MatrixMult(t3, t1, t4, dim1);
		MatrixMult(t4, g1L, A, dim1); 		//A = G0L

		MatrixMult(gL1, V10, t3, dim1);
		MatrixMult(t3, A, t4, dim1);
		MatrixAdd(gLL, t4, newSF, dim1);	//new surface Green function

		if(dim2 != 0)
		{
			double _Complex **B = createNonSquareMatrix(dim1, dim2);
			double _Complex **tn1 = createNonSquareMatrix(dim2, dim1);
			double _Complex **tn2 = createSquareMatrix(dim2);
			double _Complex **tn2a = createSquareMatrix(dim2);

			double _Complex **tn3 = createNonSquareMatrix(dim1, dim2);
			double _Complex **tn4 = createNonSquareMatrix(dim1, dim2);


			MatrixMult(t2, t1, t3, dim1);
			MatrixSubtract(unit, t3, t4, dim1);
			InvertMatrixGSL(t4, t3, dim1);
			MatrixMult(t3, t2, t4, dim1);
			MatrixMultNS(t4, g0n, B, dim1, dim1, dim2);	//B = G1n

			
			MatrixMultNS(gn0, V01, tn1, dim2, dim1, dim1);
			MatrixMultNS(tn1, B, tn2, dim2, dim1, dim2);
			MatrixAdd(gnn, tn2, tn2a, dim2);			//Gnn overwrites gnn - might cause trouble - test
			MatrixCopyPart(tn2a, gnn, 0, 0, 0, 0, dim2, dim2);	//gnL replaces gn0

			//MatrixAdd(gnn, tn2, gnn, dim2);			//Gnn overwrites gnn - might cause trouble - test
			
			MatrixMult(t2, A, t3, dim1);
			MatrixAdd(g1L, t3, t4, dim1);
			MatrixMult(V01, t4, t3, dim1);
			MatrixMultNS(gn0, t3, tn1, dim2, dim1, dim1);
			MatrixCopyPart(tn1, gn0, 0, 0, 0, 0, dim2, dim1);	//gnL replaces gn0

			MatrixMultNS(t1, B, tn3, dim1, dim1, dim2); 
			MatrixAddNS(g0n, tn3, tn4, dim1, dim2);
			MatrixMultNS(V10, tn4, tn3, dim1, dim1, dim2);
			MatrixMultNS(gL1, tn3, g0n, dim1, dim1, dim2);		//gLn replaces g0n
			

			FreeMatrix(B);
			FreeMatrix(tn1);
			FreeMatrix(tn2);
			FreeMatrix(tn2a);
			FreeMatrix(tn3);
			FreeMatrix(tn4);
		}

		FreeMatrix(A);
		FreeMatrix(t1);
		FreeMatrix(t2);
		FreeMatrix(t3);
		FreeMatrix(t4);
		FreeMatrix(unit);
	}


	
	
/*Connect LHS and RHS leads, and update matrix elements corresponding to sites within the Left hand lead (if dim2!=0)
	- matrices with a small g are updated, matrices with a large G zre previously zero 
	- essentially a more complex version of GFsites (dim2=0), but not wrappered	*/

void ConnectSides (double _Complex **gLL, double _Complex **gRR, double _Complex **GLR, double _Complex **GRL, double _Complex **gnn, double _Complex **gnL, double _Complex **gLn, double _Complex **GnR, double _Complex **GRn, double _Complex **VLR, double _Complex **VRL, int dim1, int dim2)
	{
		
		double _Complex **A = createSquareMatrix(dim1);
		double _Complex **t1 = createSquareMatrix(dim1);
		double _Complex **t2 = createSquareMatrix(dim1);
		double _Complex **t3 = createSquareMatrix(dim1);
		double _Complex **t4 = createSquareMatrix(dim1);

		double _Complex **GLL = createSquareMatrix(dim1);
		double _Complex **GRR = createSquareMatrix(dim1);
	

		double _Complex **unit= createSquareMatrix(dim1);


		int i, j, k;
		for(i=0; i<dim1; i++)
			unit[i][i] = 1.0;


		MatrixMult(gLL, VLR, t1, dim1);		//t1 = gLL VLR
		MatrixMult(gRR, VRL, t2, dim1);  	//t2 = gRR VRL

		MatrixMult(t1, t2, t3, dim1);
		MatrixSubtract(unit, t3, t4, dim1);
		InvertMatrixGSL(t4, t3, dim1);
		MatrixMult(t3, gLL, GLL, dim1);		//GLL

		MatrixMult(t2, GLL, GRL, dim1);		//GRL

		MatrixMult(t2, t1, t3, dim1);
		MatrixSubtract(unit, t3, t4, dim1);
		InvertMatrixGSL(t4, t3, dim1);
		MatrixMult(t3, gRR, GRR, dim1);		//GRR

		MatrixMult(t1, GRR, GLR, dim1);		//GLR

		if(dim2 != 0)
		{
			double _Complex **tn1 = createNonSquareMatrix(dim2, dim1);
			double _Complex **tn1a = createNonSquareMatrix(dim2, dim1);
			double _Complex **tn2 = createSquareMatrix(dim2);
			double _Complex **tn3 = createNonSquareMatrix(dim1, dim2);

			double _Complex **Gnn = createSquareMatrix(dim1);
			double _Complex **GnL = createSquareMatrix(dim1);
			double _Complex **GLn = createSquareMatrix(dim1);

			
			MatrixMult(GRR, VRL, t3, dim1);				
			MatrixMultNS(t3, gLn, GRn, dim1, dim1, dim2);		//GRn

			MatrixMultNS(gnL, VLR, tn1, dim2, dim1, dim1);		//tn1 = gnL VLR
			MatrixMultNS(tn1, GRR, GnR, dim2, dim1, dim1);		//GnR

			MatrixMultNS(tn1, GRL, tn1a, dim2, dim1, dim1);
			MatrixAddNS(gnL, tn1a, GnL, dim2, dim1);		//GnL

			MatrixMultNS(t1, GRn, tn3, dim1, dim1, dim2);
			MatrixAddNS(gLn, tn3, GLn, dim1, dim2);			//GLn

			MatrixMultNS(tn1, GRn, tn2, dim2, dim1, dim2);
			MatrixAdd(gnn, tn2, Gnn, dim2);				//Gnn
			
			MatrixCopyPart(GnL, gnL, 0, 0, 0, 0, dim2, dim1);
			MatrixCopyPart(GLn, gLn, 0, 0, 0, 0, dim1, dim2);
			MatrixCopy(Gnn, gnn, dim2);

			FreeMatrix(tn1);
			FreeMatrix(tn1a);
			FreeMatrix(tn2);
			FreeMatrix(tn3);
			FreeMatrix(Gnn);
			FreeMatrix(GnL);
			FreeMatrix(GLn);

		}

		MatrixCopy(GLL, gLL, dim1);
		MatrixCopy(GRR, gRR, dim1);


		FreeMatrix(A);
		FreeMatrix(t1);
		FreeMatrix(t2);
		FreeMatrix(t3);
		FreeMatrix(t4);
		FreeMatrix(GLL);
		FreeMatrix(GRR);
		FreeMatrix(unit);
	}



/*wrapper for ribbon GFs with disorder or big blocks in them */

	void RibbonGF (double _Complex En, double _Complex **G, void *p)
	{
		RibbonGF_params *params = (RibbonGF_params *)p;
	
		//parameters
			int N = (params->size);
			double hopping = (params->hopping);
			double error = (params->rubio_error);
	
	
			int numsites = (params->numsites);
			int *pos = (params->pos);
	
			int blockfunc = (params->blockfunc);
			int blocksize = (params->blocksize);

			int block2size = (params->block2size);
			double _Complex *DA = (params->DA);

			int numij = (params->numij);
			int *posij = (params->posij);


			//printf("%d	%lf	%lf	%d	%d	%d	%d	%d	%d\n", N, hopping, error, numsites, blockfunc, blocksize, block2size, numij, posij[0]);

		int count, i, j, k;
	

		//Matrix Declarations
			double _Complex **ginv = createSquareMatrix(2*N);
			double _Complex **g00 = createSquareMatrix(2*N);
			double _Complex **SL = createSquareMatrix(2*N);
			double _Complex **SLnew = createSquareMatrix(2*N);
			double _Complex **SLnew2 = createSquareMatrix(2*N);
			double _Complex **SR = createSquareMatrix(2*N);
			double _Complex **V01 = createSquareMatrix(2*N);
			double _Complex **V10 = createSquareMatrix(2*N);
		
	
			double _Complex **G11 = createSquareMatrix(2*N);
			double _Complex **GLL = createSquareMatrix(2*N);
			double _Complex **G1L = createSquareMatrix(2*N);
			double _Complex **GL1 = createSquareMatrix(2*N);

			double _Complex **GLR = createSquareMatrix(2*N);
			double _Complex **GRL = createSquareMatrix(2*N);
			
	
			double _Complex **Gnn; 
			double _Complex **GnL; 
			double _Complex **GLn; 
			double _Complex **GnR; 
			double _Complex **GRn; 

		if(numsites>0)
		{
			Gnn = createSquareMatrix(numsites);
			GnL = createNonSquareMatrix(numsites, 2*N);
			GLn = createNonSquareMatrix(2*N, numsites);
			GnR = createNonSquareMatrix(numsites, 2*N);
			GRn = createNonSquareMatrix(2*N, numsites);
		}

			//Left and Right hand side semi infinite leads
			(params->gfroutine) (ginv, V01, V10, N, hopping, En);
			InvertMatrixGSL(ginv, g00, 2*N);
			RubioSGF (SL, g00, V10, V01, 2*N, &count, error);
			//printf("#rubioL count: %d\n", count);
			RubioSGF (SR, g00, V01, V10, 2*N, &count, error);
			//printf("#rubioR count: %d\n", count);

		
			if(numsites > 0)
			{
				for(i=0; i<numsites; i++)
				{
					for(j=0; j<2*N; j++)
					{
						GnL[i][j] = SL[pos[i]][j];
						GLn[j][i] = SL[j][pos[i]];
					}
					for(j=0; j<numsites; j++)
					{
						Gnn[i][j] = SL[pos[i]][pos[j]];
					}
				}
			}

			//large, clean block
			if(blocksize > 0)
			{	
				if(blockfunc==0)
				{
					double _Complex *emptyDA = createCompArray(blocksize*2*N);
					SimpleBlock (G11, GLL, G1L, GL1, g00, V01, V10, blocksize, 2*N, emptyDA);
					free(emptyDA);
				}
	
				if(blockfunc==1)
				{
					SimpleRubioBlock (G11, GLL, G1L, GL1, g00, V01, V10, blocksize, 2*N);
				}

				
				ConnectBlock (SL, G11, GLL, G1L, GL1, SLnew, Gnn, GnL, GLn, V01, V10, 2*N, numsites);
				
			}
			else
				MatrixCopy(SL, SLnew, 2*N);



			//disordered block
			if(block2size > 0)
			{
				SimpleBlock (G11, GLL, G1L, GL1, g00, V01, V10, block2size, 2*N, DA);
				
				
				ConnectBlock (SLnew, G11, GLL, G1L, GL1, SLnew2, Gnn, GnL, GLn, V01, V10, 2*N, numsites);		
				
			}

			else
				MatrixCopy(SLnew, SLnew2, 2*N);

			
		
			//Connect Left and Right Hand Sides
			ConnectSides(SLnew2, SR, GLR, GRL, Gnn, GnL, GLn, GnR, GRn, V01, V10, 2*N, numsites);

			
			//return required Green functions in G
			MatrixCopyPart(Gnn, G, 0, 0, 0, 0, numsites, numsites);

			for(i=0; i<numij; i++)
			{
				for(j=0; j<numij; j++)
				{
					if(posij[i] < 2*N)
					{
						if(posij[j] < 2*N)
						{
							G[numsites+i][numsites+j] =  SLnew2[posij[i]][posij[j]] ;
						}
						else
						{
							G[numsites+i][numsites+j] = GLR[posij[i]][posij[j] - 2*N] ;
						}
					}
					else
					{
						if(posij[j] < 2*N)
						{
							G[numsites+i][numsites+j] = GRL[posij[i] - 2*N][posij[j]] ;
						}
						else
						{
							G[numsites+i][numsites+j] = SR[posij[i] - 2*N][posij[j] - 2*N] ;
						}
					}
				}

				for(j=0; j<numsites; j++)
				{
					if(posij[i] < 2*N)
					{
						G[numsites+i][j] = GLn[posij[i]][j] ;
						G[j][numsites+i] = GnL[j][posij[i]] ;
					}
					else
					{
						G[numsites+i][j] = GRn[posij[i] -2*N][j] ;
						G[j][numsites+i] = GnR[j][posij[i] - 2*N] ;
					}

				}

			}


			FreeMatrix(ginv);
			FreeMatrix(g00);
			FreeMatrix(SL);
			FreeMatrix(SLnew);
			FreeMatrix(SLnew2);
			FreeMatrix(SR);
			FreeMatrix(V01);
			FreeMatrix(V10);
			FreeMatrix(G11);
			FreeMatrix(GLL);
			FreeMatrix(GL1);
			FreeMatrix(G1L);
			FreeMatrix(GLR);
			FreeMatrix(GRL);

			if(numsites>0)
			{
				FreeMatrix(Gnn); 
				FreeMatrix(GnL); 
				FreeMatrix(GLn); 
				FreeMatrix(GnR); 
				FreeMatrix(GRn); 
			}
			

			//printEMatrix(G, numsites + numij);
	}



/*wrapper for ribbon GFs for lots of arbitrary sites along a ribbon, that can also have disorder 
  this is handy for coupling calculations that involve, e.g. multiple sites in more than one unit cell for adsorbed impurities*/

	void RibbonGF_selective (double _Complex En, double _Complex **G, void *p)
	{
		RibbonGFs_params *params = (RibbonGFs_params *)p;
	
		//parameters
			int N = (params->size);
			double hopping = (params->hopping);
			double error = (params->rubio_error);
	
			int numcells = (params->numcells);
			int numsites = (params->numsites);
			int *pos = (params->pos);
	
			double _Complex *DA = (params->DA);

		
		int count, i, j, k, numthiscell, numtodate;

		//Matrix Declarations
			double _Complex **ginv = createSquareMatrix(2*N);
			double _Complex **g00 = createSquareMatrix(2*N);
			double _Complex **SL = createSquareMatrix(2*N);
			double _Complex **SLnew = createSquareMatrix(2*N);
			double _Complex **SLnew2 = createSquareMatrix(2*N);
			double _Complex **SR = createSquareMatrix(2*N);
			double _Complex **V01 = createSquareMatrix(2*N);
			double _Complex **V10 = createSquareMatrix(2*N);
			double _Complex **Vtemp = createSquareMatrix(2*N);
			double _Complex **GLR = createSquareMatrix(2*N);
			double _Complex **GRL = createSquareMatrix(2*N);
			double _Complex **G00 = createSquareMatrix(2*N);
			double _Complex **G11 = createSquareMatrix(2*N);
			double _Complex **unit = createSquareMatrix(2*N);
			for(i=0; i<2*N; i++)
				unit[i][i] = 1.0;

	
			double _Complex **Gnn; 
			double _Complex **GnL; 
			double _Complex **GLn; 
			double _Complex **GnR; 
			double _Complex **GRn; 
			double _Complex **Gnn_old; 
			double _Complex **GnL_old; 
			double _Complex **GLn_old; 
			int *thiscellsites;



			//Left and Right hand side semi infinite leads
				(params->gfroutine) (ginv, V01, V10, N, hopping, En);
				InvertMatrixGSL(ginv, g00, 2*N);
				RubioSGF (SL, g00, V10, V01, 2*N, &count, error);
				RubioSGF (SR, g00, V01, V10, 2*N, &count, error);

		
			//FIRST CHAIN - SL with relevant onsites added with dyson
				for(i=0; i<2*N; i++)
					Vtemp[i][i] = DA[i] ;
				dyson(SL, Vtemp, G00, 2*N);


			//SECOND CHAIN - g00 dysoned
				for(i=0; i<2*N; i++)
					Vtemp[i][i] = DA[2*N + i] ;
				dyson(g00, Vtemp, G11, 2*N);


			//Connect first and second chains
				
				//call connect sides
				ConnectSides(G00, G11, GLR, GRL, Gnn_old, GnL_old, GLn_old, GnR, GRn, V01, V10, 2*N, 0);



			//Fill Gnn, GnL, GLn with any elements from old edge (00)
				numthiscell=0; numtodate=0;

				//count number of sites in cell to add to "n"
					for(i=0; i< numsites; i++)
					{
						if(pos[i] < 2*N)
							numthiscell++;
					}			
	

				//if num of sites in this cell > 0,  add relevant terms to Gnn, GnL, GLn
				if(numthiscell > 0)
				{
					Gnn = createSquareMatrix(numtodate + numthiscell);
					GnL = createNonSquareMatrix(numtodate + numthiscell, 2*N);
					GLn = createNonSquareMatrix(2*N, numtodate + numthiscell);
					
// 					if(numtodate > 0)
// 					{
// 						//update first numtodate elements of Gnn
// 					}

					thiscellsites = createIntArray(numthiscell);	//remember to free this later!

					for(i=0; i<numthiscell; i++)
					{
						thiscellsites[i] = pos[numtodate + i];	
					}	
	
					for(i=0; i<numthiscell; i++)
					{
						for(j=0; j<numthiscell; j++)
						{
							Gnn[i][j] = G00[thiscellsites[i]][thiscellsites[j]];
						}
						for(j=0; j<2*N; j++)
						{
							GnL[i][j] = GLR[thiscellsites[i]][j];
							GLn[j][i] = GRL[j][thiscellsites[i]];
						}

					}
					numtodate += numthiscell;

					
					//matrix jiggerypokery
					Gnn_old = createSquareMatrix(numtodate);
					GnL_old = createNonSquareMatrix(numtodate, 2*N);
					GLn_old = createNonSquareMatrix(2*N, numtodate);
					MatrixCopy(Gnn, Gnn_old, numtodate);
					MatrixCopyPart(GnL, GnL_old, 0, 0, 0, 0, numtodate, 2*N);
					MatrixCopyPart(GLn, GLn_old, 0, 0, 0, 0, 2*N, numtodate);
					MatrixCopy(G11, G00, 2*N);
					FreeMatrix(Gnn);
					FreeMatrix(GLn);
					FreeMatrix(GnL);
					free(thiscellsites);

				}



			for(j=1; j<numcells-1; j++)
			{
				  					    		    		    				      												

				//(j+2)^th cell - g00 dysoned
				  for(i=0; i<2*N; i++)
					  Vtemp[i][i] = DA[2*N*(j+1) + i] ;
				  dyson(g00, Vtemp, G11, 2*N);
				
				//matrix declarations  
				  if(numtodate>0)
				  {
				    GnR = createNonSquareMatrix(numtodate, 2*N);
				    GRn = createNonSquareMatrix(2*N, numtodate);
				  }
							
				//Connect (j+2)^th cell
				  ConnectSides(G00, G11, GLR, GRL, Gnn_old, GnL_old, GLn_old, GnR, GRn, V01, V10, 2*N, numtodate);
				  
				  

				  
				//count number of sites in previous edge (j+1) to add to "n"
				  numthiscell=0;

				  for(i=numtodate; i< numsites; i++)
				  {
					  if(pos[i] < 2*N*(j+1) && pos[i] >= 2*N*(j))
						  numthiscell++;
				  }			


				//if num of sites in this cell > 0,  add relevant terms to Gnn, GnL, GLn
				 //some of these should be updated regardless: i.e. Gnr -> Glr etc
				    //  if(numthiscell > 0)
				    //  {
  
					      
					  Gnn = createSquareMatrix(numtodate + numthiscell);
					  GnL = createNonSquareMatrix(numtodate + numthiscell, 2*N);
					  GLn = createNonSquareMatrix(2*N, numtodate + numthiscell);
						    
					  //update first numtodate elements of Gnn, GLn, GnL
					    if(numtodate > 0)
					    {
						    //update first numtodate elements of Gnn
						    MatrixCopyPart(Gnn_old, Gnn, 0, 0, 0, 0, numtodate, numtodate);
						    MatrixCopyPart(GnR, GnL, 0, 0, 0, 0, numtodate, 2*N);
						    MatrixCopyPart(GRn, GLn, 0, 0, 0, 0, 2*N, numtodate);
					    }

					    thiscellsites = createIntArray(numthiscell);	

					    //run a check to make sure this is selecting the right elements in each cell!
					    //i.e. comment out all the calculations and just check indices
					    for(i=0; i<numthiscell; i++)
					    {
						    thiscellsites[i] = pos[numtodate + i - 2*N*j];	
					    }	
					    
		    
						    for(i=0; i<numthiscell; i++)
						    {
							    for(j=0; j<numthiscell; j++)
							    {
								    Gnn[numtodate+i][numtodate+j] = G00[thiscellsites[i]][thiscellsites[j]];
							    }
							    for(j=0; j<2*N; j++)
							    {
								    GnL[numtodate+i][j] = GLR[thiscellsites[i]][j];
								    GLn[j][numtodate+i] = GRL[j][thiscellsites[i]];
							    }

						    }
						    
						//    if(numtodate>0)
						//    {
							FreeMatrix(Gnn_old);
							FreeMatrix(GLn_old);
							FreeMatrix(GnL_old);
						//    }
						    
						    numtodate += numthiscell;
						    
						    
						    //matrix jiggerypokery
						    Gnn_old = createSquareMatrix(numtodate);
						    GnL_old = createNonSquareMatrix(numtodate, 2*N);
						    GLn_old = createNonSquareMatrix(2*N, numtodate);
						    MatrixCopy(Gnn, Gnn_old, numtodate);
						    MatrixCopyPart(GnL, GnL_old, 0, 0, 0, 0, numtodate, 2*N);
						    MatrixCopyPart(GLn, GLn_old, 0, 0, 0, 0, 2*N, numtodate);
						    MatrixCopy(G11, G00, 2*N);
						    FreeMatrix(Gnn);
						    FreeMatrix(GLn);
						    FreeMatrix(GnL);
						    FreeMatrix(GRn);
						    FreeMatrix(GnR);
						    free(thiscellsites);

					//    }
				  


			}





			
		
			//Connect Left and Right Hand Sides
// 			ConnectSides(SLnew2, SR, GLR, GRL, Gnn, GnL, GLn, GnR, GRn, V01, V10, 2*N, numsites);
// 
// 			
// 			//return required Green functions in G
// 			MatrixCopyPart(Gnn, G, 0, 0, 0, 0, numsites, numsites);
// 
// 			for(i=0; i<numij; i++)
// 			{
// 				for(j=0; j<numij; j++)
// 				{
// 					if(posij[i] < 2*N)
// 					{
// 						if(posij[j] < 2*N)
// 						{
// 							G[numsites+i][numsites+j] =  SLnew2[posij[i]][posij[j]] ;
// 						}
// 						else
// 						{
// 							G[numsites+i][numsites+j] = GLR[posij[i]][posij[j] - 2*N] ;
// 						}
// 					}
// 					else
// 					{
// 						if(posij[j] < 2*N)
// 						{
// 							G[numsites+i][numsites+j] = GRL[posij[i] - 2*N][posij[j]] ;
// 						}
// 						else
// 						{
// 							G[numsites+i][numsites+j] = SR[posij[i] - 2*N][posij[j] - 2*N] ;
// 						}
// 					}
// 				}
// 
// 				for(j=0; j<numsites; j++)
// 				{
// 					if(posij[i] < 2*N)
// 					{
// 						G[numsites+i][j] = GLn[posij[i]][j] ;
// 						G[j][numsites+i] = GnL[j][posij[i]] ;
// 					}
// 					else
// 					{
// 						G[numsites+i][j] = GRn[posij[i] -2*N][j] ;
// 						G[j][numsites+i] = GnR[j][posij[i] - 2*N] ;
// 					}
// 
// 				}
// 
// 			}

		     // printf("ok1\n");

			FreeMatrix(ginv);
			FreeMatrix(g00);
			FreeMatrix(SL);
			FreeMatrix(SLnew);
			FreeMatrix(SLnew2);
			FreeMatrix(SR);
			FreeMatrix(V01);
			FreeMatrix(V10);
			FreeMatrix(Vtemp);
			FreeMatrix(G11);
			FreeMatrix(G00);

// 			FreeMatrix(GLL);
// 			FreeMatrix(GL1);
// 			FreeMatrix(G1L);
			FreeMatrix(GLR);
			FreeMatrix(GRL);
			FreeMatrix(unit);
			
			FreeMatrix(Gnn_old);
			FreeMatrix(GLn_old);
			FreeMatrix(GnL_old);

			if(numsites>0)
			{
				//FreeMatrix(Gnn); 
				//FreeMatrix(GnL); 
				//FreeMatrix(GLn); 
				//FreeMatrix(GnR); 
				//FreeMatrix(GRn); 
			}
			

			//printEMatrix(G, numsites + numij);
	}











	double getConduc (double _Complex En, void *p, double _Complex **V01, double _Complex **V10)
	{
		RibbonGF_params *params = (RibbonGF_params *)p;
		
		//parameters
			int N = (params->size);
			int numsites = (params->numsites);
			
	//	printf("%d	%d\n", N, numsites);

		int i, j;
		double _Complex **Gbig = createSquareMatrix(numsites + 4*N);
		double _Complex **G00tilde = createSquareMatrix(2*N);
		double _Complex **G11tilde = createSquareMatrix(2*N);
		double _Complex **G10tilde = createSquareMatrix(2*N);
		double _Complex **temp1 = createSquareMatrix(2*N);
		double _Complex **temp2 = createSquareMatrix(2*N);
		double _Complex **temp3 = createSquareMatrix(2*N);

		double ans;
	//	printf("ok\n");
		RibbonGF(En, Gbig, p);	
	//	printf("ok\n");


			for(i=0; i<2*N; i++)
			{
				for(j=0; j<2*N; j++)
				{

					G00tilde[i][j] = (Gbig[numsites + i][numsites+j] - conj(Gbig[numsites+j][numsites + i]) )/ (2*I) ;
					G11tilde[i][j] = (Gbig[numsites +2*N + i][numsites +2*N +j] - conj(Gbig[numsites + 2*N +j][numsites + 2*N + i]) )/ (2*I) ;
					G10tilde[i][j] = (Gbig[numsites +2*N + i][numsites +j] - conj(Gbig[numsites +j][numsites + 2*N + i]) )/ (2*I) ;


				}
			}
			

			MatrixMult(G00tilde, V01, temp1, 2*N);
			MatrixMult(temp1, G11tilde, temp2, 2*N);
			MatrixMult(temp2, V10, temp1, 2*N);

			MatrixMult(V01, G10tilde, temp2, 2*N);
			MatrixMult(temp2, V01, temp3, 2*N);
			MatrixMult(temp3, G10tilde, temp2, 2*N);

			MatrixSubtract(temp1, temp2, temp3, 2*N);

			ans = 4*creal(MatrixTrace(temp3, 2*N));

		FreeMatrix(Gbig);
		FreeMatrix(G00tilde);
		FreeMatrix(G11tilde);
		FreeMatrix(G10tilde);
		FreeMatrix(temp1);
		FreeMatrix(temp2);
		FreeMatrix(temp3);

		return ans;
/*
			printf("%lf 	%e	%e\n", realE, -cimag(MatrixTrace(bigG, 2*N))/M_PI,  4*creal(MatrixTrace(temp3, 2*N)));*/
	
	}




void zzedgedisorder (double _Complex *DisorderArray, double potential, double concentration, int length, int size, int depth, int seed)
{
	
	int numpossites = 2 * length * depth;	// 2 due to both edges
	int *sites = createIntArray (numpossites);
	int *sitelocs = createIntArray (numpossites);

	int numimps = (int) (numpossites * concentration) ;
	int numgen, i, j, temp;
	srand(time(NULL) + seed);

	for(i=0; i<2*size*length; i++)
		DisorderArray[i] = 0.0;

	for(i=0; i<numpossites; i++)
		sites[i] = i;

	//locations of edge sites in DisorderArray
	for(i=0; i< length; i++)
	{
		for(j=0; j<depth; j++)
		{
			sitelocs[i*2*depth + j] = i * 2 * size + j;				//near edge
			sitelocs[(i*2 +1)*depth + j] = 	2*(i +1)*size - (depth) + j;	//far edge
		}
	}

	for(i=0; i<numpossites; i++)
	{
		numgen = (rand() % numpossites);
		temp = sites[i];
		sites[i] = sites[numgen];
		sites[numgen] = temp;
	}
	for(i=0; i<numimps; i++)
	{
		DisorderArray[sitelocs[sites[i]]] = potential;
		//printf("%d\n", sitelocs[sites[i]]);
	}

	//for(i=0; i<numpossites; i++)
	//	printf("%d\n", sitelocs[i]);
	free(sites);
	free(sitelocs);

}


//returns a disorder array for edge disorder on ac ribbons

void acedgedisorder (double _Complex *DisorderArray, double potential, double concentration, int length, int size, int depth, int seed)
{
	
	int numpossites = 4 * length * depth;	// 2 due to both edges
	int *sites = createIntArray (numpossites);
	int *sitelocs = createIntArray (numpossites);

	int numimps = (int) (numpossites * concentration) ;
	int numgen, i, j, temp;
	srand(time(NULL) + seed);

	for(i=0; i<2*size*length; i++)
		DisorderArray[i] = 0.0;

	for(i=0; i<numpossites; i++)
		sites[i] = i;

	//locations of edge sites in DisorderArray
	for(i=0; i< length; i++)
	{
		for(j=0; j<depth; j++)
		{
			sitelocs[i*4*depth + j] = i * 2 * size + j;			//near edge 1
			sitelocs[(i*4 + 1)*depth + j] = (2*i +1)*size + j;		//near edge 2

			sitelocs[(i*4 + 2)*depth + j] = (2*i +1)*size - (depth) + j;	//far edge 1
			sitelocs[(i*4 + 3)*depth + j] = 2*(i +1)*size - (depth) + j;	//far edge 2

		}
	}

	for(i=0; i<numpossites; i++)
	{
		numgen = (rand() % numpossites);
		temp = sites[i];
		sites[i] = sites[numgen];
		sites[numgen] = temp;
	}
	for(i=0; i<numimps; i++)
	{
		DisorderArray[sitelocs[sites[i]]] = potential;
		//printf("%d\n", sitelocs[sites[i]]);
	}

	//for(i=0; i<numpossites; i++)
	//	printf("%d\n", sitelocs[i]);
	free(sites);
	free(sitelocs);

}


int **zzconnections (int size, int *numconnections)
{
	*numconnections = size;
	int **connections = createNonSquareIntMatrix(*numconnections, 2);
	int i, j;

	connections[0][0] = 0;
	connections[0][1] = 2*size + 1;
	j=1;

	for(i=1; i<2*size-3; i+=4)
	{
		connections[j][0] = i+2;
		connections[j][1] = 2*size + i + 1;
		j++;
		connections[j][0] = i+3;
		connections[j][1] = 2*size + i + 4;
		j++;
	}

	if(size % 2 == 0)
	{
		connections[j][0] = 2*size -1;
		connections[j][1] = 4*size -2;
	}

	return connections;
}


int **acconnections (int size, int *numconnections)
{
	if (size % 2 ==0)
		*numconnections =  (size/2);

	if (size % 2 ==1)
		*numconnections =  (size+1)/2 ;
	
	//printf("ok	%d\n", *numconnections);
	int **connections = createNonSquareIntMatrix(*numconnections, 2);
	int i, j;

	connections[0][0] = size;
	connections[0][1] = 2*size;
	j=1;
	
	for(i=1; i<size-1; i+=2)
	{
	//	printf("%d\n", i);
		connections[j][0] = size+i+1;
		connections[j][1] = 2*size + i + 1;
		//printf("%d	%d	%d\n", i, connections[j][0], connections[j][1]);

		j++;
	}
	//printf("ok\n");
	
	return connections;
}






