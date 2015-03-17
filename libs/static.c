
#include "static.h"



	double GetCoupling(double *m, double *shift, double *zeeman, double _Complex En, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams)
	{
		intcoupling_params para = {};
		para.realenergy=creal(En);
		para.m = m;
		para.shift=shift;
		para.zeeman=zeeman;
		para.hubU = hubU;
		para.gfroutine =gfroutine;
		para.gfparams=gfparams;
		double ans, error1;

		gsl_function integ;
		integ.function = &intcoupling;
		integ.params = &para;

		gsl_integration_workspace *w;
		w = gsl_integration_workspace_alloc (100000);
				//printf("intstart\n");

		gsl_integration_qagiu (&integ, cimag(En), 0.0000, 0.000001, 100000, w, &ans, &error1);

		//printf("intend\n");

		gsl_integration_workspace_free (w);
		return ans;

	}

double intcoupling (double impart, void *p)
	{
		intcoupling_params *params = (intcoupling_params *)p;
		double realen = (params->realenergy);
		double *m = (params->m);
		double *shift = (params->shift);
		double *zeeman = (params->zeeman);
		double hubU = (params->hubU);
		


		double _Complex func;
		double _Complex **g = createSquareMatrix(2);
		double _Complex **G = createSquareMatrix(2);
		double _Complex **Vbig = createSquareMatrix(2);
		double _Complex En = realen + impart*I;
		double _Complex Gup, Gdown;
		double Vex = hubU * m[0];
		double ans;

		//pristine GF matrix
		(params->gfroutine)(En, g, params->gfparams);
			//printf("ok\n");

			//Electron-electron terms
			Vbig[0][0] = -(hubU/2)*m[0] + shift[0] -zeeman[0]/2; 
		//printf("ok\n");

			Vbig[1][1] = -(hubU/2)*m[1] + shift[1] -zeeman[1]/2;

					//printf("ok\n");

			//dyson(g, Vbig, G, 2);	
			//Gup = G[0][1];
			Gup = g[0][1];

			//printEMatrix(g,2);
			//Electron-electron terms
			Vbig[0][0] = (hubU/2)*m[0] + shift[0] +zeeman[0]/2; 
			Vbig[1][1] = (hubU/2)*m[1] + shift[1] +zeeman[1]/2;

			
			//dyson(g, Vbig, G, 2);	
			//Gdown = G[1][0];
			Gdown = g[1][0];

			//func = 2*Vex*Vex*Gup*Gdown;
			//ans = creal((1.0/M_PI)*(func/(1+func)) );
			ans =  (1.0/M_PI) * log ( cabs (1.0 +  Vex * Vex * Gup * Gdown) )  ;
			//printf("%lf\n", shift);

			FreeMatrix(g);
			FreeMatrix(G);
			FreeMatrix(Vbig);
			return ans;
	}



double GetEnsembleDiff(double *m, double *shift, double _Complex En, double hubU, int Ntotal, int Nblack, int Nwhite, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams)
	{
		intEnsemble_params para = {};
		para.realenergy=creal(En);
		para.m = m;
		para.shift=shift;
		para.hubU = hubU;
		para.Ntotal = Ntotal;
		para.Nblack = Nblack;
		para.Nwhite = Nwhite;
		para.gfroutine =gfroutine;
		para.gfparams=gfparams;
		double ans, error1;

		gsl_function integ;
		integ.function = &intEnsembleDiff;
		integ.params = &para;

		gsl_integration_workspace *w;
		w = gsl_integration_workspace_alloc (100000);
				//printf("intstart\n");

		gsl_integration_qagiu (&integ, cimag(En), 1.0e-10, 1.0e-2, 100000, w, &ans, &error1);

		//printf("intend\n");

		gsl_integration_workspace_free (w);
		return ans;

	}




double intEnsembleDiff (double impart, void *p)
	{
		intEnsemble_params *params = (intEnsemble_params *)p;
		double realen = (params->realenergy);
		double *m = (params->m);
		double *shift = (params->shift);
		double hubU = (params->hubU);
		int Ntotal = (params->Ntotal);
		int Nblack = (params->Nblack);
		int Nwhite = (params->Nwhite);


		double _Complex func;
		double _Complex **g = createSquareMatrix(Ntotal);
		double _Complex **Gup = createSquareMatrix(Ntotal);
		double _Complex **Gdown = createSquareMatrix(Ntotal);
		double _Complex **Gwhite = createSquareMatrix(2*Nwhite);
		double _Complex **Vud = createSquareMatrix(Ntotal);
		double _Complex **unit = createSquareMatrix(2*Nwhite);
		double _Complex **temp1 = createSquareMatrix(2*Nwhite);
		double _Complex **temp2 = createSquareMatrix(2*Nwhite);
		double _Complex **Vbig = createSquareMatrix(2*Nwhite);


		double _Complex **G = createSquareMatrix(2);
		double _Complex En = realen + impart*I;
		double ans;
		int i, j, k;
		
		  for(i=0; i<2*Nwhite; i++)
		    unit[i][i] = 1.0;

		//pristine GF matrix 
		      (params->gfroutine)(En, g, params->gfparams);
			//printf("ok\n");

			
		//Green functions for upspin electrons
		      for(i=0; i< Ntotal; i++)
		      {
			Vud[i][i] = -(hubU/2)*m[i] + shift[i] /2; 
		      }
		      dyson(g, Vud, Gup, Ntotal);	
			
		  
		//Green functions for downspin electrons
		      for(i=0; i< Ntotal; i++)
		      {
			Vud[i][i] = (hubU/2)*m[i] + shift[i] /2; 
		      }
		      dyson(g, Vud, Gdown, Ntotal);
		  
		  
		//Green functions for white sites only - both up and down spins
		      MatrixCopyPart(Gup, Gwhite, Nblack, Nblack, 0, 0, Nwhite, Nwhite);
		      MatrixCopyPart(Gdown, Gwhite, Nblack, Nblack, Nwhite, Nwhite, Nwhite, Nwhite);
		      
		      

		//Generate matrix whose determinant appears in the Lloyd formula expression
		// (I - g V)   - where I is the unit matrix, g is the block-diagonal white sites matrix for both spin orientations
		//             - and V is the perturbation matrix changing the potential seen by electrons when the spin direction at the white sites is flipped to be opposite to those at black sites
		      for(i=0; i<Nwhite; i++)
			Vbig[i][i] = hubU*m[i];
		      
		      for(i=Nwhite; i<2*Nwhite; i++)
			Vbig[i][i] =  -hubU*m[i];
		
		
		      MatrixMult(Gwhite, Vbig, temp1, 2*Nwhite);
		      MatrixSubtract(unit, temp1, temp2, 2*Nwhite);
		      
		     // printEMatrix(Gwhite, 2*Nwhite);
		      
		      ans =  (1.0/M_PI) * log ( cabs (GetDeterminantGSL( temp2, 2*Nwhite)    )     )  ;

		      
		      
			FreeMatrix(g);
			FreeMatrix(G);
			FreeMatrix(Gup);
			FreeMatrix(Gdown);
			FreeMatrix(Gwhite);
			FreeMatrix(Vud);
			FreeMatrix(unit);
			FreeMatrix(temp1);
			FreeMatrix(temp2);

			FreeMatrix(Vbig);
			return ans;
	}
