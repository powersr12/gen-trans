
#include "hubbard.h"

void newHub (double *m, double *n, double *shift, void *hub_params, void (*gfroutine2) (double _Complex, double _Complex **, void *), void *gfparams)
{
	
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	gsl_function F;

	int i, status;
	double msumdiff=1.0, msum, x_lo, x_hi, root=0.0;


	nroots_params2 para = {};

	F.function = &nroots_new;
	F.params = &para;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);

	newHub_params *params =  (newHub_params *)hub_params;
	int N = (params->N);
	double mdthresh = (params->msumdiff_thresh);
	double max_shift = (params->max_shift);
	double min_shift = (params->min_shift);


	double *shiftnew = createDoubleArray(N);
	double *mnew = createDoubleArray(N);


	para.N = (params->N);
	para.gfroutine2 = gfroutine2;
	para.gfparams = gfparams;
	para.fermi = (params->fermi);
	para.imagEn = (params->imagEn);
	para.n = n;
	para.m = m;
	para.mnew = mnew;
	para.hubU = (params->hubU);
	para.zeeman = (params->zeeman);
	para.shift=shift;
	


	for(i=0; i<N; i++)
	{
		m[i] = (params->init_m);
	}

	//fill parameters as necessary

	//repeat until sum of absolute values of changes in magnetic moments below convergence threshold
	while (msumdiff > mdthresh)
	{
		msumdiff=0.0;
		msum=0.0;


		//loop over magnetic atoms
		for(i=0; i<N; i++)
		{
			para.i = i;

			//set up rootfinder
			x_lo = min_shift;
			x_hi = max_shift;
			gsl_root_fsolver_set (s, &F, x_lo, x_hi);
			status = GSL_CONTINUE;

			//run rootfinder until convergence
			do
			{
				status = gsl_root_fsolver_iterate (s);
				root = gsl_root_fsolver_root(s);
				x_lo = gsl_root_fsolver_x_lower (s);
				x_hi = gsl_root_fsolver_x_upper (s);
				status = gsl_root_test_interval (x_lo, x_hi, 0.00001, 0.00001);
			//printf("#	%d	%e	%e	%e	%e	%e	%lf\n", i, x_lo, x_hi, root, x_hi - x_lo, m[i], shift[i]);

			}
			while (status == GSL_CONTINUE);

			//printf("#	%d	%lf	%lf	%lf	%lf\n", i, root,  shift[i], m[i], mnew[i]);



			//store new value of band-centre shift
			if(status == GSL_SUCCESS)
			{
				shiftnew[i] = root;
			}
			
			//contribution of this atom to absolute change in moments
			msumdiff += fabs(mnew[i] - m[i]);

		}

		//update values of moment and bandcentre with those calculated in this iteration
		for(i=0; i<N; i++)
		{
			printf("#	%lf\t", mnew[i]);
			m[i] = mnew[i] ;
			shift[i] = shiftnew[i];
		}	
		printf("\n");

	}

	free(mnew);
	free(shiftnew);
	gsl_root_fsolver_free (s);

}

double nroots_new (double shifti, void *p)
{
	nroots_params2 *params = (nroots_params2 *)p;
	intdos_params2 para = {};

	

	double nup, ndown, error1, error2;
	int j;

	//DOS integrator - initialisation
	gsl_function integ;
	integ.function = &intdos_new;
	integ.params = &para;

	//parameters needed in this routine
	double *shift = (params->shift);
	int N = (params->N);
	double hubU = (params->hubU);
	double *m = (params->m);
	double *mnew = (params->mnew);
	double *n = (params->n);
	int i = (params->i);
	double *zeeman = (params->zeeman);

	//perturbation matrices for up and down spin electrons
	double _Complex **Vup = createSquareMatrix(N);
	double _Complex **Vdown = createSquareMatrix(N);

	//parameters passed to DOS integrator
	para.i = i;
	para.N = N;
	para.fermi = (params->fermi);
	para.imagEn = (params->imagEn);
	para.gfroutine2 = (params->gfroutine2);
	para.gfparams = (params->gfparams);
	


	//set-up perturbation potential
		//site being investigated has potential generated from current "guess" bandcentre shift
		//other sites use bandcentre shifts from previous iteration of loop in newHub
		//this method ensures a more symmetric end result

	for(j=0; j< N; j++)
	{
		if(j!=i)
		{
			Vup[j][j] =  -(hubU/2)*(m[j]) + shift[j] -zeeman[i]/2;  
			Vdown[j][j] = (hubU/2)*(m[j]) + shift[j] +zeeman[i]/2; 
		}
		else
		{
			Vup[j][j] =  -(hubU/2)*(m[j]) + shifti -zeeman[i]/2;  
			Vdown[j][j] = (hubU/2)*(m[j]) + shifti +zeeman[i]/2; 
		}
	}

	gsl_integration_workspace *w;
	w = gsl_integration_workspace_alloc (100000);
	para.V = Vup;
	gsl_integration_qags (&integ, 0.0, 1.0, 0.00000, 0.0001, 100000,  w, &nup, &error1);
	nup += 0.5;

	para.V = Vdown;
	gsl_integration_qags (&integ, 0.0, 1.0, 0.00000, 0.0001, 100000,  w, &ndown, &error2);
	ndown += 0.5;

	gsl_integration_workspace_free (w);


	mnew[i] = nup - ndown;

	FreeMatrix(Vup);
	FreeMatrix(Vdown);

	return n[i] - nup - ndown;
}



	double intdos_new(double impart, void *p)
	{
		intdos_params2 *params = (intdos_params2 *)p;

		double ans;

		int i = (params->i);
		int N = (params -> N);
		double imagEn = (params->imagEn);
		double fermi = (params->fermi);
		double _Complex **V = (params->V);
		void (*gfroutine2)(double _Complex, double _Complex **, void *) = (params->gfroutine2);
		void *gfparams = (params->gfparams);

		double _Complex **G = createSquareMatrix(N);
		double _Complex **Gpert = createSquareMatrix(N);

		double _Complex En = fermi + ((1 + imagEn - impart)/impart)*I;

		(*gfroutine2) (En, G, gfparams);
		
		dyson(G, V, Gpert, N);
		ans = ((1 + imagEn)/M_PI) * (creal(Gpert[i][i])) / (impart*impart) ;

		FreeMatrix(G);
		FreeMatrix(Gpert);

		return ans;
	}
