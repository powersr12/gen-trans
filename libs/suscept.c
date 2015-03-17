

#include "suscept.h"

	//returns the HF susceptibility matrix (or term) for the magnetic atoms i.e. CHI_AB , not CHI_ijkl
	void HFsusMatrix(double _Complex **HFsus, double omega, double fermi, double eta, double *zeeman, int N, double *m, double *shift, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams)
	{
		gsl_integration_workspace *wo;
		gsl_function integ;
		Iint_params para = {};
		integ.params = &para;
		wo = gsl_integration_workspace_alloc (100000);

		int i, j, k;
		double  I1r, I2r, I3r, I1i, I2i, I3i, error;
		para.omega = omega;
		para.omegaf = fermi;
		para.eta = eta;
		para.zeeman = zeeman;
		para.N = N;
		para.m = m;
		para.shift = shift;
		para.gfroutine = gfroutine;
		para.gfparams = gfparams;
		para.hubU = hubU;

		for(i=0; i<N; i++)
		{
			para.i1 = i;
			para.j1 = i;

			for(j=0; j<N; j++)
			{

				para.l1 = j;
				para.k1 = j;

				integ.function = &I1int;
				para.reim = 0;
				gsl_integration_qag (&integ, 0.0, 1.0, 1e-10, 0.001, 100000, 3, wo, &I1r, &error);
				//printf("r1\n");
				para.reim = 1;
				gsl_integration_qag (&integ, 0.0, 1.0, 1e-10, 0.001, 100000, 3, wo, &I1i, &error);
				//printf("i1\n");

				integ.function = &I2int;
				para.reim = 0;
				gsl_integration_qag (&integ, 0.0, 1.0, 1e-10, 0.001, 100000, 3, wo, &I2r, &error);
				//printf("r2\n");

				para.reim = 1;
				gsl_integration_qag (&integ, 0.0, 1.0, 1e-10, 0.001, 100000, 3, wo, &I2i, &error);
				//printf("i2\n");


				integ.function = &I3int;
				para.reim = 0;
				gsl_integration_qag (&integ, para.omegaf - omega, para.omegaf, 1e-10, 0.001, 100000, 3, wo, &I3r, &error);
				//printf("r3\n");

				para.reim = 1;
				gsl_integration_qag (&integ, para.omegaf - omega, para.omegaf, 1e-10, 0.001, 100000, 3, wo, &I3i, &error);
				//printf("i3\n");

				HFsus[i][j] = I1r + I2r + I3r + (I1i + I2i + I3i )*I;

			}
		}





		gsl_integration_workspace_free (wo);

	}


//returns individual HF susceptibility elements chi_ijkl 
double _Complex HFsusElement(double omega, double fermi, double eta, double *zeeman, int N, double *m, double *shift, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, int i1, int j1, int k1, int l1)
	{
		gsl_integration_workspace *wo;
		gsl_function integ;
		Iint_params para = {};
		integ.params = &para;
		wo = gsl_integration_workspace_alloc (10000);

		double  I1r, I2r, I3r, I1i, I2i, I3i, error;
		para.omega = omega;
		para.omegaf = fermi;
		para.eta = eta;
		para.zeeman = zeeman;
		para.N = N;
		para.m = m;
		para.shift = shift;
		para.gfroutine = gfroutine;
		para.gfparams = gfparams;
		para.hubU = hubU;

		
		para.i1 = i1;
		para.j1 = j1;

		para.l1 = k1;
		para.k1 = l1;

		integ.function = &I1int;
		para.reim = 0;
		gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I1r, &error);
		para.reim = 1;
		gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I1i, &error);

		integ.function = &I2int;
		para.reim = 0;
		gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I2r, &error);
		para.reim = 1;
		gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I2i, &error);

		integ.function = &I3int;
		para.reim = 0;
		gsl_integration_qags (&integ, para.omegaf - omega, para.omegaf, 0.00, 0.01, 10000, wo, &I3r, &error);
		para.reim = 1;
		gsl_integration_qags (&integ, para.omegaf - omega, para.omegaf, 0.00, 0.01, 10000, wo, &I3i, &error);


		gsl_integration_workspace_free (wo);
		return I1r + I2r + I3r + (I1i + I2i + I3i )*I;


	}










//returns the HF susceptibility matrix (or term) for the magnetic atoms i.e. CHI_AB , not CHI_ijkl
	void HFsusMatrix_comp(double _Complex **HFsus, double omega, double fermi, double eta, double *zeeman, int N, double *m, double *shift, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, double _Complex **I1, double _Complex **I2, double _Complex **I3 )
	{
		gsl_integration_workspace *wo;
		gsl_function integ;
		Iint_params para = {};
		integ.params = &para;
		wo = gsl_integration_workspace_alloc (100000);

		int i, j, k;
		double  I1r, I2r, I3r, I1i, I2i, I3i, error;
		para.omega = omega;
		para.omegaf = fermi;
		para.eta = eta;
		para.zeeman = zeeman;
		para.N = N;
		para.m = m;
		para.shift = shift;
		para.gfroutine = gfroutine;
		para.gfparams = gfparams;
		para.hubU = hubU;


		for(i=0; i<N; i++)
		{
			para.i1 = i;
			para.j1 = i;

			for(j=0; j<N; j++)
			{
				para.l1 = j;
				para.k1 = j;
				
				integ.function = &I1int;
				para.reim = 0;
				gsl_integration_qag (&integ, 0.0, 1.0, 1e-10, 0.001, 100000, 3, wo, &I1r, &error);
				//printf("r1\n");
				para.reim = 1;
				gsl_integration_qag (&integ, 0.0, 1.0, 1e-10, 0.001, 100000, 3, wo, &I1i, &error);
				I1[i][j] = I1r +I1i*I;
				//printf("i1\n");

				integ.function = &I2int;
				para.reim = 0;
				gsl_integration_qag (&integ, 0.0, 1.0, 1e-10, 0.001, 100000, 3, wo, &I2r, &error);
				//printf("r2\n");

				para.reim = 1;
				gsl_integration_qag (&integ, 0.0, 1.0, 1e-10, 0.001, 100000, 3, wo, &I2i, &error);
				I2[i][j] = I2r +I2i*I;
				//printf("i2\n");


				integ.function = &I3int;
				para.reim = 0;
				gsl_integration_qag (&integ, para.omegaf - omega, para.omegaf, 1e-10, 0.001, 100000, 3, wo, &I3r, &error);
				//printf("r3\n");

				para.reim = 1;
				gsl_integration_qag (&integ, para.omegaf - omega, para.omegaf, 1e-10, 0.001, 100000, 3, wo, &I3i, &error);
				I3[i][j] = I3r +I3i*I;
				//printf("i3\n");


				HFsus[i][j] = I1r + I2r + I3r + (I1i + I2i + I3i )*I;

			}
		}





		gsl_integration_workspace_free (wo);

	}



	//returns the integrand for the I1 integration
	double I1int (double intvar, void  *p)
	{
		Iint_params *params = (Iint_params *)p;
		double omegaf = (params->omegaf);
		double omega = (params->omega);
		double eta = (params->eta);
		int i1 = (params->i1);
		int j1 = (params->j1);
		int k1 = (params->k1);
		int l1 = (params->l1);
		int reim = (params->reim);
		int N = (params->N);
		double y = (1 + eta - intvar) / intvar;

		double _Complex **Gup= createSquareMatrix(N);
		double _Complex **Gdown= createSquareMatrix(N);
		double _Complex ans;


		MagGreen(omegaf + I*y, omegaf + omega + I*y, Gup, Gdown, params->m, params->shift, params->gfroutine, params->gfparams, params->hubU, N, params->zeeman);

		ans = (1 + eta) * Gup[l1][i1] * Gdown[j1][k1] / (2 * M_PI * intvar * intvar) ;
		//printf("%e	%e	%e\n", intvar, creal(ans), cimag(ans));

		
		FreeMatrix(Gup);
		FreeMatrix(Gdown);

		if(reim==0)
			return creal(ans);
		if(reim==1)
			return cimag(ans);

	}

	//returns the integrand for the I2 integration
	double I2int (double intvar, void  *p)
	{
		Iint_params *params = (Iint_params *)p;
		double omegaf = (params->omegaf);
		double omega = (params->omega);
		double eta = (params->eta);
		int i1 = (params->i1);
		int j1 = (params->j1);
		int k1 = (params->k1);
		int l1 = (params->l1);
		int reim = (params->reim);
		int N = (params->N);
		double y = (1 + eta - intvar) / intvar;

		double _Complex **Gup= createSquareMatrix(N);
		double _Complex **Gdown= createSquareMatrix(N);
		double _Complex ans;

		MagGreen(omegaf - omega + I*y, omegaf + I*y, Gup, Gdown, params->m, params->shift, params->gfroutine, params->gfparams, params->hubU, N, params->zeeman);

		ans = (1 + eta) * conj(Gdown[k1][j1]) * conj(Gup[i1][l1]) / (2 * M_PI * intvar * intvar) ;
		//printf("%e	%e	%e\n", intvar, creal(ans), cimag(ans));

		//printf("%e	%e	%e\n", intvar, creal(ans), cimag(ans));

		FreeMatrix(Gup);
		FreeMatrix(Gdown);

		if(reim==0)
			return creal(ans);
		if(reim==1)
			return cimag(ans);

		
	}

	//returns the integrand for the I3 integration
	double I3int (double intvar, void  *p)
	{
		Iint_params *params = (Iint_params *)p;
		double omegaf = (params->omegaf);
		double omega = (params->omega);
		double eta = (params->eta);
		int i1 = (params->i1);
		int j1 = (params->j1);
		int k1 = (params->k1);
		int l1 = (params->l1);
		int reim = (params->reim);
		int N = (params->N);
		double y = (1 + eta - intvar) / intvar;

		double _Complex **Gup= createSquareMatrix(N);
		double _Complex **Gdown= createSquareMatrix(N);
		double _Complex ans;

		MagGreen(intvar + eta*I, intvar + omega  + eta*I, Gup, Gdown, params->m, params->shift, params->gfroutine, params->gfparams, params->hubU, N, params->zeeman);

	//	ans = (1 + eta) * conj(Gdown[k1][j1]) * conj(Gup[i1][l1]) / (2 * M_PI * intvar * intvar) ;



		ans = (-I/(2*M_PI)) * conj(Gup[i1][l1]) * Gdown[j1][k1] ;
		//printf("%.40e	%e	%e	%e	%e	%e	%e\n", intvar, creal(ans), cimag(ans), creal(conj(Gup[i1][l1])), cimag(conj(Gup[i1][l1])), creal(Gdown[j1][k1]), cimag(Gdown[j1][k1]) );

		FreeMatrix(Gup);
		FreeMatrix(Gdown);

		//printf("%e	%e	%e\n", intvar, creal(ans), cimag(ans));
		if(reim==0)
			return creal(ans);
		if(reim==1)
			return cimag(ans);

	}


//PRISTINE VERSIONS OF THE GREEN FUNCTIONS

void HFsusMatrix_pristine(double _Complex **HFsus, double omega, double fermi, double eta, double *zeeman, int N, double *m, double *shift, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams)
{
		gsl_integration_workspace *wo;
		gsl_function integ;
		Iint_params para = {};
		integ.params = &para;
		wo = gsl_integration_workspace_alloc (10000);

		int i, j, k;
		double  I1r, I2r, I3r, I1i, I2i, I3i, error;
		para.omega = omega;
		para.omegaf = fermi;
		para.eta = eta;
		para.zeeman = zeeman;
		para.N = N;
		para.m = m;
		para.shift = shift;
		para.gfroutine = gfroutine;
		para.gfparams = gfparams;
		para.hubU = hubU;

		for(i=0; i<N; i++)
		{
			para.i1 = i;
			para.j1 = i;

			for(j=0; j<N; j++)
			{
				para.l1 = j;
				para.k1 = j;

				integ.function = &I1int_pristine;
				para.reim = 0;
				gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I1r, &error);
				para.reim = 1;
				gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I1i, &error);

				integ.function = &I2int_pristine;
				para.reim = 0;
				gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I2r, &error);
				para.reim = 1;
				gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I2i, &error);

				integ.function = &I3int_pristine;
				para.reim = 0;
				gsl_integration_qags (&integ, para.omegaf - omega, para.omegaf, 0.00, 0.01, 10000, wo, &I3r, &error);
				para.reim = 1;
				gsl_integration_qags (&integ, para.omegaf - omega, para.omegaf, 0.00, 0.01, 10000, wo, &I3i, &error);


				HFsus[i][j] = I1r + I2r + I3r + (I1i + I2i + I3i )*I;

			}
		}

		gsl_integration_workspace_free (wo);

}

void HFsusMatrix_pristine_comp(double _Complex **HFsus, double omega, double fermi, double eta, double *zeeman, int N, double *m, double *shift, double hubU, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, double _Complex **I1, double _Complex **I2, double _Complex **I3)
{
		gsl_integration_workspace *wo;
		gsl_function integ;
		Iint_params para = {};
		integ.params = &para;
		wo = gsl_integration_workspace_alloc (10000);

		int i, j, k;
		double  I1r, I2r, I3r, I1i, I2i, I3i, error;
		para.omega = omega;
		para.omegaf = fermi;
		para.eta = eta;
		para.zeeman = zeeman;
		para.N = N;
		para.m = m;
		para.shift = shift;
		para.gfroutine = gfroutine;
		para.gfparams = gfparams;
		para.hubU = hubU;

		for(i=0; i<N; i++)
		{
			para.i1 = i;
			para.j1 = i;

			for(j=0; j<N; j++)
			{
				para.l1 = j;
				para.k1 = j;

				integ.function = &I1int_pristine;
				para.reim = 0;
				gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I1r, &error);
				para.reim = 1;
				gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I1i, &error);
				I1[i][j] = I1r +I1i*I;

				integ.function = &I2int_pristine;
				para.reim = 0;
				gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I2r, &error);
				para.reim = 1;
				gsl_integration_qags (&integ, 0.0, 1.0, 0.00, 0.01, 10000, wo, &I2i, &error);
				I2[i][j] = I2r +I2i*I;

				integ.function = &I3int_pristine;
				para.reim = 0;
				gsl_integration_qags (&integ, para.omegaf - omega, para.omegaf, 0.00, 0.001, 10000, wo, &I3r, &error);
				//printf("\n");
				para.reim = 1;
				gsl_integration_qags (&integ, para.omegaf - omega, para.omegaf, 0.00, 0.001, 10000, wo, &I3i, &error);
				I3[i][j] = I3r +I3i*I;
				//printf("\n");



				HFsus[i][j] = I1r + I2r + I3r + (I1i + I2i + I3i )*I;

			}
		}

		gsl_integration_workspace_free (wo);

}

	double I1int_pristine (double intvar, void  *p)
	{
		Iint_params *params = (Iint_params *)p;
		double omegaf = (params->omegaf);
		double omega = (params->omega);
		double eta = (params->eta);
		int i1 = (params->i1);
		int j1 = (params->j1);
		int k1 = (params->k1);
		int l1 = (params->l1);
		int reim = (params->reim);
		int N = (params->N);
		double y = (1 + eta - intvar) / intvar;

		double _Complex **Gup= createSquareMatrix(N);
		double _Complex **Gdown= createSquareMatrix(N);
		double _Complex ans;

		NonMagGreen(omegaf + I*y, omegaf + omega + I*y, Gup, Gdown, params->m, params->shift, params->gfroutine, params->gfparams, params->hubU, N, params->zeeman);

		ans = (1 + eta) * Gup[l1][i1] * Gdown[j1][k1] / (2 * M_PI * intvar * intvar) ;

		
		//printf("%e	%e	%e\n", intvar, creal(ans), cimag(ans));
		FreeMatrix(Gup);
		FreeMatrix(Gdown);

		if(reim==0)
			return creal(ans);
		if(reim==1)
			return cimag(ans);

	}

double I2int_pristine (double intvar, void  *p)
	{
		Iint_params *params = (Iint_params *)p;
		double omegaf = (params->omegaf);
		double omega = (params->omega);
		double eta = (params->eta);
		int i1 = (params->i1);
		int j1 = (params->j1);
		int k1 = (params->k1);
		int l1 = (params->l1);
		int reim = (params->reim);
		int N = (params->N);
		double y = (1 + eta - intvar) / intvar;

		double _Complex **Gup= createSquareMatrix(N);
		double _Complex **Gdown= createSquareMatrix(N);
		double _Complex ans;

		NonMagGreen(omegaf - omega + I*y, omegaf + I*y, Gup, Gdown, params->m, params->shift, params->gfroutine, params->gfparams, params->hubU, N, params->zeeman);

		ans = (1 + eta) * conj(Gdown[k1][j1]) * conj(Gup[i1][l1]) / (2 * M_PI * intvar * intvar) ;
		//printf("%e	%e	%e\n", intvar, creal(ans), cimag(ans));
		FreeMatrix(Gup);
		FreeMatrix(Gdown);

		if(reim==0)
			return creal(ans);
		if(reim==1)
			return cimag(ans);

		
	}

	double I3int_pristine (double intvar, void  *p)
	{
		Iint_params *params = (Iint_params *)p;
		double omegaf = (params->omegaf);
		double omega = (params->omega);
		double eta = (params->eta);
		int i1 = (params->i1);
		int j1 = (params->j1);
		int k1 = (params->k1);
		int l1 = (params->l1);
		int reim = (params->reim);
		int N = (params->N);
		double y = (1 + eta - intvar) / intvar;

		double _Complex **Gup= createSquareMatrix(N);
		double _Complex **Gdown= createSquareMatrix(N);
		double _Complex ans;

		NonMagGreen(intvar + eta*I, intvar + omega  + eta*I, Gup, Gdown, params->m, params->shift, params->gfroutine, params->gfparams, params->hubU, N, params->zeeman);

		//ans = (1 + eta) * conj(Gdown[k1][j1]) * conj(Gup[i1][l1]) / (2 * M_PI * intvar * intvar) ;



		ans = (-I/(2*M_PI)) * conj(Gup[i1][l1]) * Gdown[j1][k1] ;
		if(i1==1 && k1==0)
		//printf("%.40e	%e	%e	%e	%e	%e	%e\n", intvar, creal(ans), cimag(ans), creal(conj(Gup[i1][l1])), cimag(conj(Gup[i1][l1])), creal(Gdown[j1][k1]), cimag(Gdown[j1][k1]) );

		FreeMatrix(Gup);
		FreeMatrix(Gdown);

		if(reim==0)
			return creal(ans);
		if(reim==1)
			return cimag(ans);

	}

