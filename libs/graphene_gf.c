

#include "graphene_gf.h"
#include "matrices.h"


	//returns the (up AND down spin) graphene GF matrices required for the integrands of I1, I2 and I3.
	void MagGreen(double _Complex upenergy, double _Complex downenergy, double _Complex **Gup, double _Complex **Gdown, double *m, double *shift, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, double hubU, int N, double *zeemanterm)
	{
		
		//Memory declarations
			double _Complex **Gun= createSquareMatrix(N);
			double _Complex **Vup= createSquareMatrix(N);
			double _Complex **Vdown= createSquareMatrix(N);
			int count, i, j;

		//define Vup and Vdown
			for(i=0; i<N; i++)
			{
				Vup[i][i] = -(hubU/2)*(m[i]) + shift[i] -zeemanterm[i]/2; 
				Vdown[i][i] = (hubU/2)*(m[i]) + shift[i] +zeemanterm[i]/2;
			}
		//printf("ok\n");

		//Gup
			(*gfroutine)(upenergy, Gun, gfparams);			
			dyson(Gun, Vup, Gup, N);

		//Gdown
			(*gfroutine)(downenergy, Gun, gfparams);			
			dyson(Gun, Vdown, Gdown, N);

		//Free up memory
			FreeMatrix(Gun);
			FreeMatrix(Vup);
			FreeMatrix(Vdown);
	

	}

	//returns the (up AND down spin) graphene GF matrices required for the integrands of I1, I2 and I3.
	//This version returns the pristine Green functions
	void NonMagGreen(double _Complex upenergy, double _Complex downenergy, double _Complex **Gup, double _Complex **Gdown, double *m, double *shift, void (*gfroutine)(double _Complex , double _Complex **, void *), void *gfparams, double hubU, int N, double *zeemanterm)
	{
		
		//Memory declarations
			double _Complex **Gun= createSquareMatrix(N);
			double _Complex **Vup= createSquareMatrix(N);
			double _Complex **Vdown= createSquareMatrix(N);
			int count, i, j;

		//define Vup and Vdown
			for(i=0; i<N; i++)
			{
				Vup[i][i] = -(hubU/2)*(m[i]) + shift[i] -zeemanterm[i]/2; 
				Vdown[i][i] = (hubU/2)*(m[i]) + shift[i] +zeemanterm[i]/2;
			}
		//printf("ok\n");

		//Gup
			(*gfroutine)(upenergy, Gun, gfparams);	
			MatrixCopy(Gun, Gup, N);		
			//dyson(Gun, Vup, Gup, N);

		//Gdown
			(*gfroutine)(downenergy, Gun, gfparams);
			MatrixCopy(Gun, Gdown, N);		
			//dyson(Gun, Vdown, Gdown, N);

		//Free up memory
			FreeMatrix(Gun);
			FreeMatrix(Vup);
			FreeMatrix(Vdown);
	

	}



//nice wrapper routine for a matrix of Green functions for sites whose positions in form (a1, a2, sublattice) are in pos[][]
	void grapheneGFmatrix( double _Complex En, double _Complex **G, void *p)
	{
		gfmat_p *params = (gfmat_p *)p;
		int i, j;

		int **pos = (params->pos);
		int N = (params->N);
		int psym = (params->psym);

		//if(N==2)
		//	printf("%d	%d	%d\n", pos[1][0], pos[1][1], pos[1][2]);


		if(psym==0)
		{
			for(i=0; i<N; i++)
			{

				for(j=0; j<N; j++)
				{
					//printf("N%d	%d	%d	%d\n", N, pos[j][0] - pos[i][0], pos[j][1] - pos[i][1], getdiagt(pos[i][2], pos[j][2]));

					G[i][j] = graphenegf(En, pos[j][0] - pos[i][0], pos[j][1] - pos[i][1], 0, getdiagt(pos[i][2], pos[j][2]) ) + I*graphenegf(En, pos[j][0] - pos[i][0], pos[j][1] - pos[i][1], 1, getdiagt(pos[i][2], pos[j][2]));

					//printf("\n");
				}
			}
		}
		if(psym==1)
		{
			for(j=0; j<N; j++)
			{
				G[0][j] = graphenegf(En, pos[j][0] - pos[0][0], pos[j][1] - pos[0][1], 0, getdiagt(pos[0][2], pos[j][2])) + I*graphenegf(En, pos[j][0] - pos[0][0], pos[j][1] - pos[0][1], 1, getdiagt(pos[0][2], pos[j][2]));
				for(i=1; i<N; i++)
				{
					G[i][(j+1)%N] = G[0][j];
				}
			}
		}

	}

//nice wrapper routine for a matrix of Green functions for sites whose positions in form (a1, a2, sublattice) are in pos[][]
	void grapheneGFmatrix_strain( double _Complex En, double _Complex **G, void *p)
	{
		gfmats_p *params = (gfmats_p *)p;
		int i, j;

		int **pos = (params->pos);
		int N = (params->N);
		int psym = (params->psym);
		double t1 = (params->t1);
		double t2 = (params->t2);

		//if(N==2)
		//	printf("%d	%d	%d\n", pos[1][0], pos[1][1], pos[1][2]);


		if(psym==0)
		{
			for(i=0; i<N; i++)
			{

				for(j=0; j<N; j++)
				{
					//printf("N%d	%d	%d	%d\n", N, pos[j][0] - pos[i][0], pos[j][1] - pos[i][1], getdiagt(pos[i][2], pos[j][2]));

					G[i][j] = graphenegf_strain(En, pos[j][0] - pos[i][0], pos[j][1] - pos[i][1], 0, getdiagt(pos[i][2], pos[j][2]), t1, t2 ) + I*graphenegf_strain(En, pos[j][0] - pos[i][0], pos[j][1] - pos[i][1], 1, getdiagt(pos[i][2], pos[j][2]), t1, t2);
					//printf("\n");
				}
			}
		}
		if(psym==1)
		{
			for(j=0; j<N; j++)
			{
				G[0][j] = graphenegf_strain(En, pos[j][0] - pos[0][0], pos[j][1] - pos[0][1], 0, getdiagt(pos[0][2], pos[j][2]), t1, t2) + I*graphenegf_strain(En, pos[j][0] - pos[0][0], pos[j][1] - pos[0][1], 1, getdiagt(pos[0][2], pos[j][2]), t1, t2);
				for(i=1; i<N; i++)
				{
					G[i][(j+1)%N] = G[0][j];
				}
			}
		}

	}

//calculate moments from up and down spin potentials
double getlocalm (double *m, double *shift, double *zeemanterm, double hubU, void (*gfroutine2) (double _Complex, double _Complex **, void *), void *gfparams, int i, int numsites, double fermi, double imagEn)
{
	intdosother_params para = {};

	double nup, ndown, error1, error2;
	int j;

	//DOS integrator - initialisation
	gsl_function integ;
	integ.function = &intdos_other;
	integ.params = &para;

	//perturbation matrices for up and down spin electrons
	double _Complex **Vup = createSquareMatrix(numsites);
	double _Complex **Vdown = createSquareMatrix(numsites);

	//parameters passed to DOS integrator
	para.i = i;
	para.N = numsites;
	para.fermi = fermi;
	para.imagEn = imagEn;
	para.gfroutine2 = gfroutine2;
	para.gfparams = gfparams;
	

	//set-up perturbation potential
	for(j=0; j<numsites; j++)
	{
		Vup[j][j] = 	-(hubU/2)*(m[j]) + shift[j] -zeemanterm[j]/2;
		Vdown[j][j] = (hubU/2)*(m[j]) + shift[j] -zeemanterm[j]/2;
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



	return nup - ndown;
}


	//function for returning integrand for occupation (etc) at sites other than the magnetic impurity
	//useful for induced moments at neighbouring sites
	double intdos_other(double impart, void *p)
	{
		intdosother_params *params = (intdosother_params *)p;

		double ans;

		int i = (params->i);	//which element to return
		int N = (params -> N);	//matrix size
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


		FreeMatrix(Gpert);
		FreeMatrix(G);
		return ans;
	}


//nice wrapper routine for a 2x2 matrix of Green functions cacluated using SPA for the off diagonals
	void grapheneSPAmatrix( double _Complex En, double _Complex **G, void *p)
	{
		gfmat_p *params = (gfmat_p *)p;
		int i, j;

		int **pos = (params->pos);
		int N = (params->N);
		int psym = (params->psym);

		if(N > 2)
			exit(1);


		G[0][0] = graphenegf(En, 0, 0, 0, 0) +I*graphenegf(En, 0, 0, 1, 0);

		if(N==2)
		{
			G[1][1] = G[0][0];
			G[0][1] = graphene_SPA_AC(En, pos[1][0] - pos[0][0], pos[1][1] - pos[0][1], 0, 0 ) + I*graphene_SPA_AC(En, pos[1][0] - pos[0][0], pos[1][1] - pos[0][1], 1, 0 );
			G[1][0] = G[0][1];
		}
	

	}


	//green functions for graphene in armchair direction using SPA approach - USE WITH CAUTION!
	double graphene_SPA_AC (double _Complex En, int a1, int a2, int reim, int diagt)
	{
		if(a1 != a2)
			exit(1);
		if(diagt != 0)
			exit(1);
	
		double _Complex cont1, cont2;
		double _Complex kz0, A1, A2, B1, B2;

		//stationary point #1 - E>t
			kz0=0.0;
			A1 = cacos(((En*En/(t*t) - 1.0)/4) - 1.0);
			B1 = -(En*En + 3*t*t)/(2*(csqrt((t*t-En*En)*(En*En - 9*t*t))));
			if(cimag(A1) < 0.0)
				A1=-A1;
			if((creal(B1))*(creal(En)) <= 0.0)
				B1=-B1;
	
			cont1 = ((I*En)/(4*M_PI*t*t)) * ((cexp(I*A1*(a1+a2)))/(ccos(kz0) * csin(A1))) * (csqrt ((I*M_PI)/(B1*(a1+a2))));


		//stationary point #2 - E<t
			kz0=cacos((csqrt(t*t - En*En))/(2*t));
			if(fabs(creal(kz0)) > M_PI/2)
				kz0 = cacos((-csqrt(t*t - En*En))/(2*t));

			A2 = cacos((csqrt(t*t - En*En))/(t));
			B2 = -(En*En + 3*t*t)/(2*En*(csqrt(t*t-En*En)));
			if(cimag(A2) < 0.0)
				A2=-A2;
			if( (creal(B2)/creal(A2)) < 0.0)
				B2=-B2;
	
			cont2 = ((2*I*En)/(4*M_PI*t*t)) * ((cexp(I*A2*(a1+a2)))/(ccos(kz0) * csin(A2))) * (csqrt ((I*M_PI)/(B2*(a1+a2))));

		if(reim==0)
			return creal(cont1+cont2);
		if(reim==1)
			return cimag(cont1+cont2);

	}




	double graphenegf(double _Complex En, int a1, int a2, int reim, int diagt)
	{
		double error, gfn;
		size_t test;
		intg_params para = {};
		gsl_integration_workspace *wo;
		gsl_function integ;
		integ.params = &para;
		wo = gsl_integration_workspace_alloc (1000000);

		para.x = a2*dis*sqrt(3)/2;
		para.y = a1*dis + a2*dis/2;

		//printf("%d	%d	%d	%d\n", a1, a2, reim, diagt);

		para.En = En;
		integ.function = &intg;
		para.reim = reim;
		para.diagt = diagt;
		//printf("start GF int %d %d %d E %lf %lf \n ", a1, a2, diagt, creal(En), cimag(En));
		//printf("start integral\n");
		//gsl_integration_qags (&integ, -M_PI/(d*sqrt(3)), M_PI/(d*sqrt(3)) , 1.0e-7, 0.0000001, 1000000,  wo, &gfn, &error);
	//	printf("start integration	%e	%e	%d	%d	%d	%d\n", creal(En), cimag(En), a1, a2, reim, diagt);
		gsl_integration_qag (&integ, -M_PI/(dis*sqrt(3)), M_PI/(dis*sqrt(3)) , 1.0e-10, 1.0e-7, 1000000, 5, wo, &gfn, &error);

		//printf("\n");
		//printf("end integral\n");

		//printf("end GF int\n");


		gsl_integration_workspace_free (wo);

		return gfn;


	}
	
	



//returns the integrand for the calculation of the graphene GF.
	double intg(double kx, void *p)
	{
		intg_params *params = (intg_params *)p;
		double x = (params->x);
		double y = (params->y);
		double _Complex En = (params->En);
		int reim = (params->reim);
		int diagt = (params->diagt);
		double _Complex ky, ky1, ky2, kya, kyb;
		double _Complex ans, cosky1, cosky2, sinky1, sinky2, cosky, S1, S2, Sa, Sb, Sda, Sdb, Sd1, Sd2, S, ans1, ans2, ans3;
		double coskx; 

		S1=0.0;
		S2=0.0;
		
	

		coskx = cos(sqrt(3) * kx * dis / 2);
		cosky1 = (-0.5*(coskx + csqrt(En*En/(t*t) - 1 + coskx*coskx ) ));
		cosky2 = (-0.5*(coskx - csqrt(En*En/(t*t) - 1 + coskx*coskx ) ));
		

		if(diagt == 0)
		{
			getS1S2 (&S1, &S2, x, y, kx, cosky1, cosky2, En);

			ans = (S1 + S2)*En;
		}

		else if (diagt == 2)
		{
			getS1S2 (&S1, &S2, x, y, kx, cosky1, cosky2, En);
			ans1 = (S1 + S2)*t;

			getS1S2 (&S1, &S2, x, y - dis, kx, cosky1, cosky2, En);
			ans2 = (S1 + S2)*t;

			getS1S2 (&S1, &S2, x - sqrt(3)*dis/2, y - dis/2, kx, cosky1, cosky2, En);
			ans3 = (S1 + S2)*t;

			ans = ans1+ans2+ans3;
		}
		else if (diagt == 1)
		{
			
			getS1S2 (&S1, &S2, x, y, kx, cosky1, cosky2, En);
			ans1 = (S1 + S2)*t;

			getS1S2 (&S1, &S2, x, y + dis, kx, cosky1, cosky2, En);
			ans2 = (S1 + S2)*t;

			getS1S2 (&S1, &S2, x + sqrt(3)*dis/2, y + dis/2, kx, cosky1, cosky2, En);
			ans3 = (S1 + S2)*t;

			
			ans = ans1+ans2+ans3;
		}


		if (reim==0)
		{
			//printf("%lf	%lf\n", kx, creal(ans));
			return creal(ans);
		}

		if (reim==1)
		{
			//printf("%lf	%lf\n", kx, cimag(ans));
			return cimag(ans);
		}

	}
	
	
	
void getS1S2 (double _Complex *S1, double _Complex *S2, double x, double y, double kx, double _Complex cosky1, double _Complex cosky2, double _Complex En)
{
		double _Complex ky1, ky2, sinky1, sinky2, tempS1, tempS2;
		double coskx = cos(sqrt(3) * kx * dis / 2);

		if(y<0.0)
		{
			y=-y;
		}	

		ky1 = (2.0/dis)*cacos(cosky1);
		ky2 = (2.0/dis)*cacos(cosky2);
		sinky1 = csin(ky1*dis/2.0);
		sinky2 = csin(ky2*dis/2.0);


		if(cimag(ky1) < 0)
			ky1=-ky1;

		sinky1 = csin(ky1*dis/2.0);
		*S1 = (sqrt(3) * dis  * I / (8*M_PI*t*t )  * (cexp (I*(kx * x + ky1*y))) )/  (coskx * sinky1 + 2*cosky1*sinky1);

		
		if(cimag(ky2) < 0)
			ky2=-ky2;
		sinky2 = csin(ky2*dis/2.0);
		*S2 = (sqrt(3) * dis  * I / (8*M_PI*t*t )  * (cexp (I*(kx * x + ky2*y))) )/  (coskx * sinky2 + 2*cosky2*sinky2);

	
}


void graph_NNS (double _Complex **g, double _Complex En, double a1, double a2, int diagt)
{


		double tempr, tempi;
		int i, j;
		int arel[3][3];		//VECTORS CONNECTING 'a' TO IT'S NEAREST NEIGHBOURS.
		int brel[3][3];		//first index = atom num, 2nd index denotes coord a1, a2, diagt. 
					//a1, a2 are relative, diagt is absolute.
		int dt01, dt10, dt02, dt20, dt05, dt50, dt15, dt51, dt12, dt21, dt25, dt52;


		if (diagt ==0)
		{
			arel[0][0] = 0; arel[0][1] = 0; arel[0][2] = 1;
			arel[1][0] = 1; arel[1][1] = 0; arel[1][2] = 1;
			arel[2][0] = 0; arel[2][1] = 1; arel[2][2] = 1;

			brel[0][0] = 0; brel[0][1] = 0; brel[0][2] = 1;
			brel[1][0] = 1; brel[1][1] = 0; brel[1][2] = 1;
			brel[2][0] = 0; brel[2][1] = 1; brel[2][2] = 1;

			dt01 =0; dt10=0; dt02=1; dt20=2; dt05=1; dt50=2; dt15=1; dt51=2; dt12=1; dt21=2; dt25=0; dt52=0;

		}


		if (diagt ==1)
		{
			arel[0][0] = 0; arel[0][1] = 0; arel[0][2] = 1;
			arel[1][0] = 1; arel[1][1] = 0; arel[1][2] = 1;
			arel[2][0] = 0; arel[2][1] = 1; arel[2][2] = 1;

			brel[0][0] = 0; brel[0][1] = 0; brel[0][2] = 2;
			brel[1][0] =-1; brel[1][1] = 0; brel[1][2] = 2;
			brel[2][0] = 0; brel[2][1] =-1; brel[2][2] = 2;

			dt01 =1; dt10=2; dt02=1; dt20=2; dt05=0; dt50=0; dt15=2; dt51=1; dt12=0; dt21=0; dt25=2; dt52=1;

		}

		if (diagt ==2)
		{
			arel[0][0] = 0; arel[0][1] = 0; arel[0][2] = 2;
			arel[1][0] =-1; arel[1][1] = 0; arel[1][2] = 2;
			arel[2][0] = 0; arel[2][1] =-1; arel[2][2] = 2;

			brel[0][0] = 0; brel[0][1] = 0; brel[0][2] = 1;
			brel[1][0] = 1; brel[1][1] = 0; brel[1][2] = 1;
			brel[2][0] = 0; brel[2][1] = 1; brel[2][2] = 1;

			dt01 =2; dt10=1; dt02=2; dt20=1; dt05=0; dt50=0; dt15=1; dt51=2; dt12=0; dt21=0; dt25=1; dt52=2;

		}


		//build unperturbed GF matrix
			//in the following matrices '0' and '1' are the impurity sites.
			//'2' is the unit cell companion of '0' and '3' and '4' are its other nearest neighbours
			//'5' is the unit cell companion of '1' and '6' and '7' are its other NNs

			tempr = graphenegf(En, 0, 0, 0,  0);
			tempi = graphenegf(En, 0, 0, 1,  0);	

		for(i=0; i<8; i++)
		{
			g[i][i] = tempr + tempi*I;
		}

	

			//imp atoms
				tempr = graphenegf(En, a1, a2, 0, dt01);
				tempi = graphenegf(En, a1, a2, 1, dt01);
				g[0][1] = tempr + tempi*I;
			//related (25, 36, 47)
				if(diagt == 0)
				{
					g[2][5] = g[0][1];
					g[3][6] = g[0][1];
					g[4][7] = g[0][1];
				}
	
				else if(diagt == 1 || diagt ==2)
				{
					tempr = graphenegf(En, -arel[0][0] + a1 + brel[0][0], -arel[0][1] + a1 + brel[0][1] , 0, dt25);
					tempi = graphenegf(En, -arel[0][0] + a1 + brel[0][0], -arel[0][1] + a1 + brel[0][1] , 1, dt25);
					g[2][5] = tempr + tempi*I;
	
					tempr = graphenegf(En, -arel[1][0] + a1 + brel[1][0], -arel[1][1] + a1 + brel[1][1] , 0, dt25);
					tempi = graphenegf(En, -arel[1][0] + a1 + brel[1][0], -arel[1][1] + a1 + brel[1][1] , 1, dt25);
					g[3][6] = tempr + tempi*I;

					tempr = graphenegf(En, -arel[2][0] + a1 + brel[2][0], -arel[2][1] + a1 + brel[2][1] , 0, dt25);
					tempi = graphenegf(En, -arel[2][0] + a1 + brel[2][0], -arel[2][1] + a1 + brel[2][1] , 1, dt25);
					g[4][7] = tempr + tempi*I;
				}



			//inverse imp atoms
				tempr = graphenegf(En, -a1, -a2, 0, dt10);
				tempi = graphenegf(En, -a1, -a2, 1, dt10);
				g[1][0] = tempr + tempi*I;
			//related (52, 63, 74)
				if(diagt ==0)
				{
					g[5][2] = g[1][0];
					g[6][3] = g[1][0];
					g[7][4] = g[1][0];
				}
				else if(diagt == 1 || diagt ==2)
				{
					tempr = graphenegf(En, arel[0][0] - a1 - brel[0][0], arel[0][1] - a1 - brel[0][1] , 0, dt52);
					tempi = graphenegf(En, arel[0][0] - a1 - brel[0][0], arel[0][1] - a1 - brel[0][1] , 1, dt52);
					g[5][2] = tempr + tempi*I;
	
					tempr = graphenegf(En, arel[1][0] - a1 - brel[1][0], arel[1][1] - a1 - brel[1][1] , 0, dt52);
					tempi = graphenegf(En, arel[1][0] - a1 - brel[1][0], arel[1][1] - a1 - brel[1][1] , 1, dt52);
					g[6][3] = tempr + tempi*I;

					tempr = graphenegf(En, arel[2][0] - a1 - brel[2][0], arel[2][1] - a1 - brel[2][1] , 0, dt52);
					tempi = graphenegf(En, arel[2][0] - a1 - brel[2][0], arel[2][1] - a1 - brel[2][1] , 1, dt52);
					g[4][7] = tempr + tempi*I;

				}

			//unit cell companions
				tempr = graphenegf(En, arel[0][0], arel[0][1], 0, dt02);
				tempi = graphenegf(En, arel[0][0], arel[0][1], 1, dt02);
				g[0][2] = tempr + tempi*I;
				if(diagt==0)
				{		
					g[1][5] = g[0][2] ;
				}
				else if (diagt == 1 || diagt ==2)
				{
					g[5][1] = g[0][2];
				}

				tempr = graphenegf(En, -arel[0][0], -arel[0][1], 0, dt20);
				tempi = graphenegf(En, -arel[0][0], -arel[0][1], 1, dt20);
				g[2][0] = tempr + tempi*I;
				if(diagt==0)
				{		
					g[5][1] = g[2][0] ;
				}
				else if (diagt == 1 || diagt ==2)
				{
					g[1][5] = g[2][0] ;
				}


			//imp sites and their nearest (diff. cell) neighbours
				tempr = graphenegf(En, arel[1][0], arel[1][1], 0, dt02);
				tempi = graphenegf(En, arel[1][0], arel[1][1], 1, dt02);
				g[0][3] = tempr + tempi*I;
				if(diagt == 0)
				{
					g[1][6] = g[0][3] ;
				}
				else if (diagt == 1 || diagt ==2)
				{
					g[6][1] = g[0][3] ;
				}

				tempr = graphenegf(En, -arel[1][0], -arel[1][1], 0, dt20);
				tempi = graphenegf(En, -arel[1][0], -arel[1][1], 1, dt20);
				g[3][0] = tempr + tempi*I;
				if(diagt == 0)
				{
					g[6][1] = g[3][0] ;
				}
				else if (diagt == 1 || diagt ==2)
				{
					g[1][6] = g[3][0] ;
				}


				tempr = graphenegf(En, arel[2][0], arel[2][1], 0, dt02);
				tempi = graphenegf(En, arel[2][0], arel[2][1], 1, dt02);
				g[0][4] = tempr + tempi*I;
				if(diagt == 0)
				{
					g[1][7] = g[0][4] ;
				}
				else if (diagt == 1 || diagt ==2)
				{
					g[7][1] = g[0][4] ;
				}

				tempr = graphenegf(En, -arel[2][0], -arel[2][1], 0, dt20);
				tempi = graphenegf(En, -arel[2][0], -arel[2][1], 1, dt20);
				g[4][0] = tempr + tempi*I;
				if(diagt == 0)
				{
					g[7][1] = g[4][0] ;
				}
				else if (diagt == 1 || diagt ==2)
				{
					g[1][7] = g[4][0] ;
				}


			//NNs of imps interactions with each other
				tempr = graphenegf(En, -arel[0][0] + arel[1][0], -arel[0][1] + arel[1][1], 0, 0);
				tempi = graphenegf(En, -arel[0][0] + arel[1][0], -arel[0][1] + arel[1][1], 1, 0);
				g[2][3] = tempr + tempi*I;
				if(diagt == 0)
				{
					g[5][6] = g[2][3] ;
				}
				else if (diagt == 1 || diagt ==2 )
				{
					g[6][5] = g[2][3] ;
				} 

				tempr = graphenegf(En, arel[0][0] - arel[1][0], arel[0][1] - arel[1][1], 0, 0);
				tempi = graphenegf(En, arel[0][0] - arel[1][0], arel[0][1] - arel[1][1], 1, 0);
				g[3][2] = tempr + tempi*I;
				if(diagt == 0)
				{
					g[6][5] = g[3][2] ;
				}
				else if (diagt == 1 || diagt ==2 )
				{
					g[5][6] = g[3][2] ;
				} 

				tempr = graphenegf(En, -arel[0][0] + arel[2][0], -arel[0][1] + arel[2][1], 0, 0);
				tempi = graphenegf(En, -arel[0][0] + arel[2][0], -arel[0][1] + arel[2][1], 1, 0);
				g[2][4] = tempr + tempi*I;
				if(diagt == 0)
				{
					g[5][7] = g[2][4] ;
				}
				else if (diagt == 1 || diagt ==2 )
				{
					g[7][5] = g[2][4] ;
				} 


				tempr = graphenegf(En, arel[0][0] - arel[2][0], arel[0][1] - arel[2][1], 0, 0);
				tempi = graphenegf(En, arel[0][0] - arel[2][0], arel[0][1] - arel[2][1], 1, 0);
				g[4][2] = tempr + tempi*I;
				if(diagt == 0)
				{
					g[7][5] = g[4][2] ;
				}
				else if (diagt == 1 || diagt ==2 )
				{
					g[5][7] = g[4][2] ;
				} 

				tempr = graphenegf(En, -arel[1][0] + arel[2][0], -arel[1][1] + arel[2][1], 0, 0);
				tempi = graphenegf(En, -arel[1][0] + arel[2][0], -arel[1][1] + arel[2][1], 1, 0);
				g[3][4] = tempr + tempi*I;
				if(diagt == 0)
				{
					g[6][7] = g[3][4] ;
				}
				else if (diagt == 1 || diagt ==2 )
				{
					g[7][6] = g[3][4] ;
				} 

				tempr = graphenegf(En, arel[1][0] - arel[2][0], arel[1][1] - arel[2][1], 0, 0);
				tempi = graphenegf(En, arel[1][0] - arel[2][0], arel[1][1] - arel[2][1], 1, 0);
				g[4][3] = tempr + tempi*I;
				if(diagt == 0)
				{
					g[7][6] = g[4][3] ;
				}
				else if (diagt == 1 || diagt ==2 )
				{
					g[6][7] = g[4][3] ;
				}


			//connecting imp to nns of other imp
				for(i=0; i<3; i++)
				{
					tempr = graphenegf(En, a1 + brel[i][0], a2 + brel[i][1], 0, dt05);
					tempi = graphenegf(En, a1 + brel[i][0], a2 + brel[i][1], 1, dt05);
					g[0][5+i] = tempr + tempi*I;

					tempr = graphenegf(En, -a1 - brel[i][0], -a2 - brel[i][1], 0, dt50);
					tempi = graphenegf(En, -a1 - brel[i][0], -a2 - brel[i][1], 1, dt50);
					g[5+i][0] = tempr + tempi*I;

					tempr = graphenegf(En, -arel[i][0] + a1, -arel[i][1] + a2, 0, dt12);
					tempi = graphenegf(En, -arel[i][0] + a1, -arel[i][1] + a2, 1, dt12);
					g[1][2+i] = tempr + tempi*I;

					tempr = graphenegf(En, arel[i][0] - a1, arel[i][1] - a2, 0, dt21);
					tempi = graphenegf(En, arel[i][0] - a1, arel[i][1] - a2, 1, dt21);
					g[2+i][1] = tempr + tempi*I;


					//connecting nns of one imp to nns of the other (1 connection done already from earlier)
					for(j=0; j<3; j++)
					{
						if(i != j)
						{
							tempr = graphenegf(En, -arel[i][0] + a1 + brel[j][0], -arel[i][1] + a2 + brel[j][1], 0, dt25);
							tempi = graphenegf(En, -arel[i][0] + a1 + brel[j][0], -arel[i][1] + a2 + brel[j][1], 1, dt25);
							g[2+i][5+j] = tempr + tempi*I;

							tempr = graphenegf(En, -brel[i][0] - a1 + arel[j][0], -brel[i][1] -a2 + arel[j][1] , 0, dt52);
							tempi = graphenegf(En, -brel[i][0] - a1 + arel[j][0], -brel[i][1] -a2 + arel[j][1] , 1, dt52);
							g[5+i][2+j] = tempr + tempi*I;

						}
					}

				}
}


int getdiagt (int d1, int d2)
{
	//printf("%d	%d\n", d1, d2);
	if(d1 == d2)
		return 0;
	if(d1==0 && d2==1)
		return 1;
	if(d1==1 && d2==0)
		return 2;
}

void graphenebridgematrix( double _Complex En, double _Complex **G, void *p)
{
	gfmatb_p *params = (gfmatb_p *)p;
	int i, j;

	int **pos = (params->pos);
	int N = (params->N);
	double tprime = (params->tprime);	//no psym , 3rd part of pos (diagt) ignored. first 2 parts of pos used for separation vector


	double _Complex G_diag, G01, G10;
	double _Complex **Gtemp = createSquareMatrix(3*N);
	double _Complex **V = createSquareMatrix(3*N);
	double _Complex **Gtemp2 = createSquareMatrix(3*N);


	G_diag = graphenegf(En, 0, 0, 0, 0 ) + I*graphenegf(En, 0, 0, 1, 0);
	G01 = graphenegf(En, 0, 0, 0, 1 ) + I*graphenegf(En, 0, 0, 1, 1);
	G10 = graphenegf(En, 0, 0, 0, 2 ) + I*graphenegf(En, 0, 0, 1, 2);

	for(i=0; i<N; i++)
	{
		Gtemp[i*N][i*N] = G_diag;
		Gtemp[i*N+1][i*N+1] = G_diag ;
		Gtemp[i*N][i*N+1] = G01;
		Gtemp[i*N+1][i*N] = G10;

		for(j=0; j<N; j++)
		{

			if(i!=j)
			{

				Gtemp[i*N][j*N] = graphenegf(En, pos[j][0] - pos[i][0], pos[j][0] - pos[i][0], 0, 0 ) + I*graphenegf(En, pos[j][0] - pos[i][0], pos[j][0] - pos[i][0], 1, 0);
				Gtemp[i*N][j*N+1] = graphenegf(En, pos[j][0] - pos[i][0], pos[j][0] - pos[i][0], 0, 1 ) + I*graphenegf(En, pos[j][0] - pos[i][0], pos[j][0] - pos[i][0], 1, 1);
				Gtemp[i*N+1][j*N] = graphenegf(En, pos[j][0] - pos[i][0], pos[j][0] - pos[i][0], 0, 2 ) + I*graphenegf(En, pos[j][0] - pos[i][0], pos[j][0] - pos[i][0], 1, 2);
				Gtemp[i*N+1][j*N+1] = Gtemp[i*N][j*N] ;

			}

		}

	}

	//Gtemp[2*N][2*N] = 1 / En ;
	//Gtemp[2*N+1][2*N+1] = 	1 / En ;
	
	for(i=0; i<N; i++)
	{
		Gtemp[2*N+i][2*N+i] = 1 / En ;
		V[2*N +i][2*i] = tprime;
		V[2*N +i][2*i + 1] = tprime;
		V[2*i][2*N +i] = tprime;
		V[2*i + 1][2*N +i] = tprime;
	}
		//printEMatrix(Gtemp, 3*N);
		//printEMatrix(V, 3*N);

	dyson(Gtemp, V, Gtemp2, 3*N);

	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			G[i][j] = Gtemp2[2*N+i][2*N+j] ;
		}
	}
	
	FreeMatrix(Gtemp);
	FreeMatrix(V);
	FreeMatrix(Gtemp2);
}


void graphenetopmatrix( double _Complex En, double _Complex **G, void *p)
{
	gfmatb_p *params = (gfmatb_p *)p;
	int i, j;

	int **pos = (params->pos);
	int N = (params->N);
	double tprime = (params->tprime);	

	double _Complex **Gtemp = createSquareMatrix(2*N);
	double _Complex **V = createSquareMatrix(2*N);
	double _Complex **Gtemp2 = createSquareMatrix(2*N);


	for(i=0; i<N; i++)
	{

		for(j=0; j<N; j++)
		{
			//printf("N%d	%d	%d	%d\n", N, pos[j][0] , pos[j][1] , getdiagt(pos[i][2], pos[j][2]));

			Gtemp[i][j] = graphenegf(En, pos[j][0] - pos[i][0], pos[j][1] - pos[i][1], 0, getdiagt(pos[i][2], pos[j][2]) ) + I*graphenegf(En, pos[j][0] - pos[i][0], pos[j][1] - pos[i][1], 1, getdiagt(pos[i][2], pos[j][2]));

			//printf("\n");
		}
	}
	
	

	for(i=0; i<N; i++)
	{
		Gtemp[N+i][N+i] = 1 / En ;
		V[N +i][i] = tprime;
		V[i][N+i] = tprime;
	}
		//printEMatrix(Gtemp, 2*N);
		//printEMatrix(V, 2*N);

	dyson(Gtemp, V, Gtemp2, 2*N);

	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			G[i][j] = Gtemp2[N+i][N+j] ;
		}
	}
	
	FreeMatrix(Gtemp);
	FreeMatrix(V);
	FreeMatrix(Gtemp2);
}


void graphenecentrematrix( double _Complex En, double _Complex **G, void *p)
{
	gfmatb_p *params = (gfmatb_p *)p;
	int i, j, k, l;

	int **pos = (params->pos);
	int N = (params->N);
	double tprime = (params->tprime);	//no psym , 3rd part of pos (diagt) ignored. first 2 parts of pos used for separation vector


	double _Complex G_diag, G01, G10;
	double _Complex **Gtemp = createSquareMatrix(7*N);
	double _Complex **V = createSquareMatrix(7*N);
	double _Complex **Gtemp2 = createSquareMatrix(7*N);
	int **temppos = createNonSquareIntMatrix(6, 3);

	//relative positions of the 6 atoms in the hexagon, with 0,0,0 the bottom left and lattice vectors | and /
	//temppos[0][0] = 0; temppos[0][1] = 0; temppos[0][2] = 0;
	//temppos[1][0] = 0; temppos[1][1] = 0; temppos[1][2] = 1;
	//temppos[2][0] = 0; temppos[2][1] = 1; temppos[2][2] = 0;
	//temppos[3][0] = 1; temppos[3][1] = 0; temppos[3][2] = 1;
	//temppos[4][0] = 1; temppos[4][1] = 0; temppos[4][2] = 0;
	//temppos[5][0] = 1; temppos[5][1] = -1; temppos[5][2] = 1;


	temppos[0][0] = 0; temppos[0][1] = 0; temppos[0][2] = 0;
	temppos[1][0] = -1; temppos[1][1] = 0; temppos[1][2] = 1;
	temppos[2][0] = -1; temppos[2][1] = 1; temppos[2][2] = 0;
	temppos[3][0] = -1; temppos[3][1] = 1; temppos[3][2] = 1;
	temppos[4][0] = 0; temppos[4][1] = 1; temppos[4][2] = 0;
	temppos[5][0] = 0; temppos[5][1] = 0; temppos[5][2] = 1;

	

	G_diag = graphenegf(En, 0, 0, 0, 0 ) + I*graphenegf(En, 0, 0, 1, 0);

	for(i=0; i<6; i++)
	{
		Gtemp[i][i] = G_diag;
		
		for(j=0; j< 6; j++)
		{
			if(i !=j)
				Gtemp[i][j] = graphenegf(En, temppos[j][0] - temppos[i][0], temppos[j][1] - temppos[i][1], 0, getdiagt(temppos[i][2], temppos[j][2]) ) + I*graphenegf(En, temppos[j][0] - temppos[i][0], temppos[j][1] - temppos[i][1], 1, getdiagt(temppos[i][2], temppos[j][2]));
		}
	}

	for(i=1; i<N; i++)
	{
		MatrixCopyPart(Gtemp, Gtemp, 0, 0, 6*i, 6*i, 6, 6);
	}

	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			if (i != j)
			{
				for(k=0; k<6; k++)
				{
					for(l=0; l<6; l++)
					{
						Gtemp[6*i +k][6*j + l] = graphenegf (En, pos[j][0] + temppos[l][0] - pos[i][0] - temppos[k][0], pos[j][1] + temppos[l][1] - pos[i][1] - temppos[k][1], 0, getdiagt(temppos[k][2], temppos[l][2]) ) + I* graphenegf (En, pos[j][0] + temppos[l][0] - pos[i][0] - temppos[k][0], pos[j][1] + temppos[l][1] - pos[i][1] - temppos[k][1], 1, getdiagt(temppos[k][2], temppos[l][2]) );
					}
				}
			}
		}
	}


	//Gtemp[2*N][2*N] = 1 / En ;
	//Gtemp[2*N+1][2*N+1] = 	1 / En ;
	
	for(i=0; i<N; i++)
	{
		Gtemp[6*N+i][6*N+i] = 1 / En ;
	
		for(j=0; j<6; j++)
		{
			V[6*N + i][6*i + j] = tprime;
			V[6*i + j][6*N +i] = tprime;
		}
	}
		//printEMatrix(Gtemp, 6);
		//printEMatrix(V, 7*N);

	dyson(Gtemp, V, Gtemp2, 7*N);
	//printEMatrix(Gtemp2, 7*N);


	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			G[i][j] = Gtemp2[6*N+i][6*N+j] ;
		}
	}
	
	free(temppos[0]);
	free(temppos);

	FreeMatrix(Gtemp);
	FreeMatrix(V);
	FreeMatrix(Gtemp2);
}

double _Complex graphenecentresum ( double _Complex En, double _Complex **G, void *p)
{
	gfmatb_p *params = (gfmatb_p *)p;
	int i, j, k, l;

	int **pos = (params->pos);
	int N = (params->N); //has to be 2
	double tprime = (params->tprime);	//no psym , 3rd part of pos (diagt) ignored. first 2 parts of pos used for separation vector


	double _Complex G_diag, G01, G10;
	int **temppos = createNonSquareIntMatrix(6, 3);
	double _Complex ans = 0.0;

	//relative positions of the 6 atoms in the hexagon, with 0,0,0 the bottom left and lattice vectors | and /
	//temppos[0][0] = 0; temppos[0][1] = 0; temppos[0][2] = 0;
	//temppos[1][0] = 0; temppos[1][1] = 0; temppos[1][2] = 1;
	//temppos[2][0] = 0; temppos[2][1] = 1; temppos[2][2] = 0;
	//temppos[3][0] = 1; temppos[3][1] = 0; temppos[3][2] = 1;
	//temppos[4][0] = 1; temppos[4][1] = 0; temppos[4][2] = 0;
	//temppos[5][0] = 1; temppos[5][1] = -1; temppos[5][2] = 1;


	temppos[0][0] = 0; temppos[0][1] = 0; temppos[0][2] = 0;
	temppos[1][0] = -1; temppos[1][1] = 0; temppos[1][2] = 1;
	temppos[2][0] = -1; temppos[2][1] = 1; temppos[2][2] = 0;
	temppos[3][0] = -1; temppos[3][1] = 1; temppos[3][2] = 1;
	temppos[4][0] = 0; temppos[4][1] = 1; temppos[4][2] = 0;
	temppos[5][0] = 0; temppos[5][1] = 0; temppos[5][2] = 1;

	
/*
	for(i=0; i<N; i++)
	{
		for(j=0; j<N; j++)
		{
			if (i != j)
			{
				for(k=0; k<6; k++)
				{
					for(l=0; l<6; l++)
					{
						Gtemp[6*i +k][6*j + l] = graphenegf (En, pos[j][0] + temppos[l][0] - pos[i][0] - temppos[k][0], pos[j][1] + temppos[l][1] - pos[i][1] - temppos[k][1], 0, getdiagt(temppos[k][2], temppos[l][2]) ) + I* graphenegf (En, pos[j][0] + temppos[l][0] - pos[i][0] - temppos[k][0], pos[j][1] + temppos[l][1] - pos[i][1] - temppos[k][1], 1, getdiagt(temppos[k][2], temppos[l][2]) );
					}
				}
			}
		}
	}*/

	for(k=0; k<N; k++)
	{
		for(l=0; l<N; l++)
		{
			ans += graphenegf (En, pos[1][0] + temppos[l][0] - pos[0][0] - temppos[k][0], pos[1][1] + temppos[l][1] - pos[0][1] - temppos[k][1], 0, getdiagt(temppos[k][2], temppos[l][2]) ) + I* graphenegf (En, pos[1][0] + temppos[l][0] - pos[0][0] - temppos[k][0], pos[1][1] + temppos[l][1] - pos[0][1] - temppos[k][1], 1, getdiagt(temppos[k][2], temppos[l][2]) );
		}
	}


	//Gtemp[2*N][2*N] = 1 / En ;
	//Gtemp[2*N+1][2*N+1] = 	1 / En ;
	return ans;
	
	
}




double graphenegf_strain(double _Complex En, int a1, int a2, int reim, int diagt, double t1, double t2)
	{
		double error, gfn;
		size_t test;
		intgstrain_params para = {};
		gsl_integration_workspace *wo;
		gsl_function integ;
		integ.params = &para;
		wo = gsl_integration_workspace_alloc (1000000);

		para.x = a1+a2;
		para.y = a2-a1;

		//printf("%d	%d	%d	%d\n", a1, a2, reim, diagt);

		para.En = En;
		integ.function = &intg_strain;
		para.reim = reim;
		para.diagt = diagt;
		para.t1 = t1;
		para.t2 = t2;
		//printf("start GF int %d %d %d E %lf %lf \n ", a1, a2, diagt, creal(En), cimag(En));
		//printf("start integral\n");
		//gsl_integration_qags (&integ, -M_PI/(d*sqrt(3)), M_PI/(d*sqrt(3)) , 1.0e-7, 0.0000001, 1000000,  wo, &gfn, &error);
	//	printf("start integration	%e	%e	%d	%d	%d	%d\n", creal(En), cimag(En), a1, a2, reim, diagt);
		gsl_integration_qag (&integ, -M_PI, M_PI, 1.0e-7, 1.0e-3, 1000000, 5, wo, &gfn, &error);

		//printf("\n");
		//printf("end integral\n");

		//printf("end GF int\n");


		gsl_integration_workspace_free (wo);

		return gfn;


	}
	
	



//returns the integrand for the calculation of the graphene GF.
	double intg_strain(double kx, void *p)
	{
		intgstrain_params *params = (intgstrain_params *)p;
		int x = (params->x);
		int y = (params->y);
		double _Complex En = (params->En);
		int reim = (params->reim);
		int diagt = (params->diagt);
		double t1 = (params->t1);
		double t2 = (params->t2);
		double _Complex ky, ky1, ky2, kya, kyb;
		double _Complex ans, cosky1, cosky2, sinky1, sinky2, cosky, S1, S2, Sa, Sb, Sda, Sdb, Sd1, Sd2, S, ans1, ans2, ans3;
		double coskx; 

		S1=0.0;
		S2=0.0;
		
	

		coskx = cos(kx);
		cosky1 = (-0.5*(t2/t1)*(coskx + csqrt(En*En/(t2*t2) - 1 + coskx*coskx ) ));
		cosky2 = (-0.5*(t2/t1)*(coskx - csqrt(En*En/(t2*t2) - 1 + coskx*coskx ) ));
		

		if(diagt == 0)
		{
			getS1S2_strain (&S1, &S2, x, y, kx, cosky1, cosky2, En, t1, t2);

			ans = (S1 + S2)*En;
		}

//test for diagonal matrix elements first
		else if (diagt == 1)
		{
			getS1S2_strain (&S1, &S2, x, y, kx, cosky1, cosky2, En, t1, t2);
			ans1 = (S1 + S2)*t2;

			getS1S2_strain (&S1, &S2, x + 1, y , kx, cosky1, cosky2, En, t1, t2);
			ans2 = (S1*cosky1 + S2*cosky2)*2*t1;

			

			ans = ans1+ans2;
		}
		else if (diagt == 2)
		{
			
			getS1S2_strain (&S1, &S2, x, y, kx, cosky1, cosky2, En, t1, t2);
			ans1 = (S1 + S2)*t2;

			getS1S2_strain (&S1, &S2, x - 1, y , kx, cosky1, cosky2, En, t1, t2);
			ans2 = (S1*cosky1 + S2*cosky2)*2*t1;

			
			ans = ans1+ans2;
		}


		if (reim==0)
		{
			//printf("%lf	%lf\n", kx, creal(ans));
			return creal(ans);
		}

		if (reim==1)
		{
			//printf("%lf	%lf\n", kx, cimag(ans));
			return cimag(ans);
		}

	}
	
	
	
void getS1S2_strain (double _Complex *S1, double _Complex *S2, double x, double y, double kx, double _Complex cosky1, double _Complex cosky2, double _Complex En, double t1, double t2)
{
		double _Complex ky1, ky2, sinky1, sinky2, tempS1, tempS2;
		double coskx = cos(kx);

 		if(y<0)
 		{
			y=-y;
 		}	

		ky1 = cacos(cosky1);
		ky2 = cacos(cosky2);
		sinky1 = csin(ky1);
		sinky2 = csin(ky2);


		if(cimag(ky1) < 0)
			ky1=-ky1;

		sinky1 = csin(ky1);
		*S1 = (I / (8*M_PI)  * (cexp (I*(kx * x + ky1*y))) )/  (t1 * t2 * coskx * sinky1 + 2 * t1 * t1 * cosky1 * sinky1);

		
		if(cimag(ky2) < 0)
			ky2=-ky2;
		sinky2 = csin(ky2);
		*S2 = (I / (8*M_PI)  * (cexp (I*(kx * x + ky2*y))) )/  (t1 * t2 * coskx * sinky2 + 2 * t1 * t1 * cosky2 * sinky2);

	
}

//direction = 0 for zigzag strain, = 1 for armchair strain
void get_t1t2 (double epsilon, double sigma, int direction, double *t1, double *t2)
{
    if(direction ==0)
    {
      *t1 = t * exp(-3.37 * (0.75*epsilon - 0.25*epsilon*sigma));
      *t2 = t * exp(-3.37 * (-epsilon*sigma));
    }
    
    if(direction ==1)
    {
      *t1 = t * exp(-3.37 * (0.25*epsilon - 0.75*epsilon*sigma));
      *t2 = t * exp(-3.37 * epsilon);

    }
  
}









