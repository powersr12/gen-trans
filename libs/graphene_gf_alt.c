

#include "graphene_gf.h"


	double graphenegf(double _Complex En, double a1, double a2, int reim, int diagt)
	{
		double error, gfn;
		size_t test;
		intg_params para = {};
		gsl_integration_workspace *wo;
		gsl_function integ;
		integ.params = &para;
		wo = gsl_integration_workspace_alloc (10000000);

		para.x = a2*d*sqrt(3)/2;
		para.y = a1*d + a2*d/2;



		para.En = En;
		integ.function = &intg;
		para.reim = reim;
		para.diagt = diagt;
		//printf("start GF int %.2lf %.2lf E %lf %lf \n ");
		gsl_integration_qags (&integ, -M_PI/(d*sqrt(3)), M_PI/(d*sqrt(3)) , 0.0000001, 0.0001, 10000000,  wo, &gfn, &error);
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
		
	

		coskx = cos(sqrt(3) * kx * d / 2);
		cosky1 = (-0.5*(coskx + csqrt(En*En/(t*t) - 1 + coskx*coskx ) ));
		cosky2 = (-0.5*(coskx - csqrt(En*En/(t*t) - 1 + coskx*coskx ) ));
		

		if(diagt == 0)
		{
			getS1S2 (&S1, &S2, x, y, kx, cosky1, cosky2, En);

			ans = (S1 + S2)*En;
		}

		else if (diagt == 1)
		{
			getS1S2 (&S1, &S2, x, y, kx, cosky1, cosky2, En);
			ans1 = (S1 + S2)*t;

			getS1S2 (&S1, &S2, x, y - d, kx, cosky1, cosky2, En);
			ans2 = (S1 + S2)*t;

			getS1S2 (&S1, &S2, x - sqrt(3)*d/2, y - d/2, kx, cosky1, cosky2, En);
			ans3 = (S1 + S2)*t;

			ans = ans1+ans2+ans3;
		}
		else if (diagt == 2)
		{
			
			getS1S2 (&S1, &S2, x, y, kx, cosky1, cosky2, En);
			ans1 = (S1 + S2)*t;

			getS1S2 (&S1, &S2, x, y + d, kx, cosky1, cosky2, En);
			ans2 = (S1 + S2)*t;

			getS1S2 (&S1, &S2, x + sqrt(3)*d/2, y + d/2, kx, cosky1, cosky2, En);
			ans3 = (S1 + S2)*t;

			
			ans = ans1+ans2+ans3;
		}


		if (reim==0)
		{
			return creal(ans);
		}

		if (reim==1)
		{
			return cimag(ans);
		}

	}

void getS1S2 (double _Complex *S1, double _Complex *S2, double x, double y, double kx, double _Complex cosky1, double _Complex cosky2, double _Complex En)
{
		double _Complex ky1, ky2, sinky1, sinky2, tempS1, tempS2;
		double coskx = cos(sqrt(3) * kx * d / 2);

		if(y<0.0)
		{
			y=-y;
		}	

		ky1 = (2.0/d)*cacos(cosky1);
		ky2 = (2.0/d)*cacos(cosky2);
		sinky1 = csin(ky1*d/2.0);
		sinky2 = csin(ky2*d/2.0);


		if(cimag(ky1) < 0)
			ky1=-ky1;

		sinky1 = csin(ky1*d/2.0);
		*S1 = (sqrt(3) * d  * I / (8*M_PI*t*t )  * (cexp (I*(kx * x + ky1*y))) )/  (coskx * sinky1 + 2*cosky1*sinky1);

		
		if(cimag(ky2) < 0)
			ky2=-ky2;
		sinky2 = csin(ky2*d/2.0);
		*S2 = (sqrt(3) * d  * I / (8*M_PI*t*t )  * (cexp (I*(kx * x + ky2*y))) )/  (coskx * sinky2 + 2*cosky2*sinky2);

	
}

void graph_NNS (double _Complex **g, double _Complex En, double a1, double a2, int diagt)
{
		double absloc[8][3];	//absolute position vectors of the 8 lattice points considered
		double rel[8][8][3] ;	//relative separation vectors of the 8 lattice points considered
		double tempr, tempi;

		int i, j;

		for(i=0; i<2; i++)
		{
			absloc[0][i] = 0.0;
			absloc[1][0] = a1;
			absloc[1][1] = a2;
		}
		
		if(diagt == 0)
		{
			absloc[0][2] = 0.0; absloc[1][2] = 0.0;
			for(i=2; i<8; i++)
			{
				absloc[i][2] = 1.0;
			}			

		}		


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
