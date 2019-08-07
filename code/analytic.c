#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <unistd.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_integration.h>
#include "../libs/cubature.h"
#include "../libs/matrices.h"

#include "analytic.h"



//graphenegfac RETURNS AN ARRAY WITH THE REAL AND IMAGINARY PARTS OF GRAPHENE GFS FOR SPECIFIED INDICES


//STANDALONE CODE TO TEST GREEN'S FUNCTIONS.
//COMPILE WITH
// gcc matrices.o cubature.o analytic.c -lm -lgsl -lgslcblas -O3  -o analytic

// main()
// {
//     double _Complex a;
//         double _Complex aa[2][2];
// 
//     double realE=0.1;
//     double eta =1.0E-6;
//     int a1=6, a2 =0, diagt=0;
//     int gfdim = 100;
//     int i;
//     int a1a[3*gfdim], a2a[3*gfdim], diagta[3*gfdim];
//     
//     
//     for(i=0; i< gfdim; i++)
//     {
//         a1a[3*i] = i; a1a[3*i+1] = i; a1a[3*i+2] = i;
//         a2a[3*i] = i; a2a[3*i+1] = i; a2a[3*i+2] = i;
//         diagta[3*i] = 0; diagta[3*i+1] = 1; diagta[3*i+2] = 2;
//     }
//         
//     
//     
//     
// 
//     
//     double _Complex *mat;//=createCompArray(3*gfdim);
//     
//     
//     for(realE = -3.0; realE<=3.0; realE +=0.05)
//     {
// //         a = graphenegfa(realE + eta*I, a1,a2, 0, diagt) + I*graphenegfa(realE + eta*I, a1,a2, 1, diagt);
//         mat = graphenegfac(realE + eta*I, a1a, a2a, diagta, 3*gfdim);
//         
//         printf("%lf\t", realE);
//         for(i=0; i< 3*gfdim; i++)
//         {
//             printf("%lf\t%lf\t", creal(mat[i]), cimag(mat[i]));
//         }
//         printf("\n");
// //         printf("%e  %e  %e  %e  %e \n", realE, creal(a), cimag(a), creal(mat[0]), cimag(mat[0]));
//     }
// 
// }


    
double graphenegfa(double _Complex En, int a1, int a2, int reim, int diagt)
{
        double error, gfn;
        double dis=1.0;
        double t=-1.0;
        
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
        integ.function = &intga;
        para.reim = reim;
        para.diagt = diagt;
        
        gsl_integration_qag (&integ, -M_PI/(dis*sqrt(3)), M_PI/(dis*sqrt(3)) , 1.0e-10, 1.0e-7, 1000000, 5, wo, &gfn, &error);


        gsl_integration_workspace_free (wo);

        return gfn;


}
    

double _Complex *graphenegfac(double _Complex En, int *a1, int *a2, int *diagt, int dim)
{
        int i;
        double *error= createDoubleArray(2*dim) , *gfn= createDoubleArray(2*dim);
        double _Complex *out = createCompArray(dim);
        double dis=1.0;
        double t=-1.0;

        intgac_params para = {};


        para.x = createDoubleArray(dim);
        para.y = createDoubleArray(dim);
        para.diagt = createIntArray(dim);
        para.En = En;
        
        for(i=0; i<dim; i++)
        {
            para.x[i] = a2[i]*dis*sqrt(3)/2;
            para.y[i] = a1[i]*dis + a2[i]*dis/2;
            para.diagt[i] = diagt[i];
//             printf("%lf	%lf	%d\n", para.x[i], para.y[i], para.diagt[i]);
        }

        double *kxmins = createDoubleArray(1);
        double *kxmaxs = createDoubleArray(1);
        kxmins[0]=-M_PI/(dis*sqrt(3));
        kxmaxs[0]=M_PI/(dis*sqrt(3));
        
        
        adapt_integrate (2*dim, &intgac, &para, 1, kxmins, kxmaxs, 0, 1.0e-10, 1.0e-7, gfn, error);

        for(i=0; i<dim; i++)
        {
            out[i] = gfn[2*i] + I*gfn[2*i+1] ;
        }

        free(para.x); free(para.y); free(para.diagt); free(kxmins); free(kxmaxs);
        return out;


}


//checks existing library for GFs, calculates and adds if not already there
double _Complex *graphenegfac_lib(double _Complex En, int *a1, int *a2, int *diagt, int dim, char *PGFlibloc)
{
        int i;
        
        char fullname[200], direcname[300];
        int *indices_to_run = createIntArray(dim);
        int num_to_run = 0;
        double retemp, imtemp;
        double _Complex *out = createCompArray(dim);
        double _Complex *outs;
        char command[400];
        
        int *a1s, *a2s, *ds;
        
        FILE *file;
        
        for (i=0; i < dim; i++)
        {
            sprintf(direcname, "%s%.0e", PGFlibloc, cimag(En));
            sprintf(fullname, "%s/%d_%d_%d_E_%+.8lf.dat", direcname, a1[i], a2[i], diagt[i], creal(En));
            
            if( access( fullname, F_OK ) != -1 ) 
            {
                //printf("##file exists!\n");
                file = fopen(fullname, "r");
                if (fscanf(file, "%lf   %lf", &retemp, &imtemp) < 1)
                {
                    printf("##file not read properly!\n");
                    indices_to_run[num_to_run] = i;
                    num_to_run ++;
                }
                else
                {
                    out[i] = retemp + imtemp*I;
                }
                fclose(file);

            } 
            else 
            {
                //printf("## does not exist!\n");
                indices_to_run[num_to_run] = i;
                num_to_run ++;
            }
        }
        //end of loading saved GFs?
        
        if (num_to_run > 0)
        {
            printf("## new GFs to calculate: %d\n", num_to_run);
            a1s = createIntArray(num_to_run);
            a2s = createIntArray(num_to_run);
            ds = createIntArray(num_to_run);
            outs = createCompArray(num_to_run);
            
            for (i=0; i<num_to_run; i++)
            {
                a1s[i] = a1[indices_to_run[i]];
                a2s[i] = a2[indices_to_run[i]];
                ds[i] = diagt[indices_to_run[i]];
            }
            
            outs = graphenegfac(En, a1s, a2s, ds, num_to_run);
            
            for (i=0; i<num_to_run; i++)
            {
                out[indices_to_run[i]] = outs[i];
            }
            
            //printf("GFs made!\n");
            //save GFS to library
            for (i=0; i<num_to_run; i++)
            {
                sprintf(direcname, "%s%.0e", PGFlibloc, cimag(En));
                sprintf(command, "mkdir -p %s", direcname);
                //printf("%s\n", command);
                
                system(command);
                
                sprintf(fullname, "%s/%d_%d_%d_E_%+.8lf.dat", direcname, a1s[i], a2s[i], ds[i], creal(En));
                
                file = fopen(fullname, "w");
                fprintf(file, "%.10e    %.10e\n", creal(outs[i]), cimag(outs[i]));
                fclose(file);
                
                
                //out[indices_to_run[i]] = outs[i];
            }
            
            free(a1s); free(a2s); free(ds); free(outs);
        }
        

        free(indices_to_run); 
        return out;


}


    
    
//returns the vectorised integrand for the calculation of the graphene GF.
void intgac(unsigned ndim, const double *kx, void *p, unsigned fdim, double *fval)
{
        intgac_params *params = (intgac_params *)p;
        double *x = (params->x);
        double *y = (params->y);
        double _Complex En = (params->En);
        int *diagt = (params->diagt);
        double _Complex ky, ky1, ky2, kya, kyb;
        double _Complex cosky1, cosky2, sinky1, sinky2, cosky, S1, S2, Sa, Sb, Sda, Sdb, Sd1, Sd2, S, ans, ans1, ans2, ans3;
        double coskx; 
        double t=-1.0;
        double dis=1.0;
        int i;
        
       
        S1=0.0;
        S2=0.0;
        
        int dim = fdim/2;

        coskx = cos(sqrt(3) * kx[0] * dis / 2);
        cosky1 = (-0.5*(coskx + csqrt(En*En/(t*t) - 1 + coskx*coskx ) ));
        cosky2 = (-0.5*(coskx - csqrt(En*En/(t*t) - 1 + coskx*coskx ) ));
        
        for(i=0; i<dim; i++)
        {
        
            ans = 0.0;
            if(diagt[i] == 0)
            {
                    getS1S2a (&S1, &S2, x[i], y[i], kx[0], cosky1, cosky2, En);
                    ans = (S1 + S2)*En;
            }

            else if (diagt[i] == 2)
            {
                    getS1S2a (&S1, &S2, x[i], y[i], kx[0], cosky1, cosky2, En);
                    ans1 = (S1 + S2)*t;

                    getS1S2a (&S1, &S2, x[i], y[i] - dis, kx[0], cosky1, cosky2, En);
                    ans2 = (S1 + S2)*t;

                    getS1S2a (&S1, &S2, x[i] - sqrt(3)*dis/2, y[i] - dis/2, kx[0], cosky1, cosky2, En);
                    ans3 = (S1 + S2)*t;

                    ans = ans1+ans2+ans3;
            }
            else if (diagt[i] == 1)
            {
                    
                    getS1S2a (&S1, &S2, x[i], y[i], kx[0], cosky1, cosky2, En);
                    ans1 = (S1 + S2)*t;

                    getS1S2a (&S1, &S2, x[i], y[i] + dis, kx[0], cosky1, cosky2, En);
                    ans2 = (S1 + S2)*t;

                    getS1S2a (&S1, &S2, x[i] + sqrt(3)*dis/2, y[i] + dis/2, kx[0], cosky1, cosky2, En);
                    ans3 = (S1 + S2)*t;

                    
                    ans = ans1+ans2+ans3;
            }
            
            fval[2*i] = creal(ans);
            fval[2*i+1] = cimag(ans);
            
        }
}


//returns the integrand for the calculation of the graphene GF.
double intga(double kx, void *p)
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
        double t=-1.0;
        double dis=1.0;

        S1=0.0;
        S2=0.0;
        


        coskx = cos(sqrt(3) * kx * dis / 2);
        cosky1 = (-0.5*(coskx + csqrt(En*En/(t*t) - 1 + coskx*coskx ) ));
        cosky2 = (-0.5*(coskx - csqrt(En*En/(t*t) - 1 + coskx*coskx ) ));
        

        if(diagt == 0)
        {
                getS1S2a (&S1, &S2, x, y, kx, cosky1, cosky2, En);

                ans = (S1 + S2)*En;
        }

        else if (diagt == 2)
        {
                getS1S2a (&S1, &S2, x, y, kx, cosky1, cosky2, En);
                ans1 = (S1 + S2)*t;

                getS1S2a (&S1, &S2, x, y - dis, kx, cosky1, cosky2, En);
                ans2 = (S1 + S2)*t;

                getS1S2a (&S1, &S2, x - sqrt(3)*dis/2, y - dis/2, kx, cosky1, cosky2, En);
                ans3 = (S1 + S2)*t;

                ans = ans1+ans2+ans3;
        }
        else if (diagt == 1)
        {
                
                getS1S2a (&S1, &S2, x, y, kx, cosky1, cosky2, En);
                ans1 = (S1 + S2)*t;

                getS1S2a (&S1, &S2, x, y + dis, kx, cosky1, cosky2, En);
                ans2 = (S1 + S2)*t;

                getS1S2a (&S1, &S2, x + sqrt(3)*dis/2, y + dis/2, kx, cosky1, cosky2, En);
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
    
    
    
void getS1S2a (double _Complex *S1, double _Complex *S2, double x, double y, double kx, double _Complex cosky1, double _Complex cosky2, double _Complex En)
{
        double _Complex ky1, ky2, sinky1, sinky2, tempS1, tempS2;
        double dis=1.0;
        double t=-1.0;
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