#include "devices.h"
#include "disorder.h"
#include "../libs/matrices.h"
#include <math.h>
#include <time.h>
#include "useful.h"


//complex absorbing potentials at the top and bottom edges of nanoribbon devices.
//based on arxiv:1805.10220 (and ref[47] within)

void cap_potential(RectRedux *DeviceCell, double width)
{
    int Ntot = *(DeviceCell->Ntot);
    int Nrem = *(DeviceCell->Nrem);
    int length1 = (DeviceCell->length);
    int length2 = (DeviceCell->length2);
    int geo = (DeviceCell->geo);
    double topedge, bottomedge, topdist, botdist, xi, xf, x;
    (DeviceCell->cap_pots) = createDoubleArray(Ntot);
    double magn = (0.23 / width) * (4*M_PI*M_PI) * (4 / (2.62*2.62)); //0.23 = hbar^2 / 2 m renormalised for del x in a

    int i, j, k;
    
    //top and bottom of ribbon (offset)
         if(geo==0)
        {
            topedge=length1*sqrt(3)/2 - 1 / (2*sqrt(3)) +0.00001;
            bottomedge=1 / (2*sqrt(3)) - 0.00001;
        }
        if(geo==1 )
        {
            topedge=(length1-1)*0.5 +0.00001;
            bottomedge = -0.00001;
        }
    
    
    
    for(j=0; j<Ntot; j++)
    {
        if( (DeviceCell->siteinfo[j][0]) == 0)
        {
            //top edge
            if( (DeviceCell->pos)[j][1] > (topedge - width) )
            {
                x = (DeviceCell->pos)[j][1];
                xi = (topedge - width);
                xf = topedge;
                
                (DeviceCell->cap_pots)[j] = magn * ( pow(( width / (xf - 2*xi + x) ),2) + pow(( width / (xf - x) ),2) -2 );
            }
            
            //bottom edge
            if( (DeviceCell->pos)[j][1] < (bottomedge + width) )
            {
                x = (DeviceCell->pos)[j][1];
                xi = bottomedge + width;
                xf = bottomedge;
                
                (DeviceCell->cap_pots)[j] = magn * ( pow(( width / (xf - 2*xi + x) ),2) + pow(( width / (xf - x) ),2) -2 );
            }
        }
        //printf("%lf %lf\n", (DeviceCell->pos)[j][1], (DeviceCell->cap_pots)[j]);
    }
            
    
    
}

//absorbing potential on all 4 sides -- be VERY CAREFUL HERE!
void cap_potential2(RectRedux *DeviceCell, double width)
{
    int Ntot = *(DeviceCell->Ntot);
    int Nrem = *(DeviceCell->Nrem);
    int length1 = (DeviceCell->length);
    int length2 = (DeviceCell->length2);
    int geo = (DeviceCell->geo);
    double topedge, bottomedge, leftedge, rightedge, topdist, bottomdist, leftdist, rightdist, xi, xf, x;
    (DeviceCell->cap_pots) = createDoubleArray(Ntot);
    double magn = (0.23 / width) * (4*M_PI*M_PI) * (4 / (2.62*2.62)); //0.23 = hbar^2 / 2 m renormalised for del x in a
    
    
    int i, j, k;
    
    //top, bottom, left, right of ribbon (offset)
        if(geo==0)
        {
            topedge=length1*sqrt(3)/2 - 1 / (2*sqrt(3)) +0.00001;
            bottomedge=1 / (2*sqrt(3)) - 0.00001;
            leftedge= -0.00001;
            rightedge=length2*1.0 -0.5 +0.00001 ;
        }
        if(geo==1 )
        {
            topedge=(length1-1)*0.5 +0.00001;
            bottomedge = -0.00001;
            leftedge = - 0.00001;
            rightedge = length2*sqrt(3) - 1 / (sqrt(3)) +0.00001;
        }
    
    
    
    for(j=0; j<Ntot; j++)
    {
        if( (DeviceCell->siteinfo[j][0]) == 0)
        {
            topdist = topedge - (DeviceCell->pos)[j][1];
            bottomdist = (DeviceCell->pos)[j][1] - bottomedge;
            rightdist = rightedge - (DeviceCell->pos)[j][0];
            leftdist = (DeviceCell->pos)[j][0] - leftedge; 
            
            //top edge
            if( (topdist < width) && (topdist <= leftdist) && (topdist <= rightdist) )
            //if( (DeviceCell->pos)[j][1] > (topedge - width) )
            {
                x = (DeviceCell->pos)[j][1];
                xi = (topedge - width);
                xf = topedge;
                
                (DeviceCell->cap_pots)[j] = magn * ( pow(( width / (xf - 2*xi + x) ),2) + pow(( width / (xf - x) ),2) -2 );
            }
            
            //bottom edge
            //if( (DeviceCell->pos)[j][1] < (bottomedge + width) )
            if( (bottomdist < width) && (bottomdist <= leftdist) && (bottomdist <= rightdist) )
            {
                x = (DeviceCell->pos)[j][1];
                xi = bottomedge + width;
                xf = bottomedge;
                
                (DeviceCell->cap_pots)[j] = magn * ( pow(( width / (xf - 2*xi + x) ),2) + pow(( width / (xf - x) ),2) -2 );
            }
            
            
            //left edge
            if( (leftdist < width) && (leftdist < topdist) && (leftdist < bottomdist) )
            {
                x = (DeviceCell->pos)[j][0];
                xi = leftedge + width;
                xf = leftedge;
                
                (DeviceCell->cap_pots)[j] = magn * ( pow(( width / (xf - 2*xi + x) ),2) + pow(( width / (xf - x) ),2) -2 );
            }
            
            //right edge
            if( (rightdist < width) && (rightdist < topdist) && (rightdist < bottomdist) )
            {
                x = (DeviceCell->pos)[j][0];
                xi = rightedge - width;
                xf = rightedge;
                
                (DeviceCell->cap_pots)[j] = magn * ( pow(( width / (xf - 2*xi + x) ),2) + pow(( width / (xf - x) ),2) -2 );
            }
            
            
        }
        //printf("%lf %lf %lf\n", (DeviceCell->pos)[j][0], (DeviceCell->pos)[j][1], (DeviceCell->cap_pots)[j]);
    }
            
    
    
}




//form of general additional disorder routines: read in DeviceCell info and a disorder params structure and fileout info
//additional info written ADDED to what is already in site_pots and other DeviceCell arrays
//no way of adding offdiagonal disorder directly - this could possibly by varying atomic positions here and using a clever hopping rule at the GF stage


//this adds potential disorder in the form of electron-hole puddles or anderson-type disorder
//parameters are:
// 'conc'   - concentration of scattering centres
// 'delta'  - range of scatterer strength between - and + delta, in units of t0
// 'xi'     - range parameter determining gaussial falloff of parameter, in units of a


void potentialDisorder (RectRedux *DeviceCell, void *p, int disprof_out, char *filename )
{
  potDis_para *params = (potDis_para *)p;
  srand(time(NULL) + params->seed);

  int Ntot = *(DeviceCell->Ntot);
  int Nrem = *(DeviceCell->Nrem);
  int Nimp = (int)(Nrem * (params->conc));
  
  int i, j, k;
  int *possible_sites = createIntArray(Nrem);
  double *Uns = createDoubleArray(Nimp);
  double *dispots = createDoubleArray(Ntot);
  
  j=0;
  for(i=0; i<Ntot; i++)
  {
    if( (DeviceCell->siteinfo)[i][0] == 0)
    {
      possible_sites[j] = i;
      j++;
    }
  }
  
  for(i=0; i<Nrem; i++)
  {
    j=myRandInt(i, Nrem-1);
    k=possible_sites[i];
    possible_sites[i] = possible_sites[j];
    possible_sites[j] = k;
    
    if(i<Nimp)
    {
      Uns[i] = myRandNum(- (params->delta), (params->delta));
//       printf("%lf\n ", Uns[i]);

    }
  }
  
  double kappa=0, mfp;
  for(i=0; i<Nimp; i++)
  {
    for(j=0; j<Ntot; j++)
    {
      if( (DeviceCell->siteinfo[j][0]) == 0)
      {
	dispots[j] += Uns[i] * exp( -  ( pow((DeviceCell->pos[j][0]) - (DeviceCell->pos[possible_sites[i]][0]), 2) + pow((DeviceCell->pos[j][1]) - (DeviceCell->pos[possible_sites[i]][1]), 2)) / (2*(params->xi)*(params->xi)));
	kappa += (1.0/Nimp) * exp( -  ( pow((DeviceCell->pos[j][0]) - (DeviceCell->pos[possible_sites[i]][0]), 2) + pow((DeviceCell->pos[j][1]) - (DeviceCell->pos[possible_sites[i]][1]), 2)) / (2*(params->xi)*(params->xi)));
// 	 printf("#s: %e\n", kappa);
      }
    }
  }
    
  
  for(i=0; i<Ntot; i++)
  {
    (DeviceCell->site_pots[i]) += dispots[i];
  }
  
  
  FILE *out;
  if(disprof_out != 0)
  {
    out = fopen(filename, "w");
    
    for(i=0; i<Ntot; i++)
    {
      fprintf(out, "%lf %lf %e\n", (DeviceCell->pos)[i][0], (DeviceCell->pos)[i][1], dispots[i]);
    }
    
    fclose(out);
  }
  
  
  char filename_add[300];
  sprintf(filename_add, "%s.kap", filename);
  
  out = fopen(filename_add, "w");
  
  fprintf(out, "#kappa: %e\n", kappa);
  
  fclose(out);
  free(possible_sites);
  free(Uns);
  free(dispots);
  
  
}


//general edge potentials
void customEdgePots (RectRedux *DeviceCell, void *p)
{
  cedgepot_para *params = (cedgepot_para *)p;

  int Ntot = *(DeviceCell->Ntot);
  int Nrem = *(DeviceCell->Nrem);
  int length1 = (DeviceCell->length);
  int length2 = (DeviceCell->length2);
  int geo = (DeviceCell->geo);
  int i, j, k;
  double *dispots = createDoubleArray(Ntot);
  int type = params->type;
  
  double topedge, bottomedge, topdist, botdist, ribcenter;
  double subs_eps=3.9;
  double subs_thick=0.1E-6, n0, nmax, ny, gate_voltage;
  double edge_cut=5.0;
  
	//these are the exact positions of the top and bottom atomic sites
	if(geo==0)
	{
		topedge=length1*sqrt(3)/2 - 1.0/(2*sqrt(3));
		bottomedge=1.0/(2*sqrt(3));
		
		//efetov type
		if(type == 4)
		{
			topedge=length1*sqrt(3)/2;
			bottomedge=0.0;
			ribcenter = topedge/2;
		}
	}

	//armchair offset slightly to avoid division by zero for very edge sites
	if(geo==1 )
	{
		topedge=(length1 -1)*0.5;
		bottomedge = 0.0;
		
		//efetov type
		if(type == 4)
		{
			topedge=(length1)*0.5;
			bottomedge = -0.5;
			ribcenter = topedge/2 -0.25;
		}
	}

  
 
    for(j=0; j<Ntot; j++)
    {
      if( (DeviceCell->siteinfo[j][0]) == 0 && (DeviceCell->siteinfo[j][1]) == 0)
      {
	      
	if(type==0)
	{
	  if( (topedge - (params->AT3)) - (DeviceCell->pos)[j][1]  < (params->AT2) )
	  {
		if ((topedge - (params->AT3)) - (DeviceCell->pos)[j][1]  > 0 )
			dispots[j] += (params->AT1);
	  }
	  if( ((DeviceCell->pos)[j][1] - (params->AB3) - bottomedge) < (params->AB2) )
	  {	  
		if( ((DeviceCell->pos)[j][1] - (params->AB3) - bottomedge) > 0 )
		dispots[j] += (params->AB1);
		
	  }
	}
	
	if(type==1)
	{
  	  dispots[j] += (params->AT1)*exp(- (params->AT2) * fabs((DeviceCell->pos)[j][1] + (params->AT3) - topedge) ) +  (params->AB1)*exp(- (params->AB2) * fabs((DeviceCell->pos)[j][1] - (params->AB3) - bottomedge) );
	}
	
	if(type==2)
	{
  	  dispots[j] += (params->AT1)*pow(fabs( ((DeviceCell->pos)[j][1] + (params->AT3) - topedge)), -(params->AT2)) + (params->AB1)*pow(fabs( ((DeviceCell->pos)[j][1] - (params->AB3) - bottomedge)), -(params->AB2)) ;
	}
	
	if(type==3)
	{
		if((DeviceCell->pos)[j][1] > (bottomedge + (params->AB3)) && (DeviceCell->pos)[j][1] < (topedge - (params->AT3)))
		{
			dispots[j] += (params->AB1) + ((params->AT1)-(params->AB1)) * ( (DeviceCell->pos)[j][1] - (params->AB3) - bottomedge) / (topedge - (params->AT3) - (params->AB3) - bottomedge);
		}
		if((DeviceCell->pos)[j][1] < (bottomedge + (params->AB3)))
		{
			dispots[j] += (params->AB1);
		}
		if((DeviceCell->pos)[j][1] > (topedge - (params->AT3)))
		{
			dispots[j] += (params->AT1);

		}
		
	}
	
      }
      
      if( (DeviceCell->siteinfo[j][0]) == 0 && (DeviceCell->siteinfo[j][1]) == 1)
      {
	      
	if(type==0)
	{
	  if( (topedge - (params->BT3)) - (DeviceCell->pos)[j][1]  < (params->BT2) )
	  {
		if( (topedge - (params->BT3)) - (DeviceCell->pos)[j][1] > 0 )  
			dispots[j] += (params->BT1);
	  }
	  if( ((DeviceCell->pos)[j][1] - (params->BB3) - bottomedge) < (params->BB2) )
	  {
		if( ((DeviceCell->pos)[j][1] - (params->BB3) - bottomedge) > 0 )
			dispots[j] += (params->BB1);
	  }
	}
	
	if(type==1)
	{
  	  dispots[j] += (params->BT1)*exp(- (params->BT2) * fabs((DeviceCell->pos)[j][1] + (params->BT3) - topedge) ) +  (params->BB1)*exp(- (params->BB2) * fabs((DeviceCell->pos)[j][1] - (params->BB3) - bottomedge) );
	}
	
	if(type==2)
	{
  	  dispots[j] += (params->BT1)*pow(fabs( ((DeviceCell->pos)[j][1] + (params->BT3) - topedge)), -(params->BT2)) + (params->BB1)*pow(fabs( ((DeviceCell->pos)[j][1] - (params->BB3) - bottomedge)), -(params->BB2)) ;
	}
	
	if(type==3)
	{
	  //dispots[j] += (params->BB1) + ((params->BT1)-(params->BB1)) * ( (DeviceCell->pos)[j][1] - (params->BB3) - bottomedge) / (topedge - (params->BT3) - (params->BB3) - bottomedge);
	  
		if((DeviceCell->pos)[j][1] > (bottomedge + (params->BB3)) && (DeviceCell->pos)[j][1] < (topedge - (params->BT3)))
		{
			dispots[j] += (params->BB1) + ((params->BT1)-(params->BB1)) * ( (DeviceCell->pos)[j][1] - (params->BB3) - bottomedge) / (topedge - (params->BT3) - (params->BB3) - bottomedge);
		}
		if((DeviceCell->pos)[j][1] < (bottomedge + (params->BB3)))
		{
			dispots[j] += (params->BB1);
		}
		if((DeviceCell->pos)[j][1] > (topedge - (params->BT3)))
		{
			dispots[j] += (params->BT1);
		}
	  
	  
	  
	  
	}
	
      }
      
      
      //basic efefov: AT1=VG, AT2=edge_cut, AT3=subs_thick, subs_eps=3.9 (SiO2)
	if(type==4)
	{
		//allow non-default values of edge_cut and subs_thick
		if( (params->AT2) != 0.0)
		{
			edge_cut = (params->AT2);
		}
		
		if( (params->AT3) != 0.0)
		{
			subs_thick = (params->AT3);
		}
		
		gate_voltage = (params->AT1);
		
		n0= (subs_eps)*(8.85E-12)*fabs(gate_voltage)/(2*subs_thick * 1.6E-19);
		nmax = n0 * ((topedge-bottomedge)/2) / sqrt( pow( (topedge-bottomedge)/2, 2)  - pow((topedge-ribcenter) - edge_cut, 2) );
		
		
		ny = n0 * ((topedge-bottomedge)/2)  / sqrt( pow((topedge-bottomedge)/2, 2)  - pow((DeviceCell->pos)[j][1] - ribcenter, 2) );
	
		if(ny>nmax)
			ny=nmax;
		
		
		if(gate_voltage != 0.0)
		{
			dispots[j] += (- gate_voltage/ fabs(gate_voltage)) * (1.05E-28) * (ny/fabs(ny))* sqrt(M_PI * fabs(ny) ) / (2.7 * 1.6E-19) ;
		}
		
	
	}
      
      
      
      
    }
  
    
  
  for(i=0; i<Ntot; i++)
  {
    (DeviceCell->site_pots[i]) += dispots[i];
  }
  
  
  
  free(dispots);
  
  
}






