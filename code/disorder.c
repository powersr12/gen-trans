#include "disorder.h"


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
  
  double topedge, bottomedge, topdist, botdist;
  
	//these are the exact positions of the top and bottom atomic sites
	if(geo==0)
	{
		topedge=length1*sqrt(3)/2 - 1.0/(2*sqrt(3));
		bottomedge=1.0/(2*sqrt(3));
	}

	//armchair offset slightly to avoid division by zero for very edge sites
	if(geo==1 )
	{
		topedge=(length1 -1)*0.5;
		bottomedge = 0.0;
	}

  
 
    for(j=0; j<Ntot; j++)
    {
      if( (DeviceCell->siteinfo[j][0]) == 0 && (DeviceCell->siteinfo[j][1]) == 0)
      {
	      
	if(type==0)
	{
	  if( fabs((DeviceCell->pos)[j][1] + (params->AT3) - topedge) < (params->AT2) )
		dispots[j] += (params->AT1);
	  
	  if( fabs((DeviceCell->pos)[j][1] - (params->AB3) - bottomedge) < (params->AB2) )
		dispots[j] += (params->AB1);
	}
	
	if(type==1)
	{
  	  dispots[j] += (params->AT1)*exp(- (params->AT2) * fabs((DeviceCell->pos)[j][1] + (params->AT3) - topedge) ) +  (params->AB1)*exp(- (params->AB2) * fabs((DeviceCell->pos)[j][1] - (params->AB3) - bottomedge) );
	}
	
	if(type==2)
	{
  	  dispots[j] += (params->AT1)*pow(fabs( ((DeviceCell->pos)[j][1] + (params->AT3) - topedge)), -(params->AT2)) + (params->AB1)*pow(fabs( ((DeviceCell->pos)[j][1] - (params->AB3) - bottomedge)), -(params->AB2)) ;
	}
	
      }
      
      if( (DeviceCell->siteinfo[j][0]) == 0 && (DeviceCell->siteinfo[j][1]) == 1)
      {
	      
	if(type==0)
	{
	  if( fabs((DeviceCell->pos)[j][1] + (params->BT3) - topedge) < (params->BT2) )
		dispots[j] += (params->BT1);
	  
	  if( fabs((DeviceCell->pos)[j][1] - (params->BB3) - bottomedge) < (params->BB2) )
		dispots[j] += (params->BB1);
	}
	
	if(type==1)
	{
  	  dispots[j] += (params->BT1)*exp(- (params->BT2) * fabs((DeviceCell->pos)[j][1] + (params->BT3) - topedge) ) +  (params->BB1)*exp(- (params->BB2) * fabs((DeviceCell->pos)[j][1] - (params->BB3) - bottomedge) );
	}
	
	if(type==2)
	{
  	  dispots[j] += (params->BT1)*pow(fabs( ((DeviceCell->pos)[j][1] + (params->BT3) - topedge)), -(params->BT2)) + (params->BB1)*pow(fabs( ((DeviceCell->pos)[j][1] - (params->BB3) - bottomedge)), -(params->BB2)) ;
	}
	
      }
      
      
      
      
    }
  
    
  
  for(i=0; i<Ntot; i++)
  {
    (DeviceCell->site_pots[i]) += dispots[i];
  }
  
  
  
  free(dispots);
  
  
}






