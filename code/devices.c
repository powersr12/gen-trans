

#include "devices.h"


//Routines for generating various devices (RectRedux form)

//creates the simplest RectRedux array for a finite number of ribbon chains
//note that this format should be the general layout for "generateDevice" routines
//with other disorder or antidot parameters placed into parameter structure p 
void **simpleRibbonGeo (RectRedux *SiteArray, void *p, int struc_out, char *filename)
{
	  int length = SiteArray->length;
	  int length2 = SiteArray->length2;
	  int geo = SiteArray->geo;
	    
	  double smalldist;
	  FILE *out;
	  
	  if(struc_out != 0)
	  {
	    out = fopen(filename, "w");
	    
	  }
	    
	 // srand(time(NULL) + seed);
	  

	  //atomic coordinates and the atoms that are in and out, similar to previous cases
	  int tot_sites = 2*length*length2;
	  double **site_coords = createNonSquareDoubleMatrix(tot_sites, 3);
	  double *site_pots = createDoubleArray(tot_sites);
	  int **siteinfo = createNonSquareIntMatrix(tot_sites, 2);
	  double xstart, ystart;
	  int isclean, l, m;
	  int *Nrem = (SiteArray->Nrem);
	  int *Ntot = (SiteArray->Ntot);
	  
	  *Ntot = tot_sites;

	  if(geo==0)
	  {
	    for(l=0; l<length2; l++)
	    {
	      xstart=l*1.0;
	      for(m=0; m<length; m++)
	      {
		ystart= m*sqrt(3)/2 + 1/(2*sqrt(3));
		
		if((m%2) == 0)
		{
		    site_coords[l*2*length + 2*m][0] = xstart+0.5;
		    site_coords[l*2*length + 2*m + 1][0] = xstart;
		}
		
		if((m%2) == 1)
		{
		    site_coords[l*2*length + 2*m][0] = xstart;
		    site_coords[l*2*length + 2*m + 1][0] = xstart+0.5;
		}
		
		    site_coords[l*2*length + 2*m][1] = ystart;
		    site_coords[l*2*length + 2*m + 1][1] = ystart + 1/(2*sqrt(3));
		    siteinfo[l*2*length + 2*m][1]=0;
		    siteinfo[l*2*length + 2*m +1][1]=1;

		  
	      }
	    }
	  }
	  
	  if(geo==1)
	  {
	      for(l=0; l<length2; l++)
	      {
		xstart = l*sqrt(3);
		for(m=0; m<length; m++)
		{
		  if((m%2) == 0)
		  {
		    site_coords[l*2*length + m][0] = xstart;
		    site_coords[l*2*length + length + m][0] = xstart +2 / sqrt(3);
		    siteinfo[l*2*length + m][1] = 0;
		    siteinfo[l*2*length + length + m][1] = 1;
		  }
		  if((m%2) == 1)
		  {
		    site_coords[l*2*length + m][0] = xstart +  1/(2*sqrt(3));
		    site_coords[l*2*length + length + m][0] = xstart + (sqrt(3))/2;
		    siteinfo[l*2*length + m][1] = 1;
		    siteinfo[l*2*length + length + m][1] = 0;
		  }
		  
		  site_coords[l*2*length + m][1] = m*0.5;
		  site_coords[l*2*length + length + m][1] = m*0.5;

		  
		}
		
		
	      }
		
	  }
	      
	          
	      	      
	      //chaininfo needed for conductance calcs (if atoms are missing)
		      (SiteArray->chaininfo) = createNonSquareIntMatrix(length2, 4);
		      int tempint=0, tempint2=0;
	 
			for(l=0; l<length2; l++)
			{
			  tempint=0;
			  for(m=0; m<2*length; m++)
			  {
			    if(siteinfo[l*2*length +m][0] == 0)
			    {
			      tempint ++;
			      tempint2++;
			    }
			  }
			  (SiteArray->chaininfo)[l][0] = tempint;
			}
			*Nrem = tempint2;
			
			
			(SiteArray->chaininfo)[0][1] = 0;
			for(l=1; l<length2; l++)
			{
			  (SiteArray->chaininfo)[l][1] = (SiteArray->chaininfo)[l-1][1] + (SiteArray->chaininfo)[l-1][0];
			}


			//are first and last atoms in each chain present?
			for(l=1; l<length2; l++)
			{
			  if(siteinfo[l*2*length][0] == 1)
			    (SiteArray->chaininfo)[l][2] = 1;
			  
			  if(siteinfo[(l+1)*2*length -1][0] == 1)
			    (SiteArray->chaininfo)[l][3] = 1;
			}
	  
	  
	  
		      
			if(struc_out == 1)
			{
			  for(l=0; l<2*length*length2; l++)
			  {
			    if(siteinfo[l][0] == 0 )
			    {  
			      fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			    }
			  }
			  fprintf(out, "\n");
			  

			  
			}
			
			  if(struc_out != 0)
			  {
			    fclose(out);
			  }
			    
			
			  
			  //Fill the array of data structures describing the system
			  (SiteArray->pos) = site_coords;
			  (SiteArray->site_pots) = site_pots;
			  (SiteArray->siteinfo) = siteinfo;
  
  
}

//generate lead geometries and rect_reduxes to match a device
//mode=0 - basic left-right setup, same size 
void genLeads (RectRedux *SiteArray, RectRedux **Leads, int numleads, int mode, lead_para *params)
{
  int i, j;
  if(mode==0)
  {
   
    
    for(i=0; i < numleads; i++)
    {
      Leads[i] = (RectRedux *)malloc(sizeof(RectRedux));
      (Leads[i]->geo) = (SiteArray->geo);
      (Leads[i]->length) = (SiteArray->length);
      (Leads[i]->length2) = 1;
      (Leads[i]->Nrem) = (int *)malloc(sizeof(int));
      (Leads[i]->Ntot) = (int *)malloc(sizeof(int));
      simpleRibbonGeo (Leads[i], NULL, 0, NULL);
    }
    
    
    (params->shift_vecs) = createNonSquareDoubleMatrix(numleads, 3);
    (params->shift_vecs)[0][0] = -((SiteArray->pos)[2*(SiteArray->length)][0] - (SiteArray->pos)[0][0]);
    (params->shift_vecs)[1][0] = (SiteArray->pos)[2*(SiteArray->length)][0] - (SiteArray->pos)[0][0];
    (params->shift_vecs)[0][1] = 0;
    (params->shift_vecs)[1][1] = 0;
    
    
    
    //shift leads accordingly to desired positions
    for(i=0; i<*(Leads[0]->Ntot); i++)
    {
      (Leads[0]->pos)[i][0] += (params->shift_vecs)[0][0];      
    }
    
    for(i=0; i<*(Leads[1]->Ntot); i++)
    {
      (Leads[1]->pos)[i][0] += (SiteArray->length2)*(params->shift_vecs)[1][0];      
    }
    
  }
  
// 	for(i=0; i<numleads; i++)
// 	{
// 	  for(j=0; j<*(Leads[i]->Ntot); j++)
// 	  {
// 	    printf("%lf	%lf\n", (Leads[i]->pos)[j][0], (Leads[i]->pos)[j][1]);
// 	  }
// 	  printf("\n");
// 	}
  
  
  
}

//very general - creates ribbon with sublattice dependent potential
void genSublatticeDevice (RectRedux *SiteArray, int buffer_rows, double a_conc, double a_pot, double b_conc, double b_pot, int seed, int struc_out, char *filename)
{  
	  
	  int length = SiteArray->length;
	    int length2 = SiteArray->length2;
	    int geo = SiteArray->geo;
	    
	  double smalldist;
	  FILE *out;
	  
	  if(struc_out != 0)
	  {
	    out = fopen(filename, "w");
	    
	  }
	    
	  srand(time(NULL) + seed);
	  
	  
	  //atomic coordinates and the atoms that are in and out, similar to previous cases
	  int tot_sites = 2*length*length2;
	  double **site_coords = createNonSquareDoubleMatrix(tot_sites, 3);
	  double *site_pots = createDoubleArray(tot_sites);
	  int **siteinfo = createNonSquareIntMatrix(tot_sites, 2);
	  double xstart, ystart;
	  int isclean, l, m;
	  int *Nrem = (SiteArray->Nrem);
	  int *Ntot = (SiteArray->Ntot);
	  
	  *Ntot = tot_sites;
	  
	  if(geo==0)
	  {
	    for(l=0; l<length2; l++)
	    {
	      xstart=l*1.0;
	      for(m=0; m<length; m++)
	      {
		ystart= m*sqrt(3)/2 + 1/(2*sqrt(3));
		
		if((m%2) == 0)
		{
		    site_coords[l*2*length + 2*m][0] = xstart+0.5;
		    site_coords[l*2*length + 2*m + 1][0] = xstart;
		}
		
		if((m%2) == 1)
		{
		    site_coords[l*2*length + 2*m][0] = xstart;
		    site_coords[l*2*length + 2*m + 1][0] = xstart+0.5;
		}
		
		    site_coords[l*2*length + 2*m][1] = ystart;
		    site_coords[l*2*length + 2*m + 1][1] = ystart + 1/(2*sqrt(3));
		    siteinfo[l*2*length + 2*m][1]=0;
		    siteinfo[l*2*length + 2*m +1][1]=1;

		  
	      }
	    }
	  }
	  
	  if(geo==1)
	  {
	      for(l=0; l<length2; l++)
	      {
		xstart = l*sqrt(3);
		for(m=0; m<length; m++)
		{
		  if((m%2) == 0)
		  {
		    site_coords[l*2*length + m][0] = xstart;
		    site_coords[l*2*length + length + m][0] = xstart +2 / sqrt(3);
		    siteinfo[l*2*length + m][1] = 0;
		    siteinfo[l*2*length + length + m][1] = 1;
		  }
		  if((m%2) == 1)
		  {
		    site_coords[l*2*length + m][0] = xstart +  1/(2*sqrt(3));
		    site_coords[l*2*length + length + m][0] = xstart + (sqrt(3))/2;
		    siteinfo[l*2*length + m][1] = 1;
		    siteinfo[l*2*length + length + m][1] = 0;
		  }
		  
		  site_coords[l*2*length + m][1] = m*0.5;
		  site_coords[l*2*length + length + m][1] = m*0.5;

		  
		}
		
		
	      }
		
	  }
	      
	      //disordery stuff here
	      double temprandnum;
	      for(l=buffer_rows*2*length; l<2*length*length2 - 2*buffer_rows*length; l++)
	      {
		site_pots[l] = 0.0;
		temprandnum = myRandNum(0.0, 1.0);
		
		if(siteinfo[l][1] == 0)
		{
		  if(temprandnum < a_conc)
		    site_pots[l] = a_pot;
		}
		
		if(siteinfo[l][1] == 1)
		{
		  if(temprandnum < b_conc)
		    site_pots[l] = b_pot;
		}
		//printf("%d	%lf	%d	%lf\n", l, temprandnum, siteinfo[l][1], site_pots[l]);
	      }
	      
	      	      
	      //chaininfo needed for conductance calcs (if atoms are missing)
		      (SiteArray->chaininfo) = createNonSquareIntMatrix(length2, 4);
		      int tempint=0, tempint2=0;
	 
			for(l=0; l<length2; l++)
			{
			  tempint=0;
			  for(m=0; m<2*length; m++)
			  {
			    if(siteinfo[l*2*length +m][0] == 0)
			    {
			      tempint ++;
			      tempint2++;
			    }
			  }
			  (SiteArray->chaininfo)[l][0] = tempint;
			}
			*Nrem = tempint2;
			
			
			(SiteArray->chaininfo)[0][1] = 0;
			for(l=1; l<length2; l++)
			{
			  (SiteArray->chaininfo)[l][1] = (SiteArray->chaininfo)[l-1][1] + (SiteArray->chaininfo)[l-1][0];
			}


			//are first and last atoms in each chain present?
			for(l=1; l<length2; l++)
			{
			  if(siteinfo[l*2*length][0] == 1)
			    (SiteArray->chaininfo)[l][2] = 1;
			  
			  if(siteinfo[(l+1)*2*length -1][0] == 1)
			    (SiteArray->chaininfo)[l][3] = 1;
			}
	  
	  
	  
		      
			if(struc_out == 1)
			{
			  for(l=0; l<2*length*length2; l++)
			  {
			    if(siteinfo[l][0] == 0 && siteinfo[l][1] == 0 && site_pots[l] == 0.0)
			    {  
			      fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			    }
			  }
			  fprintf(out, "\n");
			  
			  for(l=0; l<2*length*length2; l++)
			  {
			    if(siteinfo[l][0] == 0 && siteinfo[l][1] == 1 && site_pots[l] == 0.0)
			    {  
			      fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			    }
			  }
			  fprintf(out, "\n");
			  
			  if(a_pot!=0.0)
			  {
			    for(l=0; l<2*length*length2; l++)
			    {
			      if(siteinfo[l][0] == 0  && site_pots[l] == a_pot)
			      {  
				fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			      }
			    }
			  }
			  fprintf(out, "\n");
			  
			  if(b_pot!=0.0)
			  {
			    for(l=0; l<2*length*length2; l++)
			    {
			      if(siteinfo[l][0] == 0  && site_pots[l] == b_pot)
			      {  
				fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			      }
			    }
			  }
			  fprintf(out, "\n");
			  
			  
			}
			
			  if(struc_out != 0)
			  {
			    fclose(out);
			  }
			    
			
			  
			  //Fill the array of data structures describing the system
			  (SiteArray->pos) = site_coords;
			  (SiteArray->site_pots) = site_pots;
			  (SiteArray->siteinfo) = siteinfo;

	
}


//A general (non-disordered) antidot barrier-type device 
//(circular ALs in triangular lattice for the moment - should be generalised later)
//Different to the routine used in the disordered antidot paper
//Modified to work with generic ribbons devices
//Layout based on sublattice device routine
//The first and last buffer_rows chains of the device will not be altered from pristine graphene
//This version requires more aligning of the graphene and antidot sheets
//i.e. ribbon width is not automatically an integer multiple of GAL cell width
void genAntidotDevice (RectRedux *SiteArray, int buffer_rows, int AD_length, double AD_rad, int lat_width, int lat_length, int seed, int struc_out, char *filename)
{  
	  
	  int length = SiteArray->length;
	  int length2 = SiteArray->length2;
	  int geo = SiteArray->geo;
	    
	  double smalldist;
	  FILE *out;
	  
	  if(struc_out != 0)
	  {
	    out = fopen(filename, "w");
	    
	  }
	    
	  srand(time(NULL) + seed);
	  
	  
	  //atomic coordinates and the atoms that are in and out, similar to previous cases
	  int tot_sites = 2*length*length2;
	  double **site_coords = createNonSquareDoubleMatrix(tot_sites, 3);
	  double *site_pots = createDoubleArray(tot_sites);
	  int **siteinfo = createNonSquareIntMatrix(tot_sites, 2);
	  double xstart, ystart;
	  int isclean, l, m, tempint, tempint2;
	  int *Nrem = (SiteArray->Nrem);
	  int i, j, k;
	  int *Ntot = (SiteArray->Ntot);
	  
	  *Ntot = tot_sites;
	  
	  if(geo==0)
	  {
	    for(l=0; l<length2; l++)
	    {
	      xstart=l*1.0;
	      for(m=0; m<length; m++)
	      {
		ystart= m*sqrt(3)/2 + 1/(2*sqrt(3));
		
		if((m%2) == 0)
		{
		    site_coords[l*2*length + 2*m][0] = xstart+0.5;
		    site_coords[l*2*length + 2*m + 1][0] = xstart;
		}
		
		if((m%2) == 1)
		{
		    site_coords[l*2*length + 2*m][0] = xstart;
		    site_coords[l*2*length + 2*m + 1][0] = xstart+0.5;
		}
		
		    site_coords[l*2*length + 2*m][1] = ystart;
		    site_coords[l*2*length + 2*m + 1][1] = ystart + 1/(2*sqrt(3));
		    siteinfo[l*2*length + 2*m][1]=0;
		    siteinfo[l*2*length + 2*m +1][1]=1;

		  
	      }
	    }
	  }
	  
	  if(geo==1)
	  {
	      for(l=0; l<length2; l++)
	      {
		xstart = l*sqrt(3);
		for(m=0; m<length; m++)
		{
		  if((m%2) == 0)
		  {
		    site_coords[l*2*length + m][0] = xstart;
		    site_coords[l*2*length + length + m][0] = xstart +2 / sqrt(3);
		    siteinfo[l*2*length + m][1] = 0;
		    siteinfo[l*2*length + length + m][1] = 1;
		  }
		  if((m%2) == 1)
		  {
		    site_coords[l*2*length + m][0] = xstart +  1/(2*sqrt(3));
		    site_coords[l*2*length + length + m][0] = xstart + (sqrt(3))/2;
		    siteinfo[l*2*length + m][1] = 1;
		    siteinfo[l*2*length + length + m][1] = 0;
		  }
		  
		  site_coords[l*2*length + m][1] = m*0.5;
		  site_coords[l*2*length + length + m][1] = m*0.5;

		  
		}
		
		
	      }
		
	  }
	  
	  //Antidot positions
	  int tot_holes = 2*lat_width*lat_length;
	  double holes[tot_holes][3];
	  double unitholes[2][2];
	  double xshift, yshift;
	  
	  //zigzag ribbon
	  if(geo==0)
	  {
	    xstart = buffer_rows*1.0;
	    ystart = 0.0;
	    xshift = 3.0*AD_length;
	    yshift = sqrt(3)*AD_length;
	    
	    unitholes[0][0]  = (int) (AD_length*3.0/4) -0.5* ( ((int) (AD_length/2)) % 2 ); 
	    unitholes[0][1] = (sqrt(3)/2.0) * (int) (AD_length /2.0)   ;
	    unitholes[1][0] = unitholes[0][0] + 1.5*AD_length;
	    unitholes[1][1] = unitholes[0][1] + sqrt(3)*AD_length/2.0;
	  }
	  
	  
	  //armchair ribbon - test at some point
	  if(geo==1)
	  {
	    xstart = buffer_rows*sqrt(3);
	    ystart = 0.0;
	    xshift = sqrt(3)*AD_length;
	    yshift = 3.0*AD_length;
	    
// 	    unitholes[0][1]  = (int) (AD_length*3.0/4) -0.5* ( ((int) (AD_length/2)) % 2 ); 
// 	    unitholes[0][0] = (sqrt(3)/2.0) * (int) (AD_length /2.0)   ;
// 	    unitholes[1][1] = unitholes[0][0] + 1.5*AD_length;
// 	    unitholes[1][0] = unitholes[0][1] + sqrt(3)*AD_length/2.0;
	    
	    unitholes[0][1] =3.0*AD_length - 0.5 - ((int) (AD_length*3.0/4) -0.5* ( ((int) (AD_length/2)) % 2 ));
            unitholes[0][0] = (sqrt(3)/2.0) * (int) (AD_length /2.0)  - 1 /(2*sqrt(3)) ;
            unitholes[1][1] = unitholes[0][1] - 1.5*AD_length;
            unitholes[1][0] = unitholes[0][0] + sqrt(3)*AD_length/2.0;

	    
	    
	  }
	  
	  
	  for(i=0; i<lat_length; i++)
	  {
	    for(j=0; j<lat_width; j++)
	    {
	      holes[2*i*lat_width + 2*j][0] = xstart + i*xshift + unitholes[0][0];
	      holes[2*i*lat_width + 2*j][1] = ystart + j*yshift + unitholes[0][1];
	      holes[2*i*lat_width + 2*j][2] = AD_rad;
	       
	      holes[2*i*lat_width + 2*j + 1][0] = xstart + i*xshift + unitholes[1][0];
	      holes[2*i*lat_width + 2*j + 1][1] = ystart + j*yshift + unitholes[1][1];
	      holes[2*i*lat_width + 2*j + 1][2] = AD_rad;
	    }
	  }
	  
	  
	  int **sites = createNonSquareIntMatrix(tot_sites, 2); //removed, num neighbours after first sweep
	   
	  //atom removal!
	      for(i=2*length; i< tot_sites - 2*length ; i++)
	      {
		for(j=0; j< tot_holes; j++)
		{
		    if(sqrt( pow( site_coords[i][0] - holes[j][0], 2) + pow( site_coords[i][1] - holes[j][1], 2)) < holes[j][2])
		    {
		      sites[i][0] = 1;
		    }
		}
	      }
	  
	      //count neighbours
	      for(i=0; i<tot_sites; i++)
	      {
		if(sites[i][0] == 0)
		{
		  tempint=0;
		  for(j=0 ; j< tot_sites ;j++)
		  {
		    if(sites[j][0] == 0)
		    {
		      if(i!=j)
		      {
			smalldist = sqrt( pow(site_coords[i][0] - site_coords[j][0], 2) +  pow(site_coords[i][1] - site_coords[j][1], 2));
		      
		      
			if (smalldist < 0.6)
			{
			  tempint++;
			}
		      }
		    }
		  }
		  sites[i][1] =tempint;
		}
		 
	      }
	      
	     
// 
// 		  //remove relevant atoms from lattice  
		    for(i=2*length; i<tot_sites-2*length; i++)
		    {
		      if(sites[i][1] < 2)
		      {
			  sites[i][0] = 1;
		      }
		    }
		
	  
		    for(i=0; i<tot_sites; i++)
		    {
		      siteinfo[i][0] = sites[i][0];
		    }
	  
		    free(sites[0]);
		    free(sites);
	      

		    
	      //chaininfo needed for conductance calcs (if atoms are missing)
		      (SiteArray->chaininfo) = createNonSquareIntMatrix(length2, 4);
		      
		      tempint=0, tempint2=0;
	 
		     
			for(l=0; l<length2; l++)
			{
			  tempint=0;
			  for(m=0; m<2*length; m++)
			  {
			    if(siteinfo[l*2*length +m][0] == 0)
			    {
			      tempint ++;
			      tempint2++;
			    }
			  }
			  (SiteArray->chaininfo)[l][0] = tempint;
			}
			*Nrem = tempint2;
			
			
			(SiteArray->chaininfo)[0][1] = 0;
			for(l=1; l<length2; l++)
			{
			  (SiteArray->chaininfo)[l][1] = (SiteArray->chaininfo)[l-1][1] + (SiteArray->chaininfo)[l-1][0];
			}


			//are first and last atoms in each chain present?
			for(l=1; l<length2; l++)
			{
			  if(siteinfo[l*2*length][0] == 1)
			    (SiteArray->chaininfo)[l][2] = 1;
			  
			  if(siteinfo[(l+1)*2*length -1][0] == 1)
			    (SiteArray->chaininfo)[l][3] = 1;
			}
	  
	  
	  
		      
			if(struc_out == 1)
			{
			  for(l=0; l<2*length*length2; l++)
			  {
			    if(siteinfo[l][0] == 0)
			      fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			    
			  }
			  
			}
			  
				    
			    
			
			  
			  //Fill the array of data structures describing the system
			
			  
			  if(struc_out != 0)
			  {
			    fclose(out);
			    
			  }
			  
			  
			  (SiteArray->pos) = site_coords;
			  (SiteArray->site_pots) = site_pots;
			  (SiteArray->siteinfo) = siteinfo;


				
	
}


void exportRectConf(RectRedux *System, char *filename)
{
  //printf("%s\n", filename);
  FILE *fileout;
  char fullname[128];
  int length = System->length;
  int length2 = System->length2;
  int geo = System->geo;
  int i, j, k;
  int tot_sites = 2*length*length2;
  
  
  sprintf(fullname, "%s.siteinfo", filename);
  fileout = fopen(fullname, "w");
    for(j=0; j<tot_sites; j++)
    {
      fprintf(fileout, "%d	%d\n", (System->siteinfo)[j][0], (System->siteinfo)[j][1]);
    }
  fclose(fileout);
  
  sprintf(fullname, "%s.pos", filename);
  fileout = fopen(fullname, "w");
    for(j=0; j<tot_sites; j++)
    {
      fprintf(fileout, "%lf	%lf\n", (System->pos)[j][0], (System->pos)[j][1]);
    }
  fclose(fileout);
  
  
  sprintf(fullname, "%s.site_pots", filename);
  fileout = fopen(fullname, "w");
    for(j=0; j<tot_sites; j++)
    {
      fprintf(fileout, "%.10lf\n", (System->site_pots)[j]);
    }
  fclose(fileout);
  
  
  sprintf(fullname, "%s.chaininfo", filename);
  fileout = fopen(fullname, "w");
    for(j=0; j<length2; j++)
    {
      fprintf(fileout, "%d	%d	%d	%d\n", (System->chaininfo)[j][0], (System->chaininfo)[j][1], (System->chaininfo)[j][2], (System->chaininfo)[j][3]);
    }
  fclose(fileout);
  
  
  sprintf(fullname, "%s.info", filename);
  fileout = fopen(fullname, "w");
       fprintf(fileout, "%d\n", *(System->Nrem));
  fclose(fileout);

}


void importRectConf(RectRedux *System, int length, int length2, char *filename)
{
  //printf("%s	%d\n", filename, numcells);
  FILE *fileout;
  char fullname[128];
 // int length = System->length;
 // int length2 = System->length2;
 // int geo = System->geo;
  int i, j, k, temp;
  int tot_sites = 2*length*length2;
  
	  System->pos = createNonSquareDoubleMatrix(tot_sites, 3);
 	  System->site_pots = createDoubleArray(tot_sites);
	  System->siteinfo = createNonSquareIntMatrix(tot_sites, 2);
	  (System->chaininfo) = createNonSquareIntMatrix(length2, 4);

	  
  
    sprintf(fullname, "%s.siteinfo", filename);
    fileout = fopen(fullname, "r");
    
    for(j=0; j<tot_sites; j++)
      {
	fscanf(fileout, "%d	%d", &(System->siteinfo)[j][0], &(System->siteinfo)[j][1]);
      }
    fclose(fileout);

    sprintf(fullname, "%s.pos", filename);
    fileout = fopen(fullname, "r");
    
    for(j=0; j<tot_sites; j++)
      {
	fscanf(fileout, "%lf	%lf", &(System->pos)[j][0], &(System->pos)[j][1]);
      }
    fclose(fileout);
    
    
    sprintf(fullname, "%s.site_pots", filename);
    fileout = fopen(fullname, "r");
   
      for(j=0; j<tot_sites; j++)
      {
	fscanf(fileout, "%lf", &(System->site_pots)[j]);
      }
   
    fclose(fileout);

    
    
    sprintf(fullname, "%s.chaininfo", filename);
    fileout = fopen(fullname, "r");
    
   
      for(j=0; j<length2; j++)
      {
	fscanf(fileout, "%d	%d	%d	%d", &(System->chaininfo)[j][0], &(System->chaininfo)[j][1], &(System->chaininfo)[j][2], &(System->chaininfo)[j][3]);
      }
    fclose(fileout);
   
  
  sprintf(fullname, "%s.info", filename);
  fileout = fopen(fullname, "r");
    
       fscanf(fileout, "%d", (System->Nrem));
    
  fclose(fileout);
  


}