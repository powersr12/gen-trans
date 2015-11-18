

#include "devices.h"


//Routines for generating various devices (RectRedux form)

//creates the simplest RectRedux array for a finite number of ribbon chains
//note that this format should be the general layout for "generateDevice" routines
//with other disorder or antidot parameters placed into parameter structure p 
void simpleRibbonGeo (RectRedux *SiteArray, void *p, int struc_out, char *filename)
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



//very general - creates ribbon with sublattice dependent potential interface
void genSublatticeInterface(RectRedux *SiteArray, void *p, int struc_out, char *filename)
{  
	  subint_para *params = (subint_para *)p;
	  double a_conc1 = (params->a_conc1);
	  double b_conc1 = (params->b_conc1);
	  double a_pot1 = (params->a_pot1);
	  double b_pot1 = (params->b_pot1);
	  double a_conc2 = (params->a_conc2);
	  double b_conc2 = (params->b_conc2);
	  double a_pot2 = (params->a_pot2);
	  double b_pot2 = (params->b_pot2);
	  int xory = (params->xory);
	  double int_pos = (params->int_pos);
	  double int_width = (params->int_width);
	  int buffer_rows = (params->buffer_rows);
	  int seed = (params->seed);
	  
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
	      //both concentrations & potential vary linearly across the interface
	      double temprandnum;
	      double a_conc_temp, b_conc_temp, a_pot_temp, b_pot_temp;
	      double side1, side2, centre, relpos;
	      double delta_conca, delta_pota, delta_concb, delta_potb;
	      side1 = int_pos - (int_width/2);
	      side2 = int_pos + (int_width/2);
	      delta_conca = a_conc2 - a_conc1;
	      delta_concb = b_conc2 - b_conc1;
	      delta_pota = a_pot2 - a_pot1;
	      delta_potb = b_pot2 - b_pot1;
	      
	      for(l=buffer_rows*2*length; l<2*length*length2 - 2*buffer_rows*length; l++)
	      {
		site_pots[l] = 0.0;
		temprandnum = myRandNum(0.0, 1.0);
		
		
		
		if(xory==0)
		{
		  relpos = (site_coords[l][0] - side1)/int_width;
		}
		
		if(xory==1)
		{
		  relpos = (site_coords[l][1] - side1)/int_width;
		}
		
		
		  
		  if(relpos <= 0)
		  {
		    a_conc_temp = a_conc1;
		    a_pot_temp = a_pot1;
		    b_conc_temp = b_conc1;
		    b_pot_temp = b_pot1;
		  }
		  
		  if(relpos > 0 && relpos < 1)
		  {
		    a_conc_temp = a_conc1 + delta_conca*relpos;
		    b_conc_temp = b_conc1 + delta_concb*relpos;
		    a_pot_temp = a_pot1 + delta_pota*relpos;
		    b_pot_temp = b_pot1 + delta_potb*relpos;
		  }
		  
		    if(relpos >=1)
		  {
		    a_conc_temp = a_conc2;
		    a_pot_temp = a_pot2;
		    b_conc_temp = b_conc2;
		    b_pot_temp = b_pot2;
		  }
		  
		  
		  
		
		
		
		if(siteinfo[l][1] == 0)
		{
		  if(temprandnum < a_conc_temp)
		    site_pots[l] = a_pot_temp;
		}
		
		if(siteinfo[l][1] == 1)
		{
		  if(temprandnum < b_conc_temp)
		    site_pots[l] = b_pot_temp;
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
	  
	  
	  
			//separates out sublattices and positive and negative potentials by size
			if(struc_out == 1)
			{

			  for(l=0; l<2*length*length2; l++)
			  {
			    if(siteinfo[l][0] == 0  && siteinfo[l][1] == 0 && site_pots[l] > 0)
			    {  
			      fprintf(out, "%lf	%lf	%lf\n", site_coords[l][0], site_coords[l][1], site_pots[l] );
			    }
			  }
			  fprintf(out, "\n");
			  
			  for(l=0; l<2*length*length2; l++)
			  {
			    if(siteinfo[l][0] == 0  && siteinfo[l][1] == 0 && site_pots[l] < 0)
			    {  
			      fprintf(out, "%lf	%lf	%lf\n", site_coords[l][0], site_coords[l][1], site_pots[l] );
			    }
			  }
			  fprintf(out, "\n");
			  
			  for(l=0; l<2*length*length2; l++)
			  {
			    if(siteinfo[l][0] == 0  && siteinfo[l][1] == 1 && site_pots[l] > 0)
			    {  
			      fprintf(out, "%lf	%lf	%lf\n", site_coords[l][0], site_coords[l][1], site_pots[l] );
			    }
			  }
			  fprintf(out, "\n");
			  
			  for(l=0; l<2*length*length2; l++)
			  {
			    if(siteinfo[l][0] == 0  && siteinfo[l][1] == 1 && site_pots[l] < 0)
			    {  
			      fprintf(out, "%lf	%lf	%lf\n", site_coords[l][0], site_coords[l][1], site_pots[l] );
			    }
			  }
			  fprintf(out, "\n");
			
			  
			  if(struc_out != 0)
			  {
			    fclose(out);
			  }
			}
			    
			
			  
			  //Fill the array of data structures describing the system
			  (SiteArray->pos) = site_coords;
			  (SiteArray->site_pots) = site_pots;
			  (SiteArray->siteinfo) = siteinfo;

	
}






//very general - creates ribbon with sublattice dependent potential
void genSublatticeDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename)
{  
	  subl_para *params = (subl_para *)p;
	  double a_conc = (params->a_conc);
	  double b_conc = (params->b_conc);
	  double a_pot = (params->a_pot);
	  double b_pot = (params->b_pot);
	  int buffer_rows = (params->buffer_rows);
	  int seed = (params->seed);
	  
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
void genAntidotDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename)
{  
	  adot_para *params = (adot_para *)p;
	  
	  int buffer_rows = (params->buffer_rows);
	  int AD_length = (params->AD_length);
	  int AD_length2 = (params->AD_length2);
	  double AD_rad = (params->AD_rad);
	  double AD_rad2 = (params->AD_rad2);

	  int lat_width = (params->lat_width);
	  int lat_length = (params->lat_length);
	  char *latgeo = (params->latgeo);
	  char *dotgeo = (params->dotgeo);
	  int isperiodic = (params->isperiodic);
	  
	  double radfluc = (params->radfluc);
	  double xyfluc = (params->xyfluc);


	  int seed = (params->seed);
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
	  int holes_per_cell, tot_holes, rad_per_hole;
	  double renorm_fac=1.0; //converts side length to "radius" in certain cases	
	  
	  if(strcmp("trig", latgeo) == 0)
	  {
	    holes_per_cell=2;
	  }
	  
	  if(strcmp("rotrig", latgeo) == 0)
	  {
	    holes_per_cell=2;
	  }
	  
	  if(strcmp("rect", latgeo) == 0)
	  {
	    holes_per_cell=1;
	  }
	  
	  if(strcmp("circ", dotgeo) == 0)
	  {
	    rad_per_hole=1;
	  }
	  if(strcmp("hexAC", dotgeo) == 0)
	  {
	    rad_per_hole=6;
	    renorm_fac = sqrt(3) / 2;
	  }
	  if(strcmp("hexZZ", dotgeo) == 0)
	  {
	    rad_per_hole=6;
	    renorm_fac = sqrt(3)/2;
	  }
	  if(strcmp("rect", dotgeo) == 0)
	  {
	    rad_per_hole=4;
	    renorm_fac = 0.5;
	  }
	  
	  
	  
	  tot_holes = holes_per_cell*lat_width*lat_length;
	  double holes[tot_holes][2+rad_per_hole];
	  double unitholes[holes_per_cell][2];
	  double xshift, yshift;
	  
	  double cellyshift;
	  
	  if(geo==0)
	    cellyshift = length * sqrt(3) / 2;
	  
	  else if (geo==1)
	    cellyshift = length * 0.5;
	  
	  //unit cell hole positions (centres) and shift vectors
	  
	  //triangular lattice
	  if(strcmp("trig", latgeo) == 0)
	  {
	  
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
		
		
		//armchair ribbon 
		if(geo==1)
		{
		  xstart = buffer_rows*sqrt(3);
		  ystart = 0.0;
		  xshift = sqrt(3)*AD_length;
		  yshift = 3.0*AD_length;
		  
	
		  unitholes[0][1] =3.0*AD_length - 0.5 - ((int) (AD_length*3.0/4) -0.5* ( ((int) (AD_length/2)) % 2 ));
		  unitholes[0][0] = (sqrt(3)/2.0) * (int) (AD_length /2.0)  - 1 /(2*sqrt(3)) ;
		  unitholes[1][1] = unitholes[0][1] - 1.5*AD_length;
		  unitholes[1][0] = unitholes[0][0] + sqrt(3)*AD_length/2.0;

		  
		  
		}
	  }
	  
	  
	  //rotated triangular lattice
	  if(strcmp("rotrig", latgeo) == 0)
	  {
	  
		//zigzag ribbon
		if(geo==0)
		{
		  xstart = buffer_rows*1.0;
		  ystart = 0.0;
		  xshift = (AD_length + 1)*1.0;
		  yshift = (AD_length+1)*sqrt(3);
		  
		  unitholes[0][0] = (2*AD_length+1)*0.5 - ((int)(AD_length/4))*1.0 - 0.5; 
		  unitholes[0][1] = (((int)(AD_length/4)) +0.5)*sqrt(3) + sqrt(3)/2;
		  unitholes[1][0] = unitholes[0][0] - (AD_length+1)*0.5  ;
		  unitholes[1][1] = unitholes[0][1] + (AD_length+1)*sqrt(3)*0.5 ;
		  
		}
		
		
		//armchair ribbon 
		if(geo==1)
		{
		  xstart = buffer_rows*sqrt(3);
		  ystart = 0.0;
		  xshift = (AD_length+1)*sqrt(3);
		  yshift = (AD_length + 1)*1.0;
		  
	
		  unitholes[0][1] = ((int)(AD_length/4))*1.0 + 0.5; 
		  unitholes[0][0] = (((int)(AD_length/4)) +0.5)*sqrt(3) + 1/(sqrt(3));;
		  unitholes[1][1] = unitholes[0][1] + (AD_length+1)*0.5  ;
		  unitholes[1][0] = unitholes[0][0] + (AD_length+1)*sqrt(3)*0.5 ;

		  
		  
		}
	  }
	  
	  //rectangular lattice
	      //ADlength is armchair direction unit vector (in units of sqrt(3)a )
	      //ADlength2 is zigzag direction unit vector  (in units of a)
	      
	  if(strcmp("rect", latgeo) == 0)
	  {
	  
		//zigzag ribbon
		if(geo==0)
		{
		  xstart = buffer_rows*1.0;
		  ystart = 0.0;
		  xshift = AD_length2;
		  yshift = sqrt(3)*AD_length;
		  
		  unitholes[0][0]  = (int) (AD_length2/2) - 0.5* (  (AD_length % 2) ); 
		  unitholes[0][1] = (sqrt(3)/2.0) * (int) (AD_length)   ;
		  
		}
		
		
		//armchair ribbon 
		if(geo==1)
		{
		  xstart = buffer_rows*sqrt(3);
		  ystart = 0.0;
		  xshift = sqrt(3)*AD_length;
		  yshift = AD_length2;
		  
	
		  unitholes[0][1]  = AD_length2 -0.5 - (int) (AD_length2/2) +0.5* (  (AD_length % 2) ); 
		  unitholes[0][0] = (sqrt(3)/2.0) * (int) (AD_length)  - 1 /(2*sqrt(3)) ;
		  
		 	  

		  
		  
		}
	  }
	  
	  
	  
	  for(i=0; i<lat_length; i++)
	  {
	    for(j=0; j<lat_width; j++)
	    {
	      for(k=0; k <holes_per_cell; k++)
	      {
		holes[holes_per_cell*i*lat_width + holes_per_cell*j + k][0] = xstart + i*xshift + unitholes[k][0] + myRandNum(-xyfluc, xyfluc);
		holes[holes_per_cell*i*lat_width + holes_per_cell*j + k][1] = ystart + j*yshift + unitholes[k][1] + myRandNum(-xyfluc, xyfluc);
		
		if(strcmp("rect", dotgeo) != 0)
		{
		  for(l=0; l<rad_per_hole; l++)
		  {
		    holes[holes_per_cell*i*lat_width + holes_per_cell*j + k][2+l] = AD_rad * renorm_fac + myRandNum(-radfluc, radfluc);
		  }
		}
		
		if(strcmp("rect", dotgeo) == 0)
		{
		  holes[holes_per_cell*i*lat_width + holes_per_cell*j + k][2] = AD_rad * renorm_fac + myRandNum(-radfluc, radfluc);
		  holes[holes_per_cell*i*lat_width + holes_per_cell*j + k][3] = AD_rad2 * renorm_fac + myRandNum(-radfluc, radfluc);
		  holes[holes_per_cell*i*lat_width + holes_per_cell*j + k][4] = AD_rad * renorm_fac + myRandNum(-radfluc, radfluc);
		  holes[holes_per_cell*i*lat_width + holes_per_cell*j + k][5] = AD_rad2 * renorm_fac + myRandNum(-radfluc, radfluc);
		}
		
	      }
	    }
	  }
	  
	  
	  int **sites = createNonSquareIntMatrix(tot_sites, 2); //removed, num neighbours after first sweep
	  //double **vertices = createNonSquareDoubleMatrix(rad_per_hole, 2);
	  double *polyx = createDoubleArray(rad_per_hole);
	  double *polyy = createDoubleArray(rad_per_hole);
	  
	   
	  //atom removal!
	      for(i=2*length; i< tot_sites - 2*length ; i++)
	      {
		for(j=0; j< tot_holes; j++)
		{
		  
		    if(strcmp("circ", dotgeo) == 0)
		    {
			if(sqrt( pow( site_coords[i][0] - holes[j][0], 2) + pow( site_coords[i][1] - holes[j][1], 2)) < holes[j][2])
			{
			  sites[i][0] = 1;
			}
		    }
		    
		    if(strcmp("rect", dotgeo) == 0)
		    {
		      polyx[0] = holes[j][0] + holes[j][2];
		      polyy[0] = holes[j][1] + holes[j][3];
		      polyx[1] = holes[j][0] - holes[j][4];
		      polyy[1] = holes[j][1] + holes[j][3];
		      polyx[2] = holes[j][0] - holes[j][4];
		      polyy[2] = holes[j][1] - holes[j][5];
		      polyx[3] = holes[j][0] + holes[j][2];
		      polyy[3] = holes[j][1] - holes[j][5];
		      
		      if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
		      {
			sites[i][0] = 1;
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
			  
			
		      
		    }
		    
		    if(strcmp("hexZZ", dotgeo) == 0)
		    {
		      if(geo == 0)
		      {
			polyx[0] = holes[j][0] + (2*holes[j][2] - holes[j][3])/sqrt(3);
			polyy[0] = holes[j][1] + holes[j][3];
			polyx[1] = holes[j][0] + (holes[j][3] - 2*holes[j][4])/sqrt(3);
			polyy[1] = holes[j][1] + holes[j][3];
			polyx[2] = holes[j][0] - (holes[j][4] + holes[j][5])/sqrt(3);
			polyy[2] = holes[j][1] + (holes[j][4] - holes[j][5]);
			polyx[3] = holes[j][0] + (holes[j][6] - 2*holes[j][5])/sqrt(3);
			polyy[3] = holes[j][1] - holes[j][6];
			polyx[4] = holes[j][0] + (2*holes[j][7] - holes[j][6])/sqrt(3);
			polyy[4] = holes[j][1] - holes[j][6];
			polyx[5] = holes[j][0] + (holes[j][2] + holes[j][7])/sqrt(3);
			polyy[5] = holes[j][1] + (holes[j][2] - holes[j][7]);
		      }
		      
		      if(geo == 1)
		      {
			polyx[0] = holes[j][0] + holes[j][3];
			polyy[0] = holes[j][1] + (2*holes[j][2] - holes[j][3])/sqrt(3);
			polyx[1] = holes[j][0] + holes[j][3];
			polyy[1] = holes[j][1] + (holes[j][3] - 2*holes[j][4])/sqrt(3);
			polyx[2] = holes[j][0] + (holes[j][4] - holes[j][5]);
			polyy[2] = holes[j][1] - (holes[j][4] + holes[j][5])/sqrt(3);
			polyx[3] = holes[j][0] - holes[j][6];
			polyy[3] = holes[j][1] + (holes[j][6] - 2*holes[j][5])/sqrt(3);
			polyx[4] = holes[j][0] - holes[j][6];
			polyy[4] = holes[j][1] + (2*holes[j][7] - holes[j][6])/sqrt(3);
			polyx[5] = holes[j][0] + (holes[j][2] - holes[j][7]);
			polyy[5] = holes[j][1] + (holes[j][2] + holes[j][7])/sqrt(3);
		      }
		      
		      
		      if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
		      {
			sites[i][0] = 1;
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
		      
		    }
		    
		    
		    if(strcmp("hexAC", dotgeo) == 0)
		    {
		      if(geo == 0)
		      {
			polyx[0] = holes[j][0] + holes[j][3];
			polyy[0] = holes[j][1] + (2*holes[j][2] - holes[j][3])/sqrt(3);
			polyx[1] = holes[j][0] + holes[j][3];
			polyy[1] = holes[j][1] + (holes[j][3] - 2*holes[j][4])/sqrt(3);
			polyx[2] = holes[j][0] + (holes[j][4] - holes[j][5]);
			polyy[2] = holes[j][1] - (holes[j][4] + holes[j][5])/sqrt(3);
			polyx[3] = holes[j][0] - holes[j][6];
			polyy[3] = holes[j][1] + (holes[j][6] - 2*holes[j][5])/sqrt(3);
			polyx[4] = holes[j][0] - holes[j][6];
			polyy[4] = holes[j][1] + (2*holes[j][7] - holes[j][6])/sqrt(3);
			polyx[5] = holes[j][0] + (holes[j][2] - holes[j][7]);
			polyy[5] = holes[j][1] + (holes[j][2] + holes[j][7])/sqrt(3);
			
		      }
		      
		      if(geo == 1)
		      {
			polyx[0] = holes[j][0] + (2*holes[j][2] - holes[j][3])/sqrt(3);
			polyy[0] = holes[j][1] + holes[j][3];
			polyx[1] = holes[j][0] + (holes[j][3] - 2*holes[j][4])/sqrt(3);
			polyy[1] = holes[j][1] + holes[j][3];
			polyx[2] = holes[j][0] - (holes[j][4] + holes[j][5])/sqrt(3);
			polyy[2] = holes[j][1] + (holes[j][4] - holes[j][5]);
			polyx[3] = holes[j][0] + (holes[j][6] - 2*holes[j][5])/sqrt(3);
			polyy[3] = holes[j][1] - holes[j][6];
			polyx[4] = holes[j][0] + (2*holes[j][7] - holes[j][6])/sqrt(3);
			polyy[4] = holes[j][1] - holes[j][6];
			polyx[5] = holes[j][0] + (holes[j][2] + holes[j][7])/sqrt(3);
			polyy[5] = holes[j][1] + (holes[j][2] - holes[j][7]);
		      }
		      
		      
		      if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
		      {
			sites[i][0] = 1;
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
		      
		    }
		    
		}
	      }
	      
	      free(polyx); free(polyy);
	  
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
		  if(isperiodic == 1)
		  {
		    for(j=0 ; j< tot_sites ;j++)
		    {
		      if(sites[j][0] == 0)
		      {
			if(i!=j)
			{
			  smalldist = sqrt( pow(site_coords[i][0] - site_coords[j][0], 2) +  pow(site_coords[i][1] - site_coords[j][1] -cellyshift, 2));
			
			
			  if (smalldist < 0.6)
			  {
			    tempint++;
			  }
			}
		      }
		    }
		  
		   for(j=0 ; j< tot_sites ;j++)
		    {
		      if(sites[j][0] == 0)
		      {
			if(i!=j)
			{
			  smalldist = sqrt( pow(site_coords[i][0] - site_coords[j][0], 2) +  pow(site_coords[i][1] - site_coords[j][1] +cellyshift, 2));
			
			
			  if (smalldist < 0.6)
			  {
			    tempint++;
			  }
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
  char fullname[200];
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
       fprintf(fileout, "%d\n", *(System->Ntot));

  fclose(fileout);

}


void importRectConf(RectRedux *System, int length, int length2, char *filename)
{
  //printf("%s	%d\n", filename, numcells);
  FILE *fileout;
  char fullname[200];
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
       fscanf(fileout, "%d", (System->Ntot));

  fclose(fileout);
  


}