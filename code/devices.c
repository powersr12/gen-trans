

#include "devices.h"


//Routines for generating various devices (RectRedux form)


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