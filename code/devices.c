
#include "devices.h"
#include "connect.h"
#include <stdio.h>
#include "../libs/matrices.h"
#include "useful.h"
#include "time.h"
#include <string.h>
#include <stdio.h>



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
// 		      (SiteArray->chaininfo) = createNonSquareIntMatrix(length2, 4);
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
// 			  (SiteArray->chaininfo)[l][0] = tempint;
			}
			*Nrem = tempint2;
			*Ntot = *Nrem;
			/*
			(SiteArray->chaininfo)[0][1] = 0;
			for(l=1; l<length2; l++)
			{
			  (SiteArray->chaininfo)[l][1] = (SiteArray->chaininfo)[l-1][1] + (SiteArray->chaininfo)[l-1][0];
			}*/


			//are first and last atoms in each chain present?
// 			for(l=1; l<length2; l++)
// 			{
// 			  if(siteinfo[l*2*length][0] == 1)
// 			    (SiteArray->chaininfo)[l][2] = 1;
// 			  
// 			  if(siteinfo[(l+1)*2*length -1][0] == 1)
// 			    (SiteArray->chaininfo)[l][3] = 1;
// 			}
// 	  
	  
	  
		      
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
                          (SiteArray->pert_pos) = NULL;
  
  
}



//simple 558GB in agnr device
void simple558GB (RectRedux *SiteArray, void *p, int struc_out, char *filename)
{
	  int length = SiteArray->length;
	  int length2 = SiteArray->length2;
	  int geo = SiteArray->geo;
	    
	  
	simple558_para *params = (simple558_para *)p;

	  double smalldist;
	  FILE *out;
	  
	  if(struc_out != 0)
	  {
	    out = fopen(filename, "w");
	    
	  }
	  int seed = (params->seed);

	  srand(time(NULL) + seed);
	  
	  int cellswide = params->cellswide;	//how many blocks of GB wide the AGNR is 
	  int GBpos = params->GBpos; 
	  int anddis= params->anddis;
	  double andD= params->andD;
	  double andW= params->andW;
	  int vacdis= params->vacdis;
	  double vacW= params->vacW;
	  double vacP= params->vacP;
	  double ycut = params->ycut;
	  
	  double ymax=(length-1)*0.5, ymin=0.0; 

	  
	  //atomic coordinates and the atoms that are in and out, similar to previous cases
	  int tot_sites = 2*length*length2 + 2*cellswide;
	  double **site_coords = createNonSquareDoubleMatrix(tot_sites, 3);
	  double *site_pots = createDoubleArray(tot_sites);
	  int **siteinfo = createNonSquareIntMatrix(tot_sites, 2);
	  double xstart, ystart;
	  int isclean, l, m;
	  int *Nrem = (SiteArray->Nrem);
	  int *Ntot = (SiteArray->Ntot);
	  double temprandnum;
	  
	  *Ntot = tot_sites;

	  if(geo==0)
	  {
		exit(1);
	  }
	  
	  if(geo==1)
	  {
	      for(l=0; l<length2; l++)
	      {
		xstart = l*sqrt(3);
		
		if(l>(GBpos -1))
		{
			xstart+= 1/sqrt(3);
		}
		
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
	      
	      for(m=0; m<cellswide; m++)
	      {
		site_coords[2*length*length2 + 2*m][0] = GBpos*sqrt(3) ;  
		site_coords[2*length*length2 + 2*m + 1 ][0] = GBpos*sqrt(3) ;
		
		site_coords[2*length*length2 + 2*m][1] = m*2.0 + 0.2;
		site_coords[2*length*length2 + 2*m + 1][1] = m*2.0 + 0.8;
		
		siteinfo[2*length*length2 + 2*m][1] = 2;
		siteinfo[2*length*length2 + 2*m + 1][1] = 2;
		
	      }
		
		//Vacancy disorder near the GB
		if(vacdis == 1)
		{
			for(l=0; l<tot_sites; l++)
			{
				if(fabs(site_coords[l][0] - GBpos*sqrt(3)) < vacW )
				{
					if(fabs(site_coords[l][1] - ymin) > ycut && fabs(site_coords[l][1] - ymax) > ycut)
					{
						temprandnum = myRandNum(0.0, 1.0);
						if(temprandnum < vacP)
						{
							siteinfo[l][0] = 1;
						}
					}
				}
			}
		}
		
		
		
		
	      
		
		
		//Anderson disorder near the GB
		if(anddis == 1)
		{
			for(l=0; l<tot_sites; l++)
			{
				if(fabs(site_coords[l][0] - GBpos*sqrt(3)) < andD && siteinfo[l][0] == 0)
				{
					if(fabs(site_coords[l][1] - ymin) > ycut && fabs(site_coords[l][1] - ymax) > ycut)
					{
						site_pots[l] = myRandNum(- andW, andW);
					}
				}
			}
		}
	  }
	      
	          
	      	      
	      //chaininfo needed for conductance calcs (if atoms are missing)
		      int tempint=0, tempint2=0;
	 
			for(l=0; l<tot_sites; l++)
			{
			  if(siteinfo[l][0] ==0)
				  tempint2++;
			}
			*Nrem = tempint2;
			

	  
		      
			if(struc_out == 1)
			{
			  for(l=0; l<tot_sites; l++)
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
                          (SiteArray->pert_pos) = NULL;
  
  
}



//generate lead geometries and rect_reduxes to match a device
//mode=0 - basic left-right setup, same size 
//mode=1 - BLG variant
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
    
    if(numleads != 2)
      exit(1);
    
    //shift leads accordingly to desired positions
    for(i=0; i<*(Leads[0]->Ntot); i++)
    {
      (Leads[0]->pos)[i][0] += (params->shift_vecs)[0][0];      
    }
    
    for(i=0; i<*(Leads[1]->Ntot); i++)
    {
      //(Leads[1]->pos)[i][0] += (SiteArray->length2)*(params->shift_vecs)[1][0];      
      
       (Leads[1]->pos)[i][0] += ((SiteArray->pos)[2*(SiteArray->length)*((SiteArray->length2) -1)][0] - (SiteArray->pos)[0][0]   + (params->shift_vecs)[1][0]); 
      //printf("## %lf\n", (Leads[1]->pos)[i][0]);
    }
    
  }
  
  
  if(mode==1)
  {
	//if this mode is called, bilayer_para are stored in "additional_params" in lead_para
	  
	for(i=0; i < numleads; i++)
	{
		Leads[i] = (RectRedux *)malloc(sizeof(RectRedux));
		(Leads[i]->geo) = (SiteArray->geo);
		(Leads[i]->length) = (SiteArray->length);
		(Leads[i]->length2) = 1;
		(Leads[i]->Nrem) = (int *)malloc(sizeof(int));
		(Leads[i]->Ntot) = (int *)malloc(sizeof(int));
		simpleBilayerGeo (Leads[i], params->additional_params, 0, NULL);
	}
	
	
	(params->shift_vecs) = createNonSquareDoubleMatrix(numleads, 3);
	(params->shift_vecs)[0][0] = -((SiteArray->pos)[2*(SiteArray->length)][0] - (SiteArray->pos)[0][0]);
	(params->shift_vecs)[1][0] = (SiteArray->pos)[2*(SiteArray->length)][0] - (SiteArray->pos)[0][0];
	(params->shift_vecs)[0][1] = 0;
	(params->shift_vecs)[1][1] = 0;
	
	if(numleads != 2)
		exit(1);
	
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
  
}
 
void genSingleRibbonLead (RectRedux *SiteArray, RectRedux *Lead, int lead_num, void *params)
{
	int i, j;
	rib_lead_para* ribpara = (rib_lead_para*)params;
	
	int devgeo = SiteArray->geo;
	int length = (SiteArray->length);
	int length2 = (SiteArray->length2);
	double ybot, ytop, y_cell_diff;
	double xleft, xright, x_cell_diff;
	
	if(devgeo == 0)
	{
		ybot = 1/(2*sqrt(3)) - sqrt(3);
		ytop = 1/(2*sqrt(3)) + length * sqrt(3)/2;
		xleft = 0.5;
		xright = length2*1.0 - 0.5;
		x_cell_diff = 1.0;
		y_cell_diff = sqrt(3);
	}
	
	if(devgeo == 1)
	{
		ybot = -1.0;
		ytop = 0.5*length;
		xleft = -1/(2*sqrt(3)) + sqrt(3)/2;
		xright = -1/(2*sqrt(3))+ sqrt(3)/2 + (length2-1)*sqrt(3);
		x_cell_diff = sqrt(3);
		y_cell_diff = 1.0;
	}
    
	
	(Lead->geo) = (ribpara->geo);
	(Lead->length) = (ribpara->width);
	(Lead->length2) = 1;
	
	simpleRibbonGeo (Lead, NULL, 0, NULL);

	//periodicity vectors for RGF leads
	(ribpara->shift_vec) = createDoubleArray(3);
	
	
	//sublattice pots
	for(i=0; i<*(Lead->Ntot); i++)
	{
		if( (Lead->siteinfo)[i][1] == 0)
		{	
			(Lead->site_pots)[i]=ribpara->homo_A_pot;
		}
		if( (Lead->siteinfo)[i][1] == 1)
		{	
			(Lead->site_pots)[i]=ribpara->homo_B_pot;
		}
	}
	
	if(ribpara->def_pos == 0)
	{
		(ribpara->shift_vec)[0] = -((SiteArray->pos)[2*(SiteArray->length)][0] - (SiteArray->pos)[0][0]);
		(ribpara->shift_vec)[1] = 0.0;
		for(i=0; i<*(Lead->Ntot); i++)
		{
			(Lead->pos)[i][0] += (ribpara->shift_vec)[0];
			(Lead->pos)[i][1] += (ribpara->start_coord)*y_cell_diff/2;
		}
	}
    
	if(ribpara->def_pos == 1)
	{
		(ribpara->shift_vec)[0] = (SiteArray->pos)[2*(SiteArray->length)][0] - (SiteArray->pos)[0][0];
		(ribpara->shift_vec)[1] = 0.0;
		for(i=0; i<*(Lead->Ntot); i++)
		{
			//(Lead->pos)[i][0] += (SiteArray->length2)*(ribpara->shift_vec)[0];
			(Lead->pos)[i][0] += ((SiteArray->pos)[2*(SiteArray->length)*((SiteArray->length2) -1)][0] - (SiteArray->pos)[0][0]   + (ribpara->shift_vec)[0]);
			
			(Lead->pos)[i][1] += (ribpara->start_coord)*y_cell_diff/2;
			//printf("## %lf, %lf\n", (Lead->pos)[i][0], (Lead->pos)[i][1]);
		}
	}
	if(ribpara->def_pos == 2)
	{
		(ribpara->shift_vec)[0] = 0.0;
		(ribpara->shift_vec)[1] = y_cell_diff;
		
		swapxy((Lead->pos), *(Lead->Ntot));
		
		for(i=0; i<*(Lead->Ntot); i++)
		{
			(Lead->pos)[i][0] += (xleft + (ribpara->start_coord)*x_cell_diff);
			(Lead->pos)[i][1] += (ytop);
		}
	}
	if(ribpara->def_pos == 3)
	{
		(ribpara->shift_vec)[0] = 0.0;
		(ribpara->shift_vec)[1] = -y_cell_diff;
		
		swapxy((Lead->pos), *(Lead->Ntot));
		
		for(i=0; i<*(Lead->Ntot); i++)
		{
			(Lead->pos)[i][0] += (xleft + (ribpara->start_coord)*x_cell_diff);
			(Lead->pos)[i][1] += (ybot);
		}
		
	}
} 


void genSingleBLGLead (RectRedux *SiteArray, RectRedux *Lead, int lead_num, void *params)
{
	int i, j;
	blg_lead_para* ribpara = (blg_lead_para*)params;
	
	int devgeo = SiteArray->geo;
	int length = (SiteArray->length);
	int length2 = (SiteArray->length2);
	double ybot, ytop, y_cell_diff;
	double xleft, xright, x_cell_diff;
	
	if(devgeo == 0)
	{
		ybot = 1/(2*sqrt(3)) - sqrt(3);
		ytop = 1/(2*sqrt(3)) + length * sqrt(3)/2;
		xleft = 0.5;
		xright = length2*1.0 - 0.5;
		x_cell_diff = 1.0;
		y_cell_diff = sqrt(3);
	}
	
	if(devgeo == 1)
	{
		ybot = -1.0;
		ytop = 0.5*length;
		xleft = -1/(2*sqrt(3)) + sqrt(3)/2;
		xright = -1/(2*sqrt(3))+ sqrt(3)/2 + (length2-1)*sqrt(3);
		x_cell_diff = sqrt(3);
		y_cell_diff = 1.0;
	}
    
	
	(Lead->geo) = (ribpara->geo);
	(Lead->length) = (ribpara->width);
	(Lead->length2) = 1;
	
	
	simpleBilayerGeo (Lead, ribpara->bilayer_para, 0, NULL);

	//periodicity vectors for RGF leads
	(ribpara->shift_vec) = createDoubleArray(3);
    
	if(ribpara->def_pos == 0)
	{
		(ribpara->shift_vec)[0] = -((SiteArray->pos)[2*(SiteArray->length)][0] - (SiteArray->pos)[0][0]);
		(ribpara->shift_vec)[1] = 0.0;
		for(i=0; i<*(Lead->Ntot); i++)
		{
			(Lead->pos)[i][0] += (ribpara->shift_vec)[0];
			(Lead->pos)[i][1] += (ribpara->start_coord)*y_cell_diff/2;
		}
	}
    
	if(ribpara->def_pos == 1)
	{
		(ribpara->shift_vec)[0] = (SiteArray->pos)[2*(SiteArray->length)][0] - (SiteArray->pos)[0][0];
		(ribpara->shift_vec)[1] = 0.0;
		for(i=0; i<*(Lead->Ntot); i++)
		{
			(Lead->pos)[i][0] += (SiteArray->length2)*(ribpara->shift_vec)[0];
			(Lead->pos)[i][1] += (ribpara->start_coord)*y_cell_diff/2;
		}
	}
	if(ribpara->def_pos == 2)
	{
		(ribpara->shift_vec)[0] = 0.0;
		(ribpara->shift_vec)[1] = y_cell_diff;
		
		swapxy((Lead->pos), *(Lead->Ntot));
		
		for(i=0; i<*(Lead->Ntot); i++)
		{
			(Lead->pos)[i][0] += (xleft + (ribpara->start_coord)*x_cell_diff);
			(Lead->pos)[i][1] += (ytop);
		}
	}
	if(ribpara->def_pos == 3)
	{
		(ribpara->shift_vec)[0] = 0.0;
		(ribpara->shift_vec)[1] = -y_cell_diff;
		
		swapxy((Lead->pos), *(Lead->Ntot));
		
		for(i=0; i<*(Lead->Ntot); i++)
		{
			(Lead->pos)[i][0] += (xleft + (ribpara->start_coord)*x_cell_diff);
			(Lead->pos)[i][1] += (ybot);
		}
		
	}
} 

//def_pos = 0,..3 - left, right, top, bottom - widths and coords in ribbonlike notation
//def_pos = 4 - full width of ribbon, coords and length in ribbon like notation
//def_pos = 5 - start_coord, start_coord2, witdth, width2 in absolute units
void genSingleMetalLead (RectRedux *SiteArray, RectRedux *Lead, int lead_num, void *params)
{
	int i, j;
	metal_lead_para* metpara = (metal_lead_para*)params;
	
	int devgeo = SiteArray->geo;
	int length = (SiteArray->length);
	int length2 = (SiteArray->length2);
	double ybot, ytop, y_cell_diff;
	double xleft, xright, x_cell_diff;
	
	if(devgeo == 0)
	{
		ybot = 0.0;
		ytop = 1/(2*sqrt(3)) + length * sqrt(3)/2;
		xleft = 0.5;
		xright = length2*1.0 - 0.5;
		x_cell_diff = 1.0;
		y_cell_diff = sqrt(3);
	}
	
	if(devgeo == 1)
	{
		ybot = 0.0;
		ytop = 0.5*length;
		xleft = -1/(2*sqrt(3)) + sqrt(3)/2;
		//xright = -1/(2*sqrt(3))+ sqrt(3)/2 + (length2-1)*sqrt(3);
		xright =  1/(2*sqrt(3)) + sqrt(3)/2 + (length2-1)*sqrt(3);

		x_cell_diff = sqrt(3);
		y_cell_diff = 1.0;
	}
    

	//Lead->pos is used here to give the max x/y of the lead
	(Lead->pos) = createNonSquareDoubleMatrix(2, 3);
    
	//lead sizes etc are such that they correspond to that of a ribbon with the same width param
	if(metpara->def_pos == 0)
	{
		(Lead->pos)[0][0] =xleft - (x_cell_diff/2);
		(Lead->pos)[0][1] = xleft - (x_cell_diff/2) + (metpara->width)*x_cell_diff -0.01;
		
		(Lead->pos)[1][0] = (metpara->start_coord2)*y_cell_diff/2;
		(Lead->pos)[1][1] = (metpara->start_coord2 + metpara->width2)*y_cell_diff/2;
	}
		
    
	if(metpara->def_pos == 1)
	{
		(Lead->pos)[0][0] =xright - (metpara->width)*x_cell_diff +0.01; // + (x_cell_diff/2);
		(Lead->pos)[0][1] =xright; // + (x_cell_diff/2) ;
		
		(Lead->pos)[1][0] = (metpara->start_coord2)*y_cell_diff/2;
		(Lead->pos)[1][1] = (metpara->start_coord2 + metpara->width2)*y_cell_diff/2;
			
	}
	
	if(metpara->def_pos == 2)
	{
		(Lead->pos)[0][0] = xleft + (metpara->start_coord)*x_cell_diff ;
		(Lead->pos)[0][1] = (Lead->pos)[0][0] + (metpara->width)*x_cell_diff/2 -x_cell_diff/2 ;
		
		(Lead->pos)[1][0] = ytop - (metpara->width2)*y_cell_diff/2 -0.001;
		(Lead->pos)[1][1] = ytop ;
		
	}
	
	if(metpara->def_pos == 3)
	{
		(Lead->pos)[0][0] = xleft + (metpara->start_coord)*x_cell_diff;
		(Lead->pos)[0][1] = (Lead->pos)[0][0] + (metpara->width)*x_cell_diff/2 -x_cell_diff/2;
		
		(Lead->pos)[1][0] = ybot -0.001;
		(Lead->pos)[1][1] = ybot  -0.001 + (metpara->width2)*y_cell_diff/2;
		
	}
	
	if(metpara->def_pos == 4)
	{
		(Lead->pos)[0][0] = xleft + (metpara->start_coord)*x_cell_diff;
		(Lead->pos)[0][1] = (Lead->pos)[0][0] + (metpara->width)*x_cell_diff/2 -x_cell_diff/2;
		
		(Lead->pos)[1][0] = ybot;
		(Lead->pos)[1][1] = ytop;
		
	}

	if(metpara->def_pos == 5)
	{
		(Lead->pos)[0][0] = (metpara->start_coord);
		(Lead->pos)[0][1] = (metpara->start_coord) + (metpara->width);
		
		(Lead->pos)[1][0] = (metpara->start_coord2);
		(Lead->pos)[1][1] = (metpara->start_coord2) + (metpara->width2);
		
	}
	

	
} 



//calls "generate" functions for each lead to turn the simple positions etc 
//in the main code into structures that can be used for calculating the startng 
//cells and later the self energies
//maybe just generate starting cells here too?)
void genCustomLeads (RectRedux *SiteArray, RectRedux **Leads, int numleads, lead_para *params)
{
  int i, j;
  leadgenfn *leadgn;
  void *indiv_params;
  
  for(i=0; i< numleads; i++)
  {
	Leads[i] = (RectRedux *)malloc(sizeof(RectRedux));
	(Leads[i]->Nrem) = (int *)malloc(sizeof(int));
	(Leads[i]->Ntot) = (int *)malloc(sizeof(int));
	*(Leads[i]->Nrem)=0;
	*(Leads[i]->Ntot)=0;
	
  }	
	
  for(i=0; i< numleads; i++)
  {	
	leadgn = ((params->multiple)[i])->indiv_gen_fn;
	indiv_params = ((params->multiple)[i])->indiv_lead_para;
	
	(leadgn) (SiteArray, Leads[i], i, indiv_params);
	
  }
  
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
                          (SiteArray->pert_pos) = NULL;

	
}



//very general - creates ribbon with two sublattice dependent potential interfaces
void genSublatticeTwoInts(RectRedux *SiteArray, void *p, int struc_out, char *filename)
{  
	  sub2int_para *params = (sub2int_para *)p;
	  double a_conc1 = (params->a_conc1);
	  double b_conc1 = (params->b_conc1);
	  double a_pot1 = (params->a_pot1);
	  double b_pot1 = (params->b_pot1);
	  double a_conc2 = (params->a_conc2);
	  double b_conc2 = (params->b_conc2);
	  double a_pot2 = (params->a_pot2);
	  double b_pot2 = (params->b_pot2);
	   double a_conc3 = (params->a_conc3);
	  double b_conc3 = (params->b_conc3);
	  double a_pot3 = (params->a_pot3);
	  double b_pot3 = (params->b_pot3);
	  int xory = (params->xory);
	  double int_pos1 = (params->int_pos1);
	  double int_width1 = (params->int_width1);
	  double int_pos2 = (params->int_pos2);
	  double int_width2 = (params->int_width2);
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
	      double side1a, side2a, side1b, side2b, centre, relposa, relposb;
	      double delta_conca, delta_pota, delta_concb, delta_potb;
	      double delta_conca2, delta_pota2, delta_concb2, delta_potb2;
	      side1a = int_pos1 - (int_width1/2);
	      side2a = int_pos1 + (int_width1/2);
	      side1b = int_pos2 - (int_width2/2);
	      side2b = int_pos2 + (int_width2/2);
	      delta_conca = a_conc2 - a_conc1;
	      delta_concb = b_conc2 - b_conc1;
	      delta_pota = a_pot2 - a_pot1;
	      delta_potb = b_pot2 - b_pot1;
	      delta_conca2 = a_conc3 - a_conc2;
	      delta_concb2 = b_conc3 - b_conc2;
	      delta_pota2 = a_pot3 - a_pot2;
	      delta_potb2 = b_pot3 - b_pot2;
	      
	      //printf("#%lf	%lf	%lf	%lf\n", side1a, side2a, side1b, side2b );
	      
	      for(l=buffer_rows*2*length; l<2*length*length2 - 2*buffer_rows*length; l++)
	      {
		site_pots[l] = 0.0;
		temprandnum = myRandNum(0.0, 1.0);
		
		
		
		if(xory==0)
		{
		  relposa = (site_coords[l][0] - side1a)/int_width1;
		  relposb = (site_coords[l][0] - side1b)/int_width2;
		}
		
		if(xory==1)
		{
		  relposa = (site_coords[l][1] - side1a)/int_width1;
		  relposb = (site_coords[l][1] - side1b)/int_width2;
		}
		
		if(relposa < 0)
		  relposa =0;
		
		if(relposa > 1)
		  relposa =1;
		
		if(relposb < 0)
		  relposb =0;
		
		if(relposb > 1)
		  relposb =1;
		
		a_conc_temp = a_conc1 + delta_conca*relposa + delta_conca2*relposb;
		b_conc_temp = b_conc1 + delta_concb*relposa + delta_concb2*relposb;
		a_pot_temp = a_pot1 + delta_pota*relposa + delta_pota2*relposb;
		b_pot_temp = b_pot1 + delta_potb*relposa + delta_potb2*relposb;
		  
// 		      if(relposa <= 0)
// 		      {
// 			a_conc_temp = a_conc1;
// 			a_pot_temp = a_pot1;
// 			b_conc_temp = b_conc1;
// 			b_pot_temp = b_pot1;
// 		      }
// 		      
// 		      if(relpos > 0 && relpos < 1)
// 		      {
// 			a_conc_temp = a_conc1 + delta_conca*relpos;
// 			b_conc_temp = b_conc1 + delta_concb*relpos;
// 			a_pot_temp = a_pot1 + delta_pota*relpos;
// 			b_pot_temp = b_pot1 + delta_potb*relpos;
// 		      }
// 		      
// 			if(relpos >=1)
// 		      {
// 			a_conc_temp = a_conc2;
// 			a_pot_temp = a_pot2;
// 			b_conc_temp = b_conc2;
// 			b_pot_temp = b_pot2;
// 		      }
		  
		  
		  
		  
		
		
		
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
		//printf("%d	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", l, site_coords[l][1], relposa, relposb, site_pots[l], a_pot_temp, b_pot_temp, a_conc_temp, b_conc_temp);
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
                          (SiteArray->pert_pos) = NULL;

	
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
// 		printf("%d	%lf	%d	%lf\n", l, temprandnum, siteinfo[l][1], site_pots[l]);
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
			      if(siteinfo[l][0] == 0  && siteinfo[l][1] == 0 && site_pots[l] == a_pot)
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
			      if(siteinfo[l][0] == 0   && siteinfo[l][1] == 1 && site_pots[l] == b_pot)
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
                          (SiteArray->pert_pos) = NULL;

	
}

void genSublatticeMoire(RectRedux *SiteArray, void *p, int struc_out, char *filename)
{  
	  submoire_para *params = (submoire_para *)p;
		double AA_mass= (params->AA_mass);
		double AA_pot= (params->AA_pot);
		double AB_mass= (params->AB_mass);
		double AB_pot= (params->AB_pot);
		double BA_mass= (params->BA_mass);
		double BA_pot= (params->BA_pot);
		int lM= (params->lM);
		double x_offset= (params->x_offset);
		double y_offset= (params->y_offset);
		int seed= (params->seed);
	  
		int i, j, k;
	 
	  
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
		
		
	      
	//potentials and masses stuff here
		//high symmetry AA, AB and BA points needed for the rectangular moire unit cell
			//these are used to divide the cell into triangles
			//the mass and potential at a point is interpolated from those at the high symm points 
			//at the vertices of its containing triangles
			
		int i1, i2; 
		if(geo==0)
		{
			i1=0; i2=1;
		}
		if(geo==1)
		{
			i1=1; i2=0;
		}
		
		double AAtriag[5][2], ABtriag[4][2], BAtriag[4][2];
		int celltriag[14][3];
		
		AAtriag[0][i1] = 0.0, 	AAtriag[0][i2] = 0.0;
		AAtriag[1][i1] = 0.0, 	AAtriag[1][i2] = lM*sqrt(3);
		AAtriag[2][i1] = lM/2.0, 	AAtriag[2][i2] = lM*sqrt(3)/2.0;
		AAtriag[3][i1] = lM, 	AAtriag[3][i2] = 0.0;
		AAtriag[4][i1] = lM, 	AAtriag[4][i2] = lM*sqrt(3.0);
		
		ABtriag[0][i1] = 0.0, 	ABtriag[0][i2] = lM/sqrt(3.0);
		ABtriag[1][i1] = lM/2.0, 	ABtriag[1][i2] = -lM/(2*sqrt(3.0));
		ABtriag[2][i1] = lM/2.0, 	ABtriag[2][i2] = 5.0*lM/(2*sqrt(3.0));
		ABtriag[3][i1] = lM, 	ABtriag[3][i2] = lM/sqrt(3.0);
		
		BAtriag[0][i1] = 0.0, 	 BAtriag[0][i2] = 2.0*lM/sqrt(3.0);
		BAtriag[1][i1] = lM/2.0, BAtriag[1][i2] = lM/(2.0*sqrt(3.0));
		BAtriag[2][i1] = lM/2.0, BAtriag[2][i2] = 7.0*lM/(2*sqrt(3.0));
		BAtriag[3][i1] = lM, BAtriag[3][i2] = 2.0*lM/(sqrt(3.0));
		
		celltriag[0][0] = 0; 	celltriag[0][1] = 1; 	celltriag[0][2]=1;
		celltriag[1][0] = 3; 	celltriag[1][1] = 1; 	celltriag[1][2]=1;
		celltriag[2][0] = 0; 	celltriag[2][1] = 0; 	celltriag[2][2]=1;
		celltriag[3][0] = 3; 	celltriag[3][1] = 3; 	celltriag[3][2]=1;
		celltriag[4][0] = 2; 	celltriag[4][1] = 0; 	celltriag[4][2]=1;
		celltriag[5][0] = 2; 	celltriag[5][1] = 3; 	celltriag[5][2]=1;
		celltriag[6][0] = 2; 	celltriag[6][1] = 0; 	celltriag[6][2]=0;
		celltriag[7][0] = 2; 	celltriag[7][1] = 3; 	celltriag[7][2]=3;
		celltriag[8][0] = 2; 	celltriag[8][1] = 2; 	celltriag[8][2]=0;
		celltriag[9][0] = 2; 	celltriag[9][1] = 2; 	celltriag[9][2]=3;
		celltriag[10][0] = 1; 	celltriag[10][1] = 2; 	celltriag[10][2]=0;
		celltriag[11][0] = 4; 	celltriag[11][1] = 2; 	celltriag[11][2]=3;
		celltriag[12][0] = 1; 	celltriag[12][1] = 2; 	celltriag[12][2]=2;
		celltriag[13][0] = 4; 	celltriag[13][1] = 2; 	celltriag[13][2]=2;
                                           

// 		for(l=0; l<14; l++)
// 		{
// 			printf("%lf	%lf\n", AAtriag[celltriag[l][0]][0],  AAtriag[celltriag[l][0]][1]);
// 			printf("%lf	%lf\n", ABtriag[celltriag[l][1]][0],  ABtriag[celltriag[l][1]][1]);
// 			printf("%lf	%lf\n", BAtriag[celltriag[l][2]][0],  BAtriag[celltriag[l][2]][1]);
// 			printf("%lf	%lf\n\n", AAtriag[celltriag[l][0]][0],  AAtriag[celltriag[l][0]][1]);
// 		}
// 		printf("\n");
// 		exit(0);
	
		//for each site, determine it's position in the rectangular moire unit cell
		double *polyx = createDoubleArray(3);
		double *polyy = createDoubleArray(3);
		double *temppolyx=createDoubleArray(3);
		double *temppolyy = createDoubleArray(3);
		double areaT, areaAA, areaAB, areaBA;
		double siteposx, siteposy;
		int mx=0, my=0;
		double *masses = createDoubleArray(tot_sites);
		double *potentials = createDoubleArray(tot_sites);
		double norm;
		
		for(i=0; i<*Ntot; i++)
		{
			//Transform coordinate into main unit cell
			siteposx = site_coords[i][0] - x_offset;
			siteposy = site_coords[i][1] - y_offset;
			
			if(geo==0)
			{
				mx = (int) (siteposx / lM);
				my = (int) (siteposy / (lM*sqrt(3)));
				
				if( mx != 0 )
				{
					siteposx = siteposx - mx * lM;
				}
				if( my != 0 )
				{
					siteposy = siteposy - my * lM * sqrt(3);
				}
				
				
				
			}
			
			if(geo==1)
			{
				my = (int) (siteposy / lM);
				mx = (int) (siteposx / (lM*sqrt(3)));
				
				if( my != 0 )
				{
					siteposy = siteposy - my * lM;
				}
				if( mx != 0 )
				{
					siteposx = siteposx - mx * lM * sqrt(3);
				}
				
				
			}
			
						
			
			//Loop over possible triangles in unit cell
			for(j=0; j<14; j++)
			{
//  			j=2;
				polyx[0] = AAtriag[celltriag[j][0]][0];
				polyy[0] = AAtriag[celltriag[j][0]][1];
				polyx[1] = ABtriag[celltriag[j][1]][0];
				polyy[1] = ABtriag[celltriag[j][1]][1];
				polyx[2] = BAtriag[celltriag[j][2]][0];
				polyy[2] = BAtriag[celltriag[j][2]][1];
				
				
				if(pntriangle(polyx, polyy, siteposx, siteposy))
				{
// 					areaT = areaTriangle(polyx, polyy);
					
					temppolyx[0] = siteposx;
					temppolyy[0] = siteposy;
					temppolyx[1] = polyx[1];
					temppolyy[1] = polyy[1];
					temppolyx[2] = polyx[2];
					temppolyy[2] = polyy[2];
					areaAA = areaTriangle(temppolyx, temppolyy);
					
					temppolyx[0] = polyx[0];
					temppolyy[0] = polyy[0];
					temppolyx[1] = siteposx;
					temppolyy[1] = siteposy;
					temppolyx[2] = polyx[2];
					temppolyy[2] = polyy[2];
					areaAB = areaTriangle(temppolyx, temppolyy);
					
					temppolyx[0] = polyx[0];
					temppolyy[0] = polyy[0];
					temppolyx[1] = polyx[1];
					temppolyy[1] = polyy[1];
					temppolyx[2] = siteposx;
					temppolyy[2] = siteposy;
					areaBA = areaTriangle(temppolyx, temppolyy);
					
					norm= areaAA*areaAA + areaAB*areaAB+ areaBA*areaBA;
					
 					masses[i] = (areaAA*areaAA*AA_mass + areaAB*areaAB*AB_mass + areaBA*areaBA*BA_mass)/norm;
					potentials[i] = (areaAA*areaAA*AA_pot + areaAB*areaAB*AB_pot + areaBA*areaBA*BA_pot)/norm;
					
					
				}
			}
			
			if(siteinfo[i][1] ==0)
				site_pots[i] = potentials[i] + 0.5*masses[i];
			
			if(siteinfo[i][1] ==1)
				site_pots[i] = potentials[i] - 0.5*masses[i];
			
						
		}
		//printf("\n");
// 		for(i=0; i<*Ntot; i++)
// 		{
// 			printf("%lf	%lf	%e\n", site_coords[i][0], site_coords[i][1], masses[i]);
// 		}
// 		      
// 		     exit(0);
		
		free(masses);
		free(potentials);
		
		
		
		
	      double temprandnum;
	      
	      
	      	      
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
				for(l=0; l<*Ntot; l++)
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
                          (SiteArray->pert_pos) = NULL;
			  
			  
}



//generates averaged sublattice potentials in left right leads (leads 0 and 1)
//this is useful for injection in and out of infinite devices
void genSublatticeLeadPots(RectRedux **Leads, void *p)
{
	int i, j;
	subl_para *params = (subl_para *)p;
	double a_conc = (params->a_conc);
	double b_conc = (params->b_conc);
	double a_pot = (params->a_pot);
	double b_pot = (params->b_pot);
	
	
	for(i=0; i<2; i++)
	{
		for(j=0; j< *(Leads[i]->Ntot); j++)
		{
			if( (Leads[i]->siteinfo)[j][0] == 0 && (Leads[i]->siteinfo)[j][1] == 0)
			{
				(Leads[i]->site_pots)[j] = a_conc * a_pot;
			}
			
			if( (Leads[i]->siteinfo)[j][0] == 0 && (Leads[i]->siteinfo)[j][1] == 1)
			{
				(Leads[i]->site_pots)[j] = b_conc * b_pot;
			}
		}
	}
	
	
}

//generates symmetric strain in leads similar to device
void customLeadStrain(RectRedux **Leads, void *p)
{
	int i, j, k;
	symstrain_para *params = (symstrain_para *)p;
	
        char *straingeo = (params->straingeo);
        double strain_mag = (params->strain_mag);
        double strain_width = (params->strain_width);
        int location = (params->location);
        
        int num_feat =0;
            if(strcmp("gaussfold", straingeo) == 0)
            {
                if(location == 0)
                {
                    num_feat =1;
                }
                if(location == 1)
                {
                    num_feat =2;
                }
            }
        
            double *featy = createDoubleArray(num_feat);
            
            if( (Leads[0]->geo) ==0)
            {
                    if(location==0)
                    {
                        featy[0] = (Leads[0]->length) * sqrt(3) /4 ;
                        
                    }
                    
                    if(location==1)
                    {
                        featy[0] = 1 / (2*sqrt(3));
                        featy[1] = (Leads[0]->length) * sqrt(3)/2 -  1 / (2*sqrt(3));
                    }
            }

            if((Leads[0]->geo)==1)
            {
                    if(location==0)
                    {
                        featy[0] = ((Leads[0]->length)-1) * 0.25;
                    }
                    
                    if(location==1)
                    {
                        featy[0] = 0.0;
                        featy[1] = ((Leads[0]->length)-1) * 0.5;
                    }
            }
	  
	double effx, effy, ux, uy, uz, u0, u1, u2;
        
	for(i=0; i<2; i++)
	{
                (Leads[i]->pert_pos) = createNonSquareDoubleMatrix(*(Leads[i]->Ntot), 3);
		for(j=0; j< *(Leads[i]->Ntot); j++)
		{
                    (Leads[i]->pert_pos)[j][0] = (Leads[i]->pos)[j][0];
                    (Leads[i]->pert_pos)[j][1] = (Leads[i]->pos)[j][1];
                    (Leads[i]->pert_pos)[j][2] = (Leads[i]->pos)[j][2];
                    
                    
                    for(k=0; k<num_feat; k++)
                    {
                        effy = (Leads[i]->pert_pos)[j][1] - featy[k];
                        uz=0; ux=0; uy=0;
                        
                        
                        //Gaussian fold
                            if(strcmp("gaussfold", straingeo) == 0)
                            {
                                //allow deformation out to 3*sigma
                                if (effy < 3.0*strain_width)
                                {
                                    uz = strain_mag * exp (- effy*effy / (strain_width*strain_width) );
                                }
                            }
                            
                            (Leads[i]->pert_pos)[j][0] += ux;
                            (Leads[i]->pert_pos)[j][1] += uy;
                            (Leads[i]->pert_pos)[j][2] += uz;
                            
                        
                        
                    }
                    
			
		}
	}
	
	
}



//A general antidot barrier-type device 
//(circular ALs in triangular lattice for the moment - should be generalised later -- HAS BEEN DONE!)
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
          int orientation = (params->orientation);

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
	  if(strcmp("triAC", dotgeo) == 0)
	  {
	    rad_per_hole=3;
	    renorm_fac = 1/ (2*sqrt(3)) ;
	  }
	  if(strcmp("triZZ", dotgeo) == 0)
	  {
	    rad_per_hole=3;
	    renorm_fac = 1/ (2*sqrt(3));
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
          double temprand;
          int tempsign=1;
          
          int *hole_int = createIntArray(tot_holes);
          
            for (i=0; i< tot_holes; i++)
            {
                if(orientation ==0)
                    hole_int[i] =1;
                            
                if(orientation ==1)
                    hole_int[i] = -1;
                        
                if(orientation ==2)
                {
                    temprand = myRandNum(-1.0, 1.0);
                    
                    if(temprand  >=0)
                        hole_int[i] =1;
                    else
                        hole_int[i] =-1;
                }
            }
          
            
	 
                      
	  
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
			
			if(isperiodic==1)
			{
			  if(sqrt( pow( site_coords[i][0] - holes[j][0], 2) + pow( site_coords[i][1] - (holes[j][1]-cellyshift), 2)) < holes[j][2])
			  {
			    sites[i][0] = 1;
			  }
			  
			  if(sqrt( pow( site_coords[i][0] - holes[j][0], 2) + pow( site_coords[i][1] - (holes[j][1]+cellyshift), 2)) < holes[j][2])
			  {
			    sites[i][0] = 1;
			  }
			  
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
			  
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 1;
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 1;
			}
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
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 1;
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 1;
			}
		      }
		      
		      
		      
		    }
		    
                    if(strcmp("triZZ", dotgeo) == 0)
		    {
                        
                        
                        tempsign=hole_int[j];
                        
		      if(geo == 0)
		      {
			polyx[0] = holes[j][0];
			polyy[0] = holes[j][1] + tempsign*3*holes[j][2]/2;
                        
			polyx[1] = holes[j][0] - sqrt(3) * holes[j][2];
			polyy[1] = holes[j][1] - tempsign*3*holes[j][2]/2;
                        
			polyx[2] = holes[j][0] + sqrt(3) * holes[j][2];
			polyy[2] = holes[j][1] - tempsign*3*holes[j][2]/2;
			
		      }
		      
		      if(geo == 1)
		      {
			polyx[0] = holes[j][0] + tempsign*3*holes[j][2]/2;
			polyy[0] = holes[j][1];
                        
			polyx[1] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[1] = holes[j][1] - sqrt(3) * holes[j][2];
                        
			polyx[2] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[2] = holes[j][1] + sqrt(3) * holes[j][2];
			
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
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 1;
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 1;
			}
		      }
		      
		      
		      
		    }
		    
		    
		    if(strcmp("triAC", dotgeo) == 0)
		    {
                        tempsign=hole_int[j];
                        
                        
                        
		      if(geo == 0)
		      {
                          			                       
                        polyx[0] = holes[j][0] + tempsign*3*holes[j][2]/2;
			polyy[0] = holes[j][1];
                        
			polyx[1] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[1] = holes[j][1] - sqrt(3) * holes[j][2];
                        
			polyx[2] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[2] = holes[j][1] + sqrt(3) * holes[j][2];
			
		      }
		      
		      if(geo == 1)
		      {
			polyx[0] = holes[j][0];
			polyy[0] = holes[j][1] + tempsign*3*holes[j][2]/2;
                        
			polyx[1] = holes[j][0] - sqrt(3) * holes[j][2];
			polyy[1] = holes[j][1] - tempsign*3*holes[j][2]/2;
                        
			polyx[2] = holes[j][0] + sqrt(3) * holes[j][2];
			polyy[2] = holes[j][1] - tempsign*3*holes[j][2]/2;
			
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
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 1;
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 1;
			}
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
		      
		         if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 1;
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 1;
			}
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
			  
				    
			    free(hole_int);
			
			  
			  //Fill the array of data structures describing the system
			
			  
			  if(struc_out != 0)
			  {
			    fclose(out);
			    
			  }
			  
			  
			  (SiteArray->pos) = site_coords;
			  (SiteArray->site_pots) = site_pots;
			  (SiteArray->siteinfo) = siteinfo;
                          (SiteArray->pert_pos) = NULL;

}


//Inverse of antidots -- creates an array of shaped regions
//On its own, will probably crash the transport routines due to lack of connectivity
//Intended to be used as part of composite multilayer systems
void genDotDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename)
{  
	  adot_para *params = (adot_para *)p;
	  
	  int buffer_rows = (params->buffer_rows);
	  int AD_length = (params->AD_length);
	  int AD_length2 = (params->AD_length2);
          int orientation = (params->orientation);

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
	  if(strcmp("triAC", dotgeo) == 0)
	  {
	    rad_per_hole=3;
	    renorm_fac = 1/ (2*sqrt(3)) ;
	  }
	  if(strcmp("triZZ", dotgeo) == 0)
	  {
	    rad_per_hole=3;
	    renorm_fac = 1/ (2*sqrt(3));
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
          double temprand;
          int tempsign=1;
          
          int *hole_int = createIntArray(tot_holes);
          
            for (i=0; i< tot_holes; i++)
            {
                if(orientation ==0)
                    hole_int[i] =1;
                            
                if(orientation ==1)
                    hole_int[i] = -1;
                        
                if(orientation ==2)
                {
                    temprand = myRandNum(-1.0, 1.0);
                    
                    if(temprand  >=0)
                        hole_int[i] =1;
                    else
                        hole_int[i] =-1;
                }
            }
          
            
	 
                      
	  
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
          //INITIALLY ALL ARE REMOVED!
          for(i=0; i<tot_sites; i++)
          {
              sites[i][0] = 1;
          }
          
	  //double **vertices = createNonSquareDoubleMatrix(rad_per_hole, 2);
	  double *polyx = createDoubleArray(rad_per_hole);
	  double *polyy = createDoubleArray(rad_per_hole);
   
	   
	  //atom removal!
	      for(i=0; i< tot_sites ; i++)
	      {
		for(j=0; j< tot_holes; j++)
		{
		  
		    if(strcmp("circ", dotgeo) == 0)
		    {
			if(sqrt( pow( site_coords[i][0] - holes[j][0], 2) + pow( site_coords[i][1] - holes[j][1], 2)) < holes[j][2])
			{
			  sites[i][0] = 0;
			}
			
			if(isperiodic==1)
			{
			  if(sqrt( pow( site_coords[i][0] - holes[j][0], 2) + pow( site_coords[i][1] - (holes[j][1]-cellyshift), 2)) < holes[j][2])
			  {
			    sites[i][0] = 0;
			  }
			  
			  if(sqrt( pow( site_coords[i][0] - holes[j][0], 2) + pow( site_coords[i][1] - (holes[j][1]+cellyshift), 2)) < holes[j][2])
			  {
			    sites[i][0] = 0;
			  }
			  
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
			sites[i][0] = 0;
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
			  
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 0;
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 0;
			}
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
			sites[i][0] = 0;
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 0;
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 0;
			}
		      }
		      
		      
		      
		    }
		    
                    if(strcmp("triZZ", dotgeo) == 0)
		    {
                        
                        
                        tempsign=hole_int[j];
                        
		      if(geo == 0)
		      {
			polyx[0] = holes[j][0];
			polyy[0] = holes[j][1] + tempsign*3*holes[j][2]/2;
                        
			polyx[1] = holes[j][0] - sqrt(3) * holes[j][2];
			polyy[1] = holes[j][1] - tempsign*3*holes[j][2]/2;
                        
			polyx[2] = holes[j][0] + sqrt(3) * holes[j][2];
			polyy[2] = holes[j][1] - tempsign*3*holes[j][2]/2;
			
		      }
		      
		      if(geo == 1)
		      {
			polyx[0] = holes[j][0] + tempsign*3*holes[j][2]/2;
			polyy[0] = holes[j][1];
                        
			polyx[1] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[1] = holes[j][1] - sqrt(3) * holes[j][2];
                        
			polyx[2] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[2] = holes[j][1] + sqrt(3) * holes[j][2];
			
		      }
		      
		      
		      if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
		      {
			sites[i][0] = 0;
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 0;
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 0;
			}
		      }
		      
		      
		      
		    }
		    
		    
		    if(strcmp("triAC", dotgeo) == 0)
		    {
                        tempsign=hole_int[j];
                        
                        
                        
		      if(geo == 0)
		      {
                          			                       
                        polyx[0] = holes[j][0] + tempsign*3*holes[j][2]/2;
			polyy[0] = holes[j][1];
                        
			polyx[1] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[1] = holes[j][1] - sqrt(3) * holes[j][2];
                        
			polyx[2] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[2] = holes[j][1] + sqrt(3) * holes[j][2];
			
		      }
		      
		      if(geo == 1)
		      {
			polyx[0] = holes[j][0];
			polyy[0] = holes[j][1] + tempsign*3*holes[j][2]/2;
                        
			polyx[1] = holes[j][0] - sqrt(3) * holes[j][2];
			polyy[1] = holes[j][1] - tempsign*3*holes[j][2]/2;
                        
			polyx[2] = holes[j][0] + sqrt(3) * holes[j][2];
			polyy[2] = holes[j][1] - tempsign*3*holes[j][2]/2;
			
		      }
		      
		      
		      if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
		      {
			sites[i][0] = 0;
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 0;
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 0;
			}
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
			sites[i][0] = 0;
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
		      
		         if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 0;
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
			  sites[i][0] = 0;
			}
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
		    for(i=0; i<tot_sites; i++)
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
			  
				    
			    free(hole_int);
			
			  
			  //Fill the array of data structures describing the system
			
			  
			  if(struc_out != 0)
			  {
			    fclose(out);
			    
			  }
			  
			  
			  (SiteArray->pos) = site_coords;
			  (SiteArray->site_pots) = site_pots;
			  (SiteArray->siteinfo) = siteinfo;
                          (SiteArray->pert_pos) = NULL;

}





//A general antidot barrier-type device 
//(circular ALs in triangular lattice for the moment - should be generalised later -- HAS BEEN DONE!)
//Different to the routine used in the disordered antidot paper
//Modified to work with generic ribbons devices
//Layout based on sublattice device routine
//The first and last buffer_rows chains of the device will not be altered from pristine graphene
//This version requires more aligning of the graphene and antidot sheets
//i.e. ribbon width is not automatically an integer multiple of GAL cell width
void genSublatticeDots(RectRedux *SiteArray, void *p, int struc_out, char *filename)
{  
	  sldot_para *params = (sldot_para *)p;
	  
	  int buffer_rows = (params->buffer_rows);
	  int AD_length = (params->SD_length);
	  int AD_length2 = (params->SD_length2);
	  double AD_rad = (params->SD_rad);
	  double AD_rad2 = (params->SD_rad2);
          int orientation = (params->orientation);

	  int lat_width = (params->lat_width);
	  int lat_length = (params->lat_length);
	  char *latgeo = (params->latgeo);
	  char *dotgeo = (params->dotgeo);
	  int isperiodic = (params->isperiodic);
	  
	  double radfluc = (params->radfluc);
	  double xyfluc = (params->xyfluc);
          
          double a_conc = (params->a_conc);
	  double a_pot = (params->a_pot);
          double b_conc = (params->b_conc);
	  double b_pot = (params->b_pot);


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
	  double xstart, ystart, temprandnum;
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
	  if(strcmp("triAC", dotgeo) == 0)
	  {
	    rad_per_hole=3;
	    renorm_fac = 1/ (2*sqrt(3)) ;
	  }
	  if(strcmp("triZZ", dotgeo) == 0)
	  {
	    rad_per_hole=3;
	    renorm_fac = 1/ (2*sqrt(3));
	  }
	  
	  
	  
	  tot_holes = holes_per_cell*lat_width*lat_length;
	  double holes[tot_holes][2+rad_per_hole];
	  double unitholes[holes_per_cell][2];
	  double xshift, yshift;
	  
	  double cellyshift;
          double temprand;
          int tempsign=1;
          
          int *hole_int = createIntArray(tot_holes);
          
            for (i=0; i< tot_holes; i++)
            {
                if(orientation ==0)
                    hole_int[i] =1;
                            
                if(orientation ==1)
                    hole_int[i] = -1;
                        
                if(orientation ==2)
                {
                    temprand = myRandNum(-1.0, 1.0);
                    
                    if(temprand  >=0)
                        hole_int[i] =1;
                    else
                        hole_int[i] =-1;
                }
            }
	  
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
	  
	  
// 	  int **sites = createNonSquareIntMatrix(tot_sites, 2); //removed, num neighbours after first sweep
	  //double **vertices = createNonSquareDoubleMatrix(rad_per_hole, 2);
	  double *polyx = createDoubleArray(rad_per_hole);
	  double *polyy = createDoubleArray(rad_per_hole);
	  
          

          

          
          
	   
	  //atom removal (from antidot code) replaced by sublattice dependent potentials
	      for(i=2*length; i< tot_sites - 2*length ; i++)
	      {
		for(j=0; j< tot_holes; j++)
		{
		  
		    if(strcmp("circ", dotgeo) == 0)
		    {
			if(sqrt( pow( site_coords[i][0] - holes[j][0], 2) + pow( site_coords[i][1] - holes[j][1], 2)) < holes[j][2])
			{
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
			}
			
			if(isperiodic==1)
			{
			  if(sqrt( pow( site_coords[i][0] - holes[j][0], 2) + pow( site_coords[i][1] - (holes[j][1]-cellyshift), 2)) < holes[j][2])
			  {
			    site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
			  }
			  
			  if(sqrt( pow( site_coords[i][0] - holes[j][0], 2) + pow( site_coords[i][1] - (holes[j][1]+cellyshift), 2)) < holes[j][2])
			  {
			    site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
			  }
			  
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
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
			  
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
			}
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
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
			}
		      }
		      
		      
		      
		    }
		    
		    if(strcmp("triZZ", dotgeo) == 0)
		    {
                        
                        
                    tempsign=hole_int[j];
                        
		      if(geo == 0)
		      {
			polyx[0] = holes[j][0];
			polyy[0] = holes[j][1] + tempsign*3*holes[j][2]/2;
                        
			polyx[1] = holes[j][0] - sqrt(3) * holes[j][2];
			polyy[1] = holes[j][1] - tempsign*3*holes[j][2]/2;
                        
			polyx[2] = holes[j][0] + sqrt(3) * holes[j][2];
			polyy[2] = holes[j][1] - tempsign*3*holes[j][2]/2;
			
		      }
		      
		      if(geo == 1)
		      {
			polyx[0] = holes[j][0] + tempsign*3*holes[j][2]/2;
			polyy[0] = holes[j][1];
                        
			polyx[1] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[1] = holes[j][1] - sqrt(3) * holes[j][2];
                        
			polyx[2] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[2] = holes[j][1] + sqrt(3) * holes[j][2];
			
		      }
		      
		      
		      if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
		      {
			
                        site_pots[i] = 0.0;
                        temprandnum = myRandNum(0.0, 1.0);
            
                        if(siteinfo[i][1] == 0)
                        {
                            if(temprandnum < a_conc)
                                site_pots[i] = a_pot;
                        }
                        
                        if(siteinfo[i][1] == 1)
                        {
                            if(temprandnum < b_conc)
                                site_pots[i] = b_pot;
                        }
                          
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			
			  if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
                            {
                                site_pots[i] = 0.0;
                                temprandnum = myRandNum(0.0, 1.0);
                    
                                if(siteinfo[i][1] == 0)
                                {
                                    if(temprandnum < a_conc)
                                        site_pots[i] = a_pot;
                                }
                                
                                if(siteinfo[i][1] == 1)
                                {
                                    if(temprandnum < b_conc)
                                        site_pots[i] = b_pot;
                                }
                            }
			
		      }
		      
		      
		      
		    }
		    
		    
		    if(strcmp("triAC", dotgeo) == 0)
		    {
                        
                        
                    tempsign=hole_int[j];
                        
		      if(geo == 0)
		      {
                          			                       
                        polyx[0] = holes[j][0] + tempsign*3*holes[j][2]/2;
			polyy[0] = holes[j][1];
                        
			polyx[1] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[1] = holes[j][1] - sqrt(3) * holes[j][2];
                        
			polyx[2] = holes[j][0] - tempsign*3*holes[j][2]/2;
			polyy[2] = holes[j][1] + sqrt(3) * holes[j][2];
			
		      }
		      
		      if(geo == 1)
		      {
			polyx[0] = holes[j][0];
			polyy[0] = holes[j][1] + tempsign*3*holes[j][2]/2;
                        
			polyx[1] = holes[j][0] - sqrt(3) * holes[j][2];
			polyy[1] = holes[j][1] - tempsign*3*holes[j][2]/2;
                        
			polyx[2] = holes[j][0] + sqrt(3) * holes[j][2];
			polyy[2] = holes[j][1] - tempsign*3*holes[j][2]/2;
			
		      }
		      
		      
		      
		      if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
		      {
			
                        site_pots[i] = 0.0;
                        temprandnum = myRandNum(0.0, 1.0);
            
                        if(siteinfo[i][1] == 0)
                        {
                            if(temprandnum < a_conc)
                                site_pots[i] = a_pot;
                        }
                        
                        if(siteinfo[i][1] == 1)
                        {
                            if(temprandnum < b_conc)
                                site_pots[i] = b_pot;
                        }
                          
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
		      
		      if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			
			  if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
                            {
                                site_pots[i] = 0.0;
                                temprandnum = myRandNum(0.0, 1.0);
                    
                                if(siteinfo[i][1] == 0)
                                {
                                    if(temprandnum < a_conc)
                                        site_pots[i] = a_pot;
                                }
                                
                                if(siteinfo[i][1] == 1)
                                {
                                    if(temprandnum < b_conc)
                                        site_pots[i] = b_pot;
                                }
                            }
			
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
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
		      }
		      
		      if(struc_out == 1 && i== 2*length)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  fprintf(out, "%lf	%lf\n", polyx[k], polyy[k]);
			}
			fprintf(out, "%lf	%lf\n\n", polyx[0], polyy[0]);
		      }
		      
		         if(isperiodic==1)
		      {
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] - cellyshift;
			}
			
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
			}
		      
			for(k=0; k<rad_per_hole; k++)
			{
			  polyy[k] = polyy[k] + 2*cellyshift;
			}
		      
			if(pnpoly(rad_per_hole, polyx, polyy, site_coords[i][0], site_coords[i][1]))
			{
                            site_pots[i] = 0.0;
                            temprandnum = myRandNum(0.0, 1.0);
		
                            if(siteinfo[i][1] == 0)
                            {
                                if(temprandnum < a_conc)
                                    site_pots[i] = a_pot;
                            }
                            
                            if(siteinfo[i][1] == 1)
                            {
                                if(temprandnum < b_conc)
                                    site_pots[i] = b_pot;
                            }
			}
		      }
		      
		    }
		    
		}
	      }
	      
	      free(polyx); free(polyy);
	  


		    
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
			  
				    
			    
			
			  
			  //Fill the array of data structures describing the system
			
			  
			  if(struc_out != 0)
			  {
			    fclose(out);
			    
			  }
			  
			  
			  (SiteArray->pos) = site_coords;
			  (SiteArray->site_pots) = site_pots;
			  (SiteArray->siteinfo) = siteinfo;
                          (SiteArray->pert_pos) = NULL;

}








void genBubbleDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename)
{  
        bubble_para *params = (bubble_para *)p;
        
        int buffer_rows = (params->buffer_rows);
        double bub_rad = (params->bub_rad);
        double bub_height = (params->bub_height);
        int lat_width = (params->lat_width);
        int lat_length = (params->lat_length);
        int cell_length = (params->cell_length);
        int cell_length2 = (params->cell_length2);
        char *latgeo = (params->latgeo);
        char *bubgeo = (params->bubgeo);
        int isperiodic = (params->isperiodic);
        
        double radfluc = (params->radfluc);
        double xyfluc = (params->xyfluc);
        double zfluc = (params->zfluc);
        double ifluc = (params->ifluc);

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
        double **pert_coords = createNonSquareDoubleMatrix(tot_sites, 3);
        double *site_pots = createDoubleArray(tot_sites);
        int **siteinfo = createNonSquareIntMatrix(tot_sites, 2);
        double xstart, ystart;
        int isclean, l, m, tempint, tempint2;
        int *Nrem = (SiteArray->Nrem);
        int i, j, k;
        int *Ntot = (SiteArray->Ntot);
	  
        
        
        
        *Ntot = tot_sites;
	  
        
        //BASE XY GEOMETRIES
        
            //zigzag ribbon
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
            
            //armchair ribbon
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
	  
        //Centre of bubble locations and radii (accounts for possible fluctuations)
        
            //Based on antidot locations for different lattice types ("holes")
        
            int holes_per_cell, tot_holes, rad_per_hole;
            double renorm_fac=1.0; //how far outside "radius" to consider, e.g. in Gaussian types
            
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
            
            //how many "radius" parameters are needed (probably 1 for most bubble types!)
            //(possibly Gaussian needs cut-off...?)
            rad_per_hole = 1;
            
            //params for different bubble geometries can be set here if they involve some renormalisation
            
            if(strcmp("triaxgauss", bubgeo) == 0)
            {
                renorm_fac=3.0;
            }
            if(strcmp("gaussbump", bubgeo) == 0)
            {
                renorm_fac = 3.0;
            }
            
	  
	  
	  
            tot_holes = holes_per_cell*lat_width*lat_length;
            double holes[tot_holes][3+rad_per_hole];
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
                        xshift = 3.0*cell_length;
                        yshift = sqrt(3)*cell_length;
                        
                        unitholes[0][0]  = (int) (cell_length*3.0/4) -0.5* ( ((int) (cell_length/2)) % 2 ); 
                        unitholes[0][1] = (sqrt(3)/2.0) * (int) (cell_length /2.0)   ;
                        unitholes[1][0] = unitholes[0][0] + 1.5*cell_length;
                        unitholes[1][1] = unitholes[0][1] + sqrt(3)*cell_length/2.0;
                        }
                        
                        
                        //armchair ribbon 
                        if(geo==1)
                        {
                        xstart = buffer_rows*sqrt(3);
                        ystart = 0.0;
                        xshift = sqrt(3)*cell_length;
                        yshift = 3.0*cell_length;
                        
                
                        unitholes[0][1] =3.0*cell_length - 0.5 - ((int) (cell_length*3.0/4) -0.5* ( ((int) (cell_length/2)) % 2 ));
                        unitholes[0][0] = (sqrt(3)/2.0) * (int) (cell_length /2.0)  - 1 /(2*sqrt(3)) ;
                        unitholes[1][1] = unitholes[0][1] - 1.5*cell_length;
                        unitholes[1][0] = unitholes[0][0] + sqrt(3)*cell_length/2.0;

                        
                        
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
                        xshift = (cell_length + 1)*1.0;
                        yshift = (cell_length+1)*sqrt(3);
                        
                        unitholes[0][0] = (2*cell_length+1)*0.5 - ((int)(cell_length/4))*1.0 - 0.5; 
                        unitholes[0][1] = (((int)(cell_length/4)) +0.5)*sqrt(3) + sqrt(3)/2;
                        unitholes[1][0] = unitholes[0][0] - (cell_length+1)*0.5  ;
                        unitholes[1][1] = unitholes[0][1] + (cell_length+1)*sqrt(3)*0.5 ;
                        
                        }
                        
                        
                        //armchair ribbon 
                        if(geo==1)
                        {
                        xstart = buffer_rows*sqrt(3);
                        ystart = 0.0;
                        xshift = (cell_length+1)*sqrt(3);
                        yshift = (cell_length + 1)*1.0;
                        
                
                        unitholes[0][1] = ((int)(cell_length/4))*1.0 + 0.5; 
                        unitholes[0][0] = (((int)(cell_length/4)) +0.5)*sqrt(3) + 1/(sqrt(3));;
                        unitholes[1][1] = unitholes[0][1] + (cell_length+1)*0.5  ;
                        unitholes[1][0] = unitholes[0][0] + (cell_length+1)*sqrt(3)*0.5 ;

                        
                        
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
                        xshift = cell_length2;
                        yshift = sqrt(3)*cell_length;
                        
                        unitholes[0][0]  = (int) (cell_length2/2) - 0.5* (  (cell_length % 2) ); 
                        unitholes[0][1] = (sqrt(3)/2.0) * (int) (cell_length)   ;
                        
                        }
                        
                        
                        //armchair ribbon 
                        if(geo==1)
                        {
                        xstart = buffer_rows*sqrt(3);
                        ystart = 0.0;
                        xshift = sqrt(3)*cell_length;
                        yshift = cell_length2;
                        
                
                        unitholes[0][1]  = cell_length2 -0.5 - (int) (cell_length2/2) +0.5* (  (cell_length % 2) ); 
                        unitholes[0][0] = (sqrt(3)/2.0) * (int) (cell_length)  - 1 /(2*sqrt(3)) ;
                        
                                

                        
                        
                        }
                }
	  
	  
	  
	  
	  
            //centre positions, heights, radii
                for(i=0; i<lat_length; i++)
                {
                    for(j=0; j<lat_width; j++)
                    {
                        for(k=0; k <holes_per_cell; k++)
                        {
                            holes[holes_per_cell*i*lat_width + holes_per_cell*j + k][0] = xstart + i*xshift + unitholes[k][0] + myRandNum(-xyfluc, xyfluc);
                            holes[holes_per_cell*i*lat_width + holes_per_cell*j + k][1] = ystart + j*yshift + unitholes[k][1] + myRandNum(-xyfluc, xyfluc);
                            holes[holes_per_cell*i*lat_width + holes_per_cell*j + k][2] = bub_height + myRandNum(-zfluc, zfluc); 
                            
                            for(l=0; l<rad_per_hole; l++)
                            {
                                holes[holes_per_cell*i*lat_width + holes_per_cell*j + k][3+l] = bub_rad * renorm_fac + myRandNum(-radfluc, radfluc);
                            }
                        }
                    }
                }
	  
	
	  
	   
	  //atomic restructuring!
	  //this assumes no overlapping bubbles -- i.e no site is in two bubble regions
	  //be  careful with bubbles with ill-defined radii! (Gaussian-type)
	  
                double effx, effy, effr, effr1, effth, ux, uy, ur, uth, uz, u0, u1, u2;
                
                for(i=0; i< tot_sites ; i++)
                {
                    pert_coords[i][0] = site_coords[i][0];
                    pert_coords[i][1] = site_coords[i][1];
                    pert_coords[i][2] = site_coords[i][2];
                }
          
                for(i=2*length; i< tot_sites - 2*length ; i++)
                {
                    
                    for(j=0; j< tot_holes; j++)
                    {
                        effx = site_coords[i][0] - holes[j][0];
                        effy = site_coords[i][1] - holes[j][1];
                        effr = sqrt(effx*effx + effy*effy);
                        effth = atan2(effy, effx);
                        
                        ur=0; uth=0; uz=0; ux=0; uy=0;
                        
                        
                        //calculate ur and uth for various different bubble types
                        //convert to ux and uy at bottom, and account for periodicity then
                        
                            //triaxial in-plane with gaussian smoothening
                                if(strcmp("triaxgauss", bubgeo) == 0)
                                {
                                    u0=holes[j][2];
                                    //allow deformation out to 3*sigma
                                    if (effr < renorm_fac*holes[j][3])
                                    {
                                        ur = u0 * effr * effr * sin (3*effth) * exp ( - (effr*effr)/(2*holes[j][3]*holes[j][3]));
                                        uth = u0 * effr * effr * sin (3*effth) * exp ( - (effr*effr)/(2*holes[j][3]*holes[j][3]));
                                    }
                                }
                            
                            //Gaussian bump
                                if(strcmp("gaussbump", bubgeo) == 0)
                                {
                                    u0=holes[j][2];
                                    //allow deformation out to 3*sigma
                                    if (effr < renorm_fac*holes[j][3])
                                    {
                                        uz = u0 * exp ( - (effr*effr)/(2*holes[j][3]*holes[j][3]));
                                    }
                                }
                                
                                
                            //Membrane model
                                if(strcmp("membrane", bubgeo) == 0)
                                {
                                    u0= 1.136 * pow(holes[j][2], 2) / holes[j][3];
                                    if (effr < renorm_fac*holes[j][3])
                                    {
                                        ur = u0 * (effr / holes[j][3]) * ( 1.0 -  (effr / holes[j][3])  );
                                        uz = holes[j][2] * (1.0 -  pow((effr / holes[j][3]), 2) ) ;
                                    }
                                }
                                
                             //Non-linear plate model
                                if(strcmp("nlplate", bubgeo) == 0)
                                {
                                    u1= 1.308 * pow(holes[j][2], 2) / pow(holes[j][3], 3);
                                    u2= -1.931 * pow(holes[j][2], 2) / pow(holes[j][3], 4);
                                    if (effr < renorm_fac*holes[j][3])
                                    {
                                        ur = effr * (holes[j][3] - effr) * (u1 + u2 *effr);
                                        uz = holes[j][2] * pow((1.0 -  pow((effr / holes[j][3]), 2) ), 2)  ;
                                    }
                                }   
                            
                            
                            
                            
                                ux = (1.0 + myRandNum(-ifluc, ifluc))*(ur * cos (effth) - uth * sin (effth));
                                uy = (1.0 + myRandNum(-ifluc, ifluc))*(ur * sin (effth) + uth * cos (effth));
                                //uz=(1.0 + myRandNum(-ifluc, ifluc))*uz;
                                
                                pert_coords[i][0] += ux;
                                pert_coords[i][1] += uy;
                                pert_coords[i][2] += uz;
                        
                        
                                
                    
                                //periodic images of bubbles in upper and lower cells, which overlap with this cell
                                if(isperiodic==1)
                                {
                                    effr1 = sqrt(effx*effx + pow(effy+cellyshift,2)   );
                                    if(effr1 < renorm_fac*holes[j][3])
                                    {
                                        pert_coords[i][0] += ux;
                                        pert_coords[i][1] += uy;
                                        pert_coords[i][2] += uz;
                                    }
                                    
                                    effr1 = sqrt(effx*effx + pow(effy-cellyshift,2)   );
                                    if(effr1 < renorm_fac*holes[j][3])
                                    {
                                        pert_coords[i][0] += ux;
                                        pert_coords[i][1] += site_coords[j][1] + uy;
                                        pert_coords[i][2] += uz;
                                    }
                                }
                        }
                    }
                
	      
	      
		    
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
                            for(j=0; j< tot_holes; j++)
                            {
                                fprintf(out, "# BUB %d : %lf	%lf  %lf   %lf\n", j, holes[j][0], holes[j][1], holes[j][2], holes[j][3] );
                            }
                            
                            
			  for(l=0; l<2*length*length2; l++)
			  {
			    if(siteinfo[l][0] == 0)
			      fprintf(out, "%lf	%lf  %lf\n", pert_coords[l][0], pert_coords[l][1], pert_coords[l][2]);
			    
			  }
			  
			}
			  
				    
			    
			
			  
			  //Fill the array of data structures describing the system
			
			  
			  if(struc_out != 0)
			  {
			    fclose(out);
			    
			  }
			  
			  
			  (SiteArray->pos) = site_coords;
                          (SiteArray->pert_pos) = pert_coords;
			  (SiteArray->site_pots) = site_pots;
			  (SiteArray->siteinfo) = siteinfo;


				
	
}

//update customLeadStrain in parallel to this, or else trouble will arise!
void genSymStrain(RectRedux *SiteArray, void *p, int struc_out, char *filename)
{  
        symstrain_para *params = (symstrain_para *)p;
        
        char *straingeo = (params->straingeo);
                
        double strain_mag = (params->strain_mag);
        double strain_width = (params->strain_width);
        int location = (params->location);

        int length = SiteArray->length;
        int length2 = SiteArray->length2;
        int geo = SiteArray->geo;
        
	    
        double smalldist;
        
        FILE *out;
        
        if(struc_out != 0)
        {
            out = fopen(filename, "w");
        }

        //atomic coordinates and the atoms that are in and out, similar to previous cases
        int tot_sites = 2*length*length2;
        double **site_coords = createNonSquareDoubleMatrix(tot_sites, 3);
        double **pert_coords = createNonSquareDoubleMatrix(tot_sites, 3);
        double *site_pots = createDoubleArray(tot_sites);
        int **siteinfo = createNonSquareIntMatrix(tot_sites, 2);
        double xstart, ystart;
        int isclean, l, m, tempint, tempint2;
        int *Nrem = (SiteArray->Nrem);
        int i, j, k;
        int *Ntot = (SiteArray->Ntot);
	  
        *Ntot = tot_sites;
	  
        
        //BASE XY GEOMETRIES
        
            //zigzag ribbon
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
            
            //armchair ribbon
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
	  
        //Strain locations
            
            int num_feat =0;
            if(strcmp("gaussfold", straingeo) == 0)
            {
                if(location == 0)
                {
                    num_feat =1;
                }
                if(location == 1)
                {
                    num_feat =2;
                }
            }
        
            double *featy = createDoubleArray(num_feat);
            if(geo==0)
            {
                    if(location==0)
                    {
                        featy[0] = length * sqrt(3) /4 ;

                       
                    }
                    
                    if(location==1)
                    {
                        featy[0] = 1 / (2*sqrt(3));
                        featy[1] = length * sqrt(3)/2 -  1 / (2*sqrt(3));
                    }
            }

            if(geo==1)
            {
                    if(location==0)
                    {
                        featy[0] = (length-1) * 0.25;
                    }
                    
                    if(location==1)
                    {
                        featy[0] = 0.0;
                        featy[1] = (length-1) * 0.5;
                    }
            }
	  
	  //atomic restructuring!
	  //this assumes no overlapping features -- i.e no site is in two regions
	  //be  careful with features with ill-defined radii! (Gaussian-type)
	  
                double effx, effy, ux, uy, uz, u0, u1, u2;
                
                for(i=0; i< tot_sites ; i++)
                {
                    pert_coords[i][0] = site_coords[i][0];
                    pert_coords[i][1] = site_coords[i][1];
                    pert_coords[i][2] = site_coords[i][2];
                }
          
                for(i=0; i< tot_sites; i++)
                {
                    
                    for(j=0; j< num_feat; j++)
                    {
                        effy = site_coords[i][1] - featy[j];
                        
                        uz=0; ux=0; uy=0;
                        
                        
                        //various feature types
                        //convert to ux and uy at bottom, and account for periodicity then
                        
                            //Gaussian fold
                                if(strcmp("gaussfold", straingeo) == 0)
                                {
                                    //allow deformation out to 3*sigma
                                    if (effy < 3.0*strain_width)
                                    {
                                        uz = strain_mag * exp (- effy*effy / (strain_width*strain_width) );
                                    }
                                }
                            
                                
                                pert_coords[i][0] += ux;
                                pert_coords[i][1] += uy;
                                pert_coords[i][2] += uz;
                    
                        }
                    }
                
	      
	      
		    
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
			      fprintf(out, "%lf	%lf  %lf\n", pert_coords[l][0], pert_coords[l][1], pert_coords[l][2]);
			    
			  }
			  
			}
			  
				    
			    
			
			  
			  //Fill the array of data structures describing the system
			
			  
			  if(struc_out != 0)
			  {
			    fclose(out);
			    
			  }
			  
			  
			  (SiteArray->pos) = site_coords;
                          (SiteArray->pert_pos) = pert_coords;
			  (SiteArray->site_pots) = site_pots;
			  (SiteArray->siteinfo) = siteinfo;


				
	
}



//update customLeadStrain in parallel to this, or else trouble will arise!
void genRandStrain(RectRedux *SiteArray, void *p, int struc_out, char *filename)
{  
        randstrain_para *params = (randstrain_para *)p;
        
                
        double strain_mag = (params->strain_mag);
        double strain_width = (params->strain_width);
        int location = (params->location);
        int buffer_rows = (params->buffer_rows);
        int seed = (params->seed);


        int length = SiteArray->length;
        int length2 = SiteArray->length2;
        int geo = SiteArray->geo;
            srand(time(NULL) + seed);

	    
        double smalldist;
        
        FILE *out;
        
        if(struc_out != 0)
        {
            out = fopen(filename, "w");
        }

        //atomic coordinates and the atoms that are in and out, similar to previous cases
        int tot_sites = 2*length*length2;
        double **site_coords = createNonSquareDoubleMatrix(tot_sites, 3);
        double **pert_coords = createNonSquareDoubleMatrix(tot_sites, 3);
        double *site_pots = createDoubleArray(tot_sites);
        int **siteinfo = createNonSquareIntMatrix(tot_sites, 2);
        double xstart, ystart;
        int isclean, l, m, tempint, tempint2;
        int *Nrem = (SiteArray->Nrem);
        int i, j, k;
        int *Ntot = (SiteArray->Ntot);
	  
        *Ntot = tot_sites;
	  
        
        //BASE XY GEOMETRIES
        
            //zigzag ribbon
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
            
            //armchair ribbon
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
	  
        //Strain locations
            
            int num_feat =0;
           
            if(location == 0)
            {
                num_feat =1;
            }
            if(location == 1)
            {
                num_feat =2;
            }
            
        
            double *featy = createDoubleArray(num_feat);
            if(geo==0)
            {
                                       
                    if(location==1)
                    {
                        featy[0] = 1 / (2*sqrt(3));
                        featy[1] = length * sqrt(3)/2 -  1 / (2*sqrt(3));
                    }
            }

            if(geo==1)
            {
                    
                    if(location==1)
                    {
                        featy[0] = 0.0;
                        featy[1] = (length-1) * 0.5;
                    }
            }
	  
	  //atomic restructuring!
	  //this assumes no overlapping features -- i.e no site is in two regions
	  //be  careful with features with ill-defined radii! (Gaussian-type)
	  
                double ux, uy, uz, effy;
                
                for(i=0; i< tot_sites ; i++)
                {
                    pert_coords[i][0] = site_coords[i][0];
                    pert_coords[i][1] = site_coords[i][1];
                    pert_coords[i][2] = site_coords[i][2];
                }
          
                for(i=0; i< tot_sites; i++)
                {
                    
                    ux=0; uy=0, uz=0;
                    if(location == 0)
                    {
                        ux = myRandNum(-strain_mag, strain_mag);
                        uy = myRandNum(-strain_mag, strain_mag);
                    }
                        
                    if(location ==1)
                    {
                        for(j=0; j< num_feat; j++)
                        {
                                                 
                            effy = site_coords[i][1] - featy[j];
                            
                            if( fabs(effy) < strain_width )
                            {
                                ux = myRandNum(-strain_mag, strain_mag);
                                uy = myRandNum(-strain_mag, strain_mag);
                            }
                        }
                    }
                        
                     
                    pert_coords[i][0] += ux;
                    pert_coords[i][1] += uy;
                    pert_coords[i][2] += uz;
                    
                        
                }
                
	      
	      
		    
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
			      fprintf(out, "%lf	%lf  %lf\n", pert_coords[l][0], pert_coords[l][1], pert_coords[l][2]);
			    
			  }
			  
			}
			  
				    
			    
			
			  
			  //Fill the array of data structures describing the system
			
			  
			  if(struc_out != 0)
			  {
			    fclose(out);
			    
			  }
			  
			  
			  (SiteArray->pos) = site_coords;
                          (SiteArray->pert_pos) = pert_coords;
			  (SiteArray->site_pots) = site_pots;
			  (SiteArray->siteinfo) = siteinfo;


				
	
}





//generates a ribbon device with a combination of smooth and rough edge disorder
//the smooth disorder is generated by a superposition of random sinusoids 
//the rough edge disorder by random vacancy introduction
void genEdgeDisorderedDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename)
{
    edgedis_para *params = (edgedis_para *)p;
  
  
    int buffer_rows = (params->buffer_rows);
    int sruns = (params->sruns);	//number of smooth sinuisoidal fns to superimpose on edge. 
    double smax = (params->smax);	//maximum amplitude of each random sinusoid
    double minper = (params->minper);	//minimum period of the random sinusoids
    int vacruns = (params->vacruns);	//number of random edge vacancy creation runs
    double vacprob = (params->vacprob);	//probability of removing atom in first vacancy creation runs
					//vacprob * (1.0 - l*1.0/vacruns) for lth run
    int seed = (params->seed);
    


    int length = (SiteArray->length);
    int length2 = (SiteArray->length2);
    int geo = (SiteArray->geo);
    
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
    int isclean, l, m, i, j, tempint;
    int *Nrem = (SiteArray->Nrem);
    int *Ntot = (SiteArray->Ntot);
    double topy, bottomy, leftx, rightx, widthrib, lengthrib, smalldist;
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


    
    bottomy = site_coords[0][1];
    topy = site_coords[2*length -1][1];
    leftx=site_coords[2*length*buffer_rows][0];
    rightx = site_coords[2*length*(length2-buffer_rows-1)][0];
    widthrib = topy-bottomy;
    lengthrib = rightx - leftx;
    
    int startloop, endloop;
    int *sites = createIntArray(tot_sites); //  num neighbours of each atom
    int *orig_sites = createIntArray(tot_sites); //  original num neighbours of each atom
    
    //SMOOTH EDGE DISORDER
    
	double *svals = createDoubleArray(2*sruns);
	double *periods = createDoubleArray(2*sruns);
	double *startx = createDoubleArray(2*sruns);
	double *endx = createDoubleArray(2*sruns);
	double *phases = createDoubleArray(2*sruns);
	double tempstart, thisminy, thismaxy;
	
	
	//set sinusoid parameters
	  //2 because of top AND bottom
	    for(i=0; i<2*sruns; i++)
	    {
	      svals[i] = myRandNum(0, smax);
	      periods[i] = myRandNum(minper, 2*(lengthrib));
	      startx[i] = myRandNum(leftx, rightx);
	      endx[i] = myRandNum(leftx, rightx);
	      phases[i] = myRandNum(0, 2*M_PI);
	      
	      if(endx[i] < startx[i])
	      {
		tempstart = endx[i];
		endx[i] = startx[i];
		startx[i] = tempstart;
		
		//not sure what this is.... sets maximum of three periods of any one oscillation?
		if(endx[i] -startx[i] > 3*periods[i])
		{
		  endx[i] = startx[i] + 3*periods[i];
		}
	      }
	    }
	
	     int **neighbours = createNonSquareIntMatrix(tot_sites, 3);
	     int neighindex;
	     double temprandnum;
	//indexes neighbours of each site - Q: why?
	      //A: having index of neighbours makes it easier to count neighbours later
	      //this needs to be done to kill danglers
	      for(i=0; i<tot_sites; i++)
	      {
		neighindex = 0;
		//only check for neighbours within a sensible range of current site index
		startloop = max(0, i - 2*length -4);
		endloop = mymin(tot_sites, i + 2*length +4);

		if(siteinfo[i][0] == 0)
		{
		  tempint=0;
		  for(j=startloop ; j< endloop ;j++)
		  {
		    if(siteinfo[j][0] == 0)
		    {
		      if(i!=j)
		      {
			smalldist = sqrt( pow(site_coords[i][0] - site_coords[j][0], 2) +  pow(site_coords[i][1] - site_coords[j][1], 2));
			if (smalldist < 0.6)
			{
			  tempint++;
			  neighbours[i][neighindex] = j;
			  neighindex++;
			}
		      }
		    }
		  }
		  sites[i] = tempint;
		  orig_sites[i] = tempint;
		}
	      }
	      
	      double rescaler, rescaler2; // these are essentially additional x (1-x) type envelopes to make things smoother
		  
	//determine if the atom is removed by determining the min and max y values allowed at it's x-coordinate
	      for(i=2*length*buffer_rows; i<tot_sites - 2*length*buffer_rows; i++)
	      {
		  thisminy=bottomy;
		  thismaxy=topy;
		  
		  for(l=0; l<sruns; l++)
		  {
		    if(site_coords[i][0] > startx[l] && site_coords[i][0] < endx[l])
		    {
			rescaler = ((site_coords[i][0] - startx[l]) * (endx[l] - site_coords[i][0]) / pow( (endx[l] -startx[l])/2, 2));
		    }
		    else
			rescaler = 0.0;
		    
		    if(site_coords[i][0] > startx[sruns+l] && site_coords[i][0] < endx[sruns+l])
		    {
			rescaler2 = ((site_coords[i][0] - startx[sruns+l]) * (endx[sruns+l] - site_coords[i][0]) / pow((endx[sruns+l] -startx[sruns+l])/2, 2));
		    }
		    else
			rescaler2 = 0.0;
		    
		    thisminy = thisminy + svals[l]*rescaler * pow(sin(M_PI * site_coords[i][0]/periods[l] + phases[l]), 2);
		    thismaxy = thismaxy - svals[sruns+l]*rescaler2 * pow(sin(M_PI * site_coords[i][0]/periods[sruns+l] + phases[sruns+l]), 2);
		  }
		  
		  if(site_coords[i][1] < thisminy || site_coords[i][1] > thismaxy)
		  {
		    siteinfo[i][0] = 1;
		  }
		
	      } 
	      
	//new neighbour count
		for(i=0; i<tot_sites; i++)
		{
		  if(siteinfo[i][0] == 0)
		  {
		    tempint =0;
		    
		    for(j=0; j<orig_sites[i]; j++)
		    {
		      if(siteinfo[neighbours[i][j]][0] ==0)
		      {
			tempint++;
		      }
		    }
		    sites[i] = tempint;
		  }
		}
		
		
		
	//vacancy disorder
	for(l=0; l<vacruns; l++)
	{
	  
	  //remove edge atoms with a certain probability
	  for(i=2*length*buffer_rows; i<tot_sites - 2*length*buffer_rows; i++)
	  {
	    if(sites[i] < 3 && siteinfo[i][0] == 0)
	    {
	      temprandnum = myRandNum(0.0, 1.0);
	      if(temprandnum < vacprob * (1.0 - l*1.0/vacruns))
	      {
		siteinfo[i][0] = 1;

	      }
	    }
	   }
	    
	  
	  //determine new edge atoms
	  for(i=0; i<tot_sites; i++)
	  {
	    if(siteinfo[i][0] == 0)
	    {
	      tempint =0;
	      
	      for(j=0; j<orig_sites[i]; j++)
	      {
		if(siteinfo[neighbours[i][j]][0] ==0)
		{
		  tempint++;
		}
	      }
	      sites[i] = tempint;
	    }
	  }
	  
	  //kill danglers
	  for(i=2*length*(buffer_rows); i<tot_sites - 2*length*(buffer_rows); i++)
	  {
	    if(sites[i] < 2)
	    {
	      siteinfo[i][0] = 1;
	    }
	  }
	  
	  //determine new edge atoms
	  for(i=0; i<tot_sites; i++)
	  {
	    if(siteinfo[i][0] == 0)
	    {
	      tempint =0;
	      
	      for(j=0; j<orig_sites[i]; j++)
	      {
		if(siteinfo[neighbours[i][j]][0] ==0)
		{
		  tempint++;
		}
	      }
	      sites[i] = tempint;
	    }
	  }
	  
	  
	    
	    
	}
	
	
	//one final killing of danglers, this time including inner chains of buffer region (if buffer_rows > 1)
	if(buffer_rows>1)
	{
	  for(i=2*length*(buffer_rows-1); i<tot_sites - 2*length*(buffer_rows-1); i++)
	  {
	    if(sites[i] < 2)
	    {
	      siteinfo[i][0] = 1;
	    }
	  }
	}
	
	
	//Memory freeing
	  free(neighbours[0]);
	  free(neighbours);
	  free(sites);
	  free(orig_sites);
	  free(svals); 
	  free(periods);
	  free(startx);
	  free(endx);
	  free(phases);
	  
	  
	  //some mapping / statistics
	      double topedgemax=0, topedgemin=5*length, bottomedgemax=0, bottomedgemin=5*length, widthmax=0, widthmin=5*length, widthavg=0;
	      double this_width, this_top, this_bottom;
	      
	      for(i=buffer_rows; i<length2-buffer_rows; i++)
	      {
		j=0;
		while(siteinfo[i*2*length +j][0] != 0)
		{
		  j++;
		}
		this_bottom = site_coords[i*2*length +j][1];
		
		j=0;
		while(siteinfo[(i+1)*2*length-1-j][0] != 0)
		{
		  j++;
		}
		this_top = site_coords[(i+1)*2*length-1-j][1];
		
		this_width = this_top - this_bottom;
		
		topedgemax = dmax(topedgemax, this_top);
		topedgemin = dmin(topedgemin, this_top);
		
		bottomedgemax = dmax(bottomedgemax, this_bottom);
		bottomedgemin = dmin(bottomedgemin, this_bottom);
		
		widthmin = dmin(widthmin, this_width);
		widthmax = dmax(widthmax, this_width);
		
		widthavg += this_width/(length2-2*buffer_rows);


		
		//printf("%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", site_coords[i*2*length +j][0], this_bottom, this_top, topedgemin,topedgemax, bottomedgemin, bottomedgemax );
	      }
		printf("# avg_width %lf\n# min_width %lf\n# max_width %lf\n# bottom_edge_fluc %lf\n# top_edge_fluc %lf\n", widthavg, widthmin, widthmax, bottomedgemax-bottomedgemin, topedgemax-topedgemin);

	 
		
		
		if(struc_out == 1)
		{
		  fprintf(out, "# avg_width %lf\n# min_width %lf\n# max_width %lf\n# bottom_edge_fluc %lf\n# top_edge_fluc %lf\n", widthavg, widthmin, widthmax, bottomedgemax-bottomedgemin, topedgemax-topedgemin);

		  for(l=0; l<2*length*length2; l++)
		  {
		    if(siteinfo[l][0] == 0)
		    {  
		      fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
		    }
		  }
		  fprintf(out, "\n");
		}
		
		
		  (SiteArray->pos) = site_coords;
		  (SiteArray->site_pots) = site_pots;
		  (SiteArray->siteinfo) = siteinfo;
                  (SiteArray->pert_pos) = NULL;
		  
		  int tempint2;
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
			  
				    

}







void exportRectConf(RectRedux *System, char *filename)
{
  //printf("%s\n", filename);
  FILE *fileout;
  char fullname[300];
  int length = System->length;
  int length2 = System->length2;
  int geo = System->geo;
  int i, j, k;
  int tot_sites = *(System->Ntot);
  

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
      fprintf(fileout, "%lf %lf %lf\n", (System->pos)[j][0], (System->pos)[j][1], (System->pos)[j][2]);
    }
  fclose(fileout);
  
  
  sprintf(fullname, "%s.site_pots", filename);
  fileout = fopen(fullname, "w");
    for(j=0; j<tot_sites; j++)
    {
      fprintf(fileout, "%.10lf\n", (System->site_pots)[j]);
    }
  fclose(fileout);
  
  
  
  sprintf(fullname, "%s.info", filename);
  fileout = fopen(fullname, "w");
       fprintf(fileout, "%d\n", *(System->Nrem));
       fprintf(fileout, "%d\n", *(System->Ntot));

  fclose(fileout);
  
  if( System->pert_pos != NULL )
  {
      
    sprintf(fullname, "%s.pert_pos", filename);
    fileout = fopen(fullname, "w");
    for(j=0; j<tot_sites; j++)
    {
        fprintf(fileout, "%lf	%lf   %lf\n", (System->pert_pos)[j][0], (System->pert_pos)[j][1], (System->pert_pos)[j][2]);
    }
    fclose(fileout);
  }
  
  

}


void importRectConf(RectRedux *System, int length, int length2, char *filename)
{
  //printf("%s	%d\n", filename, numcells);
  FILE *fileout;
  char fullname[300];
 // int length = System->length;
 // int length2 = System->length2;
 // int geo = System->geo;
  int i, j, k, temp;
  int tot_sites ;
  
 sprintf(fullname, "%s.info", filename);
  fileout = fopen(fullname, "r");
    
       fscanf(fileout, "%d", (System->Nrem));
       fscanf(fileout, "%d", (System->Ntot));

  fclose(fileout);
  tot_sites = *(System->Ntot);
  
	  System->pos = createNonSquareDoubleMatrix(tot_sites, 3);
 	  System->site_pots = createDoubleArray(tot_sites);
	  System->siteinfo = createNonSquareIntMatrix(tot_sites, 2);
          System->pert_pos = NULL;
// 	  (System->chaininfo) = createNonSquareIntMatrix(length2, 4);

	  
  
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
	fscanf(fileout, "%lf	%lf   %lf", &(System->pos)[j][0], &(System->pos)[j][1], &(System->pos)[j][2] );
      }
    fclose(fileout);
    
    
    sprintf(fullname, "%s.site_pots", filename);
    fileout = fopen(fullname, "r");
   
      for(j=0; j<tot_sites; j++)
      {
	fscanf(fileout, "%lf", &(System->site_pots)[j]);
      }
   
    fclose(fileout);
    
    //pertpos for strained systems
    
    sprintf(fullname, "%s.pert_pos", filename);
      if((fileout = fopen(fullname, "r") ))    
        {
            System->pert_pos = createNonSquareDoubleMatrix(tot_sites, 3);
            for(j=0; j<tot_sites; j++)
            {
                fscanf(fileout, "%lf	%lf   %lf", &(System->pert_pos)[j][0], &(System->pert_pos)[j][1], &(System->pert_pos)[j][2]);
            }
            fclose(fileout);
            printf("#strained geometry imported!\n");
        }
        else
            printf("#strained geometry NOT imported!\n");
    
    
    
    
//     sprintf(fullname, "%s.chaininfo", filename);
//     fileout = fopen(fullname, "r");
//     
//    
//       for(j=0; j<length2; j++)
//       {
// 	fscanf(fileout, "%d	%d	%d	%d", &(System->chaininfo)[j][0], &(System->chaininfo)[j][1], &(System->chaininfo)[j][2], &(System->chaininfo)[j][3]);
//       }
//     fclose(fileout);
   
  

  


}

//positions hall-type probes evenly along edges
//calculates centre of leads, then edge position from width
//returns starting chain index of probe
//this_probe indexing starts at 0
int HallPositioning(int length2, int num_side_probes, int this_probe, int buffer_rows, int geo_renorm, int width)
{
	int hall_denom=2*num_side_probes;
	int hall_start = (int) (( ((length2-2*buffer_rows)/hall_denom) - (width/(2*geo_renorm)))) + buffer_rows;
	if((hall_start % 2) == 1)
		hall_start+=1;
	int hall_jump =  2 * ((int) ((length2-2*buffer_rows)/hall_denom));
	
	
// 	return 2 * (int) (( (hall_mid*(length2-2*buffer_rows)/hall_denom) - (width/(2*geo_renorm)))/2) + buffer_rows;
	
	return hall_start + this_probe*hall_jump;
		
}


void Patchify ( RectRedux *System, patch_para *ppara, cellDivision *cellinfo, int struc_out, char *filename)
{
    
    
    
    int geo = (System->geo);
    int length = (System->length);
    int length2 = (System->length2);
    int **siteinfo = (System->siteinfo);		
    double **pos = (System->pos);
    double *site_pots = (System->site_pots);
    int *Nrem = (System->Nrem);
    int *Ntot = (System->Ntot);
    
    int numpatches = (ppara->numpatches);
    int i, j, k;
    
    ppara->pNrem = createIntArray(numpatches);
    
    double a1[2], a2[2], origin[2];
    if(geo == 0)
    {
        a1[0] = 0.5; a1[1] = -sqrt(3)/2;
        a2[0] = -0.5; a2[1] = -sqrt(3)/2;
        origin[0] = 0.5; origin[1] = 1/(2*sqrt(3));
    }
    
    if(geo == 1)
    {
        a1[1] = 0.5; a1[0] = sqrt(3)/2;
        a2[1] = -0.5; a2[0] = sqrt(3)/2;
        origin[1] = 0.5; origin[0] = 1/(2*sqrt(3));
    }
    
    FILE *out;
	  
    if(struc_out != 0)
    {
        out = fopen(filename, "w");
    }
    
    
    RectRedux *Patches[numpatches];
    for(i=0; i < numpatches; i++)
    {
        Patches[i] = (RectRedux *)malloc(sizeof(RectRedux));
        (Patches[i]->geo) = geo;
        (Patches[i]->length) = ppara->pl1[i];
        (Patches[i]->length2) = ppara->pl2[i];
        (Patches[i]->Nrem) = (int *)malloc(sizeof(int));
        (Patches[i]->Ntot) = (int *)malloc(sizeof(int));
        simpleRibbonGeo (Patches[i], NULL, 0, NULL);
        
        
//         //shift probe positions
        for(j=0; j< *(Patches[i]->Ntot); j++)
        {
            (Patches[i]->pos)[j][0] += ppara->pcoords[i][0] * a1[0] + ppara->pcoords[i][1] * a2[0] ;
            (Patches[i]->pos)[j][1] += ppara->pcoords[i][0] * a1[1] + ppara->pcoords[i][1] * a2[1] ;
            //printf("%lf	%lf\n", (Patches[i]->pos)[j][0], (Patches[i]->pos)[j][1]);
        }
        
        
        ppara->pNrem[i] = *(Patches[i]->Nrem);
        
        printf("#Patch %d, %d sites\n", i, ppara->pNrem[i]);
        
        
    }
    
    //combine into single device
    
            int *cell_dims = createIntArray(numpatches);
            int tot_new_dim=0;
            int rem_new_dim=0;
            
            for(i=0; i<numpatches; i++)
            {
                tot_new_dim += (*(Patches[i]->Ntot));
                rem_new_dim += (*(Patches[i]->Nrem));
            }
            int Nremnew = *Nrem + tot_new_dim;
            int Ntotnew = *Ntot + tot_new_dim;
            
            int **newsiteinfo = createNonSquareIntMatrix(Ntotnew, 2);
            double **newpos = createNonSquareDoubleMatrix(Ntotnew, 3);
            double *newpots = createDoubleArray(Ntotnew);
            double **newpertpos;
            
            if( System->pert_pos != NULL )
            {
                newpertpos = createNonSquareDoubleMatrix(Ntotnew, 3);
            }
            
            int tempcount=0;
            
            for(i=0; i<*Ntot; i++)
            {
                newsiteinfo[i][0] = siteinfo[i][0];
                newsiteinfo[i][1] = siteinfo[i][1];
                newpos[i][0] = pos[i][0];
                newpos[i][1] = pos[i][1];
                newpos[i][2] = pos[i][2];
                newpots[i] = site_pots[i];
                
                if( System->pert_pos != NULL )
                {
                    newpertpos[i][0] = (System->pert_pos)[i][0];
                    newpertpos[i][1] = (System->pert_pos)[i][1];
                    newpertpos[i][2] = (System->pert_pos)[i][2];
                }
            }
            tempcount=*Ntot;
            
            
            for(i=0; i<numpatches; i++)
            {
                for(k=0; k<*(Patches[i]->Ntot); k++)
                {
                    newsiteinfo[tempcount][0] = (Patches[i]->siteinfo)[k][0];
                    newsiteinfo[tempcount][1] = (Patches[i]->siteinfo)[k][1];
                    
                    newpos[tempcount][0] = (Patches[i]->pos)[k][0];
                    newpos[tempcount][1] = (Patches[i]->pos)[k][1];
                    newpos[tempcount][2] = (Patches[i]->pos)[k][2];
                    
                    newpots[tempcount] = (Patches[i]->site_pots)[k];
                    
                    if( Patches[i]->pert_pos != NULL )
                    {
                        newpertpos[tempcount][0] = (Patches[i]->pert_pos)[k][0];
                        newpertpos[tempcount][1] = (Patches[i]->pert_pos)[k][1];
                        newpertpos[tempcount][2] = (Patches[i]->pert_pos)[k][2];
                    }
                    
                    if(Patches[i]->pert_pos == NULL && System->pert_pos != NULL)
                    {
                        //printf("#null patches boo hoo!\n");
                        newpertpos[tempcount][0] = (Patches[i]->pos)[k][0];
                        newpertpos[tempcount][1] = (Patches[i]->pos)[k][1];
                        newpertpos[tempcount][2] = 0.0;
                        
                        //printf("#pp %lf %lf %lf %lf %lf %lf\n", newpos[tempcount][0], newpos[tempcount][1], newpos[tempcount][2], newpertpos[tempcount][0], newpertpos[tempcount][1], newpertpos[tempcount][2]);
                    }

                    tempcount++;
                }
            }
            
            
  
            
            
            free(siteinfo[0]); free(siteinfo);
            free(pos[0]); free(pos);
            free(site_pots);
            
            *Nrem = Nremnew;
            *Ntot = Ntotnew;
            
            (System->pos) = newpos;
            (System->site_pots) = newpots;
            (System->siteinfo) = newsiteinfo;
            if( System->pert_pos != NULL )
            {
                free( System->pert_pos[0]); free( System->pert_pos);
                (System->pert_pos ) = newpertpos;
            }
    
    //find "edge" sites to add to "group" - assumes perfect rectangular ribbon patch geometries
    //find "boundary" sites (e.g. edge of frame) site locations and save
 
    int edge_count=0;
    int boundary_count=0;
    int temp2=0, temp3=0, bigcount=0;
    
    
    //this is based on NNTB patches, with well-behaved and complete edges
    //be careful when generating devices, or this might not work
    
        (ppara->boundary)= (RectRedux *)malloc(sizeof(RectRedux));
        ((ppara->boundary)->Nrem) = (int *)malloc(sizeof(int));
        ((ppara->boundary)->Ntot) = (int *)malloc(sizeof(int));
        ((ppara->boundary)->pert_pos) = NULL;

        
//         (Patches[i]->geo) = geo;
//         (Patches[i]->length) = ppara->pl1[i];
//         (Patches[i]->length2) = ppara->pl2[i];
//         (Patches[i]->Nrem) = (int *)malloc(sizeof(int));
//         (Patches[i]->Ntot) = (int *)malloc(sizeof(int));
    
    
    if (geo == 0)
    {
        
        //site counting!
        edge_count = length2*2 + (length *2) -2;
        boundary_count = 2*(length + length2);
        for(i=0; i< numpatches; i++)
        {
            edge_count += (Patches[i]->length2)*2 + (Patches[i]->length *2) -2;
            boundary_count += 2*( (Patches[i]->length) + (Patches[i]->length2));
        }
//         printf("#edgecount %d\n", edge_count);
        
        //temp2 for edge counting, temp3 for boundary counting, tempcount for device counting
        //edge sites index; boundary sites given by coordinates
        tempcount=0, temp2=0, temp3=0, bigcount=0;
        //ppara->num_boundary_sites = boundary_count;
        *((ppara->boundary)->Ntot) = boundary_count;
        //ppara->boundary_pos = createNonSquareDoubleMatrix(boundary_count, 3);
        ((ppara->boundary)->pos) = createNonSquareDoubleMatrix(boundary_count, 3);
        //ppara->boundary_subl = createIntArray(boundary_count);
        ((ppara->boundary)->siteinfo) = createNonSquareIntMatrix(boundary_count, 2);
        
        cellinfo->group_dim = edge_count;
        cellinfo->group_sites= createIntArray(edge_count);
        
        //Main device edge and boundary sites
        
            //bottom left (always an edge!)
            tempcount=bigcount + 0;
            (cellinfo->group_sites)[temp2] = tempcount;
            ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
            ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/sqrt(3));
            ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
            temp2++; temp3++;
            
            //left-most chain  (the -2 in the loop takes care of odd-integer lengths)
            for(j=0; j<2*length-2; j+=4)
            {
                tempcount = bigcount +j+1;
                (cellinfo->group_sites)[temp2] = tempcount;
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] - 0.5;
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/(2*sqrt(3)));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp2++; temp3++;
                
                tempcount = bigcount +j+2;
                (cellinfo->group_sites)[temp2] = tempcount;
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] - 0.5;
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/(2*sqrt(3)));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp2++; temp3++;
            }
            
            //top left - different edge for odd and even lengths
            tempcount = bigcount + 2*length -1;
            (cellinfo->group_sites)[temp2] = tempcount;
            temp2++;
           
            ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
            ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/sqrt(3));
            ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
            temp3++;
            if( (length % 2) == 1)
            {
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] - 0.5;;
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/(2*sqrt(3)));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp3++;
            }
            
            
            
            //loop over intermediate chains)
            for(k=1; k<length2-1; k++)
            {
                tempcount = bigcount + 2*length*(k) ;
                (cellinfo->group_sites)[temp2] = tempcount;
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/sqrt(3));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp2++; temp3++;
                
                tempcount = bigcount + 2*length*(k+1) - 1 ;
                (cellinfo->group_sites)[temp2] = tempcount;
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/sqrt(3));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp2++; temp3++;
            }
            
            
            //bottom right - edge always connects to two boundary sites
            tempcount = bigcount + 2*length*(length2-1) ;
            (cellinfo->group_sites)[temp2] = tempcount;
            temp2++; 
            ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
            ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/sqrt(3));
            ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
            temp3++;
            ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] + 0.5;
            ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/(2*sqrt(3)));
            ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
            temp3++;
            
            //right-most chain
            for(j=3; j<2*length-1; j+=4)
            {
                tempcount = bigcount + 2*length*(length2-1) +j;
                (cellinfo->group_sites)[temp2] = tempcount;
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] + 0.5;
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/(2*sqrt(3)));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp2++; temp3++;
                
                tempcount = bigcount + 2*length*(length2-1) +j+1;
                (cellinfo->group_sites)[temp2] = tempcount;
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] + 0.5;
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/(2*sqrt(3)));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp2++; temp3++;
            }
            
            //top right -- different for odd and even
            tempcount = bigcount + 2*length*length2-1 ;
            (cellinfo->group_sites)[temp2] = tempcount;
            temp2++; 
            if( (length % 2) == 0)
            {            
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] + 0.5;
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/(2*sqrt(3)));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp3++;
            }
            ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
            ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/sqrt(3));
            ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
            temp3++;
            
            
            
            bigcount = 2*length*length2;
        //Boundaries and edges of the other patches
            for(i=0; i<numpatches; i++)
            {
                
                
                //bottom left (always an edge!)
                tempcount=bigcount + 0;
                (cellinfo->group_sites)[temp2] = tempcount;
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/sqrt(3));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp2++; temp3++;
                
                //left-most chain  (the -2 in the loop takes care of odd-integer lengths)
                for(j=0; j<2*(Patches[i]->length)-2; j+=4)
                {
                    tempcount = bigcount +j+1;
                    (cellinfo->group_sites)[temp2] = tempcount;
                    ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] - 0.5;
                    ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/(2*sqrt(3)));
                    ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                    temp2++; temp3++;
                    
                    tempcount = bigcount +j+2;
                    (cellinfo->group_sites)[temp2] = tempcount;
                    ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] - 0.5;
                    ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/(2*sqrt(3)));
                    ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                    temp2++; temp3++;
                }
                
                //top left - different edge for odd and even lengths
                tempcount = bigcount + 2*(Patches[i]->length) -1;
                //if( ((Patches[i]->length) % 2) == 0)
                //{
                    (cellinfo->group_sites)[temp2] = tempcount;
                    temp2++;
                //}
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/sqrt(3));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp3++;
                
                if( ((Patches[i]->length) % 2) == 1)
                {
                    ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] - 0.5;;
                    ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/(2*sqrt(3)));
                    ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                    temp3++;
                }
            
                
                
                
                
                //loop over intermediate chains)
                for(k=1; k<(Patches[i]->length2)-1; k++)
                {
                    tempcount = bigcount + 2*(Patches[i]->length)*(k) ;
                    (cellinfo->group_sites)[temp2] = tempcount;
                    ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
                    ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/sqrt(3));
                    ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                    temp2++; temp3++;
                    
                    tempcount = bigcount + 2*(Patches[i]->length)*(k+1) - 1 ;
                    (cellinfo->group_sites)[temp2] = tempcount;
                    ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
                    ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/sqrt(3));
                    ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                    temp2++; temp3++;
                }
                
                
                //bottom right - edge always connects to two boundary sites
                tempcount = bigcount + 2*(Patches[i]->length)*((Patches[i]->length2)  -1) ;
                (cellinfo->group_sites)[temp2] = tempcount;
                temp2++; 
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/sqrt(3));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp3++;
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] + 0.5;
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/(2*sqrt(3)));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp3++;
                
                //right-most chain
                for(j=3; j<2*(Patches[i]->length)-1; j+=4)
                {
                    tempcount = bigcount + 2*(Patches[i]->length)*( (Patches[i]->length2)-1) +j;
                    (cellinfo->group_sites)[temp2] = tempcount;
                    ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] + 0.5;
                    ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/(2*sqrt(3)));
                    ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                    temp2++; temp3++;
                    
                    tempcount = bigcount + 2*(Patches[i]->length)*((Patches[i]->length2)-1) +j+1;
                    (cellinfo->group_sites)[temp2] = tempcount;
                    ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] + 0.5;
                    ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/(2*sqrt(3)));
                    ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                    temp2++; temp3++;
                }
                
                //top right -- different for odd and even
                tempcount = bigcount + 2*(Patches[i]->length)*(Patches[i]->length2)-1 ;
                (cellinfo->group_sites)[temp2] = tempcount;
                temp2++; 
                if( ((Patches[i]->length) % 2) == 0)
                {            
                    ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0] + 0.5;
                    ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] - (1.0/(2*sqrt(3)));
                    ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                    temp3++;
                }
                ((ppara->boundary)->pos)[temp3][0] = (System->pos)[tempcount][0];
                ((ppara->boundary)->pos)[temp3][1] = (System->pos)[tempcount][1] + (1.0/sqrt(3));
                ((ppara->boundary)->siteinfo)[temp3][1] = 1 - (System->siteinfo)[tempcount][1];
                temp3++;
            
                
                
                bigcount += 2*(Patches[i]->length)*(Patches[i]->length2);
                
                
                
            }
            (cellinfo->group_dim) = temp2;

            
    }
//                     printf("#edgecount %d\n", cellinfo->group_dim);

    
    // calculate number of unique GF matrix elements required
    int max_conn = (boundary_count + edge_count)*boundary_count; 
    int **max_sep_indices = createNonSquareIntMatrix(max_conn, 4); //m, n, diagt, num_instances
    //indexes the connections between boundary, device sites in terms of unique GF indices
    int **boundary_device_mat = createNonSquareIntMatrix(boundary_count, boundary_count+edge_count);
    int conn_count=0; //counts the unique GF elements so far
    int already_calc=0;
    
    double **testpos = createNonSquareDoubleMatrix(2, 3);
    int testout[3];
    

    for(i=0; i<boundary_count; i++)
    {
        testpos[0][0] = ((ppara->boundary)->pos)[i][0];
        testpos[0][1] = ((ppara->boundary)->pos)[i][1];
        testpos[0][2] = ((ppara->boundary)->siteinfo)[i][1] * 1.0;
        
        for(j=0; j< boundary_count; j++)
        {
            testpos[1][0] = ((ppara->boundary)->pos)[j][0];
            testpos[1][1] = ((ppara->boundary)->pos)[j][1];
            testpos[1][2] = ((ppara->boundary)->siteinfo)[j][1] * 1.0;
            
            simplifyIndices(testpos, testout, origin, a1, a2 );
            
//             if(testout[0]==1 && testout[1] == 1 && testout[2] ==2)
//             {
//                 printf("#%d, %d SUBL check %lf, %lf\n", i, j, testpos[0][1], testpos[1][1]);
//             }
            
//             if(i==0 && j==10)
//             {
//                 printf("%d, %d, %d\n", testout[0], testout[1], testout[2]);
//             }
            
            already_calc=0;
            for(k=0; k<conn_count; k++)
            {
                //print("start\n");
                if(testout[0] == max_sep_indices[k][0]
                    && testout[1] == max_sep_indices[k][1]
                    && testout[2] == max_sep_indices[k][2]) 
                    {
                        already_calc = 1;
                        boundary_device_mat[i][j] = k;
                        max_sep_indices[k][3] ++;
                    }
            }
            if (already_calc == 0)
            {
                max_sep_indices[conn_count][0] = testout[0];
                max_sep_indices[conn_count][1] = testout[1];
                max_sep_indices[conn_count][2] = testout[2];
                max_sep_indices[conn_count][3] = 1;
                boundary_device_mat[i][j] = conn_count;
                conn_count++;
            }
            
//             if(i==0 && j==10)
//             {
//                 printf("bdm: %d\n", boundary_device_mat[i][j] );
//                 
//                 printf("%d, %d, %d\n", max_sep_indices[10][0], max_sep_indices[10][1], max_sep_indices[10][2]);
//                 
//             }
            
                
//             if(i==10 && j==0)
//             {
//                 printf("#10,0 test: %d  %d  %d\n", testout[0], testout[1], testout[2]);
//             }
        }
        
        for(j=0; j< edge_count; j++)
        {
            testpos[1][0] = newpos[(cellinfo->group_sites)[j]][0];
            testpos[1][1] = newpos[(cellinfo->group_sites)[j]][1];
            testpos[1][2] = newsiteinfo[(cellinfo->group_sites)[j]][1]*1.0;
            
            simplifyIndices(testpos, testout, origin, a1, a2 );
            
            already_calc=0;
            for(k=0; k<conn_count; k++)
            {
                //print("start\n");
                if(testout[0] == max_sep_indices[k][0]
                    && testout[1] == max_sep_indices[k][1]
                    && testout[2] == max_sep_indices[k][2]) 
                    {
                        already_calc = 1;
                        boundary_device_mat[i][boundary_count+j] = k;
                        max_sep_indices[k][3] ++;
                    }
            }
            if (already_calc == 0)
            {
                max_sep_indices[conn_count][0] = testout[0];
                max_sep_indices[conn_count][1] = testout[1];
                max_sep_indices[conn_count][2] = testout[2];
                max_sep_indices[conn_count][3] = 1;
                boundary_device_mat[i][boundary_count+j] = conn_count;
                conn_count++;
            }
            
        }
    }
    
//     printf("bdm: %d\n", boundary_device_mat[0][10] );
//                 
//     printf("%d, %d, %d\n", max_sep_indices[10][0], max_sep_indices[10][1], max_sep_indices[10][2]);
    

    double dx2, dy2;
    double ABvec[2];
    ABvec[0] = (a1[0] + a2[0])/3;
    ABvec[1] = (a1[1] + a2[1])/3;
    printf("# unique GFs: %d\n", conn_count);
    for(k=0; k<conn_count; k++)
    {
        dx2 = max_sep_indices[k][0]*a1[0] + max_sep_indices[k][1]*a2[0] - (1.5*max_sep_indices[k][2]*max_sep_indices[k][2] - 2.5*max_sep_indices[k][2])*ABvec[0];
        dy2 = max_sep_indices[k][0]*a1[1] + max_sep_indices[k][1]*a2[1] - (1.5*max_sep_indices[k][2]*max_sep_indices[k][2] - 2.5*max_sep_indices[k][2])*ABvec[1];
    
        
        
//         printf("# (%d, %d, %d) ## %d instances %lf sep\n", max_sep_indices[k][0], max_sep_indices[k][1], max_sep_indices[k][2], max_sep_indices[k][3], sqrt(dx2*dx2+dy2*dy2));
    }
    
    
    (ppara->sep_indices) = createNonSquareIntMatrix(max_conn, 3); //m, n, diagt,
    
    for(i=0; i<conn_count; i++)
    {
        (ppara->sep_indices)[i][0] = max_sep_indices[i][0];
        (ppara->sep_indices)[i][1] = max_sep_indices[i][1];
        (ppara->sep_indices)[i][2] = max_sep_indices[i][2];
    }
    free(max_sep_indices[0]); free(max_sep_indices);
    free(testpos[0]); free(testpos);
    (ppara->conn_count) = conn_count;
    (ppara->boundary_device) = boundary_device_mat;
    
    
    
    	if(struc_out != 0)
	{
	  for(i=0; i<Ntotnew; i++)
	  {
	    if(newsiteinfo[i][0] == 0 && newsiteinfo[i][1] == 0)
	    {
	      fprintf(out, "%lf	%lf\n", newpos[i][0], newpos[i][1]);
	    }
	  }
	  fprintf(out, "\n");
	  for(i=0; i<Ntotnew; i++)
	  {
	    if(newsiteinfo[i][0] == 0 && newsiteinfo[i][1] == 1)
	    {
	      fprintf(out, "%lf	%lf\n", newpos[i][0], newpos[i][1]);
	    }
	  }
	  fprintf(out, "\n");
	  
          for(i=0; i<cellinfo->group_dim; i++)
	  {
            if(newsiteinfo[(cellinfo->group_sites)[i]][0] == 0 && newsiteinfo[(cellinfo->group_sites)[i]][1] == 0)
	    {
	      fprintf(out, "%lf	%lf\n", newpos[(cellinfo->group_sites)[i]][0], newpos[(cellinfo->group_sites)[i]][1]);
	    }
	  }
	  fprintf(out, "\n");
          for(i=0; i<cellinfo->group_dim; i++)
	  {
            if(newsiteinfo[(cellinfo->group_sites)[i]][0] == 0 && newsiteinfo[(cellinfo->group_sites)[i]][1] == 1)
	    {
	      fprintf(out, "%lf	%lf\n", newpos[(cellinfo->group_sites)[i]][0], newpos[(cellinfo->group_sites)[i]][1]);
	    }
	  }
              
         fprintf(out, "\n");   
         for(i=0; i<*((ppara->boundary)->Ntot); i++)
         {
             if( ((ppara->boundary)->siteinfo)[i][1] == 0)
             {
                 fprintf(out, "%lf	%lf\n", ((ppara->boundary)->pos)[i][0], ((ppara->boundary)->pos)[i][1]  );
             }
         }
         fprintf(out, "\n");   
         for(i=0; i<*((ppara->boundary)->Ntot); i++)
         {
             if( ((ppara->boundary)->siteinfo)[i][1] == 1)
             {
                 fprintf(out, "%lf	%lf\n", ((ppara->boundary)->pos)[i][0], ((ppara->boundary)->pos)[i][1]  );
             }
         }     

	  
	  
	}

        if(struc_out != 0)
        {
            fclose(out);
        }
        //exit(1);

}


//Returns the simplest unique indices for lattice separations in graphene
void simplifyIndices(double **input, int *output, double *origin, double *a1, double *a2 )
{
    int m, n, diagt;
    int m1, n1;
    
    int a = (int)round(input[0][2]), b = (int)round(input[1][2]);
    int diagt_in = b*(3*b -3*a -2) +2*a;   //0 for AA,BB, 1 for AB, 2 for BA
//    printf("#SUBs: %d, %d, %d\n", a, b, diagt_in);

    //separation vector between A and B atoms within sublattice
    double ABvec[2];
    ABvec[0] = (a1[0] + a2[0])/3;
    ABvec[1] = (a1[1] + a2[1])/3;
//     printf("#AB %lf %lf\n", ABvec[0], ABvec[1]);
    
    double delxi = input[1][0] - input[0][0], delyi = input[1][1] - input[0][1];
    double alpha1[2], alpha2[2], angle=0;
    double delx = delxi, dely=delyi;
    double dx2, dy2;
    int dorig= diagt_in;
//     printf("#input:  %lf, %lf     %lf, %lf\n", input[0][0], input[0][1], input[1][0], input[1][1]);
//     printf("#dell: %lf, %lf\n", delx, dely);
    
    if(diagt_in ==1)
    {
        delx -= ABvec[0];
        dely -= ABvec[1];
    }
    
    if(diagt_in ==2)
    {
        delx += ABvec[0];
        dely += ABvec[1];
    }   
    
    
    
    m1 = (int) round ( (delx/a2[0] - dely/a2[1]) / (a1[0]/a2[0] - a1[1]/a2[1]) ) ;
    n1 = (int) round ( (delx/a1[0] - dely/a1[1]) / (a2[0]/a1[0] - a2[1]/a1[1]) ) ;
    diagt = diagt_in;
    
    m=m1; n=n1;
    //reduce and simplify!
    int loopcount=0;
    
//     if(m==-1 && n == 2 && dorig ==2)
//     {
//         printf("this iteration... %lf, %lf (%lf) \n", delxi, delyi, sqrt(delxi*delxi+delyi*delyi));
//     }
    
    double init_sep, final_sep;
    init_sep = sqrt(delxi*delxi +delyi*delyi);

    if( m1!=0 || n1!=0 || diagt !=0)
    {
//         printf("orig (%d, %d, %d) (%lf, %lf) %lf \n", m1 , n1, diagt, delxi, delyi, sqrt(delxi*delxi +delyi*delyi));
//         printf("#del %lf %lf\n", delx, dely);
        //diagt=2 -> diagt=1
        if (diagt_in == 2)
        {
            m1 = -m1; n1 = -n1; diagt=1;
            diagt_in = diagt;
//             printf("->inversion: (%d, %d, %d)\t", m1 , n1 , diagt );
            delx=-delx;
            dely=-dely;
        }
        
        //Bring everything to Z0
        while( (m1 <= 0 || n1 < 0 || diagt > 1) && loopcount < 4 )
        {
         
            if (diagt_in == 2)
            {
                m1 = -m1; n1 = -n1; diagt=1;
                diagt_in = diagt;
    //             printf("->inversion: (%d, %d, %d)\t", m1 , n1 , diagt );
                delx=-delx;
                dely=-dely;
            }
            
            //rotate from zone 1 back to zone 0
            if( (m1<=0) && (n1>0) && (abs(m1) < n1))
            {
                angle=-M_PI/3;
                alpha1[0] = a1[0] * cos(angle) - a1[1] * sin (angle);
                alpha1[1] = a1[0] * sin(angle) + a1[1] * cos (angle);
                alpha2[0] = a2[0] * cos(angle) - a2[1] * sin (angle);
                alpha2[1] = a2[0] * sin(angle) + a2[1] * cos (angle);
                ABvec[0] = (alpha1[0] + alpha2[0])/3;
                ABvec[1] = (alpha1[1] + alpha2[1])/3;
                
                
                //delx = delxi; 
                //dely = delyi;
                   
                 if(diagt_in ==1)
                {
                     delx += alpha1[0];
                     dely += alpha1[1];
                     diagt = 2;
                }
                if(diagt_in ==2)
                {
                    delx -= alpha1[0];
                    dely -= alpha1[0];
                    diagt = 1;
                }
                
                
                
                m1 = (int) round ( (delx/alpha2[0] - dely/alpha2[1]) / (alpha1[0]/alpha2[0] - alpha1[1]/alpha2[1]) ) ;
                n1 = (int) round ( (delx/alpha1[0] - dely/alpha1[1]) / (alpha2[0]/alpha1[0] - alpha2[1]/alpha1[1]) ) ;
//                 printf("->Z1\t (%d, %d, %d)\t", m1 , n1 , diagt );
                diagt_in = diagt;
                delx = m1*a1[0] + n1*a2[0] ;
                dely = m1*a1[1] + n1*a2[1] ;
            }
            
            //rotate from zone 2 back to zone 0
            if( (m1<0) && (n1>0) && (abs(m1) >= n1))
            {
                angle=-2*M_PI/3;
                alpha1[0] = a1[0] * cos(angle) - a1[1] * sin (angle);
                alpha1[1] = a1[0] * sin(angle) + a1[1] * cos (angle);
                alpha2[0] = a2[0] * cos(angle) - a2[1] * sin (angle);
                alpha2[1] = a2[0] * sin(angle) + a2[1] * cos (angle);
                ABvec[0] = (alpha1[0] + alpha2[0])/3;
                ABvec[1] = (alpha1[1] + alpha2[1])/3;
                
                                
                if(diagt_in ==1)
                {
                     delx -= alpha2[0];
                     dely -= alpha2[1];
                }
                if(diagt_in ==2)
                {
                    delx += alpha2[0];
                    dely += alpha2[0];
                }
                
                m1 = (int) round ( (delx/alpha2[0] - dely/alpha2[1]) / (alpha1[0]/alpha2[0] - alpha1[1]/alpha2[1]) ) ;
                n1 = (int) round ( (delx/alpha1[0] - dely/alpha1[1]) / (alpha2[0]/alpha1[0] - alpha2[1]/alpha1[1]) ) ;
//                 printf("->Z2\t (%d, %d, %d)\t", m1 , n1 , diagt );
                diagt_in = diagt;
                delx = m1*a1[0] + n1*a2[0] ;
                dely = m1*a1[1] + n1*a2[1] ;
//                 printf("delx: %lf, dely: %lf\n", delx, dely);

            }
            
            
            //rotate from zone 3 back to zone 0
            if( (m1<0) && (n1<=0))
            {

                m1 = -m1;
                n1=-n1;
                if(diagt_in==1)
                    diagt=2;
                if(diagt_in==2)
                    diagt=1;
                
//                 printf("->Z3\t (%d, %d, %d)\t", m1 , n1 , diagt );
                diagt_in = diagt;
                delx = m1*a1[0] + n1*a2[0] ;
                dely = m1*a1[1] + n1*a2[1] ;
            }
            
            //rotate from zone 4 back to zone 0
            if( (m1>=0) && (n1<0) && (abs(n1) > m1))
            {
                angle=-4*M_PI/3;
                alpha1[0] = a1[0] * cos(angle) - a1[1] * sin (angle);
                alpha1[1] = a1[0] * sin(angle) + a1[1] * cos (angle);
                alpha2[0] = a2[0] * cos(angle) - a2[1] * sin (angle);
                alpha2[1] = a2[0] * sin(angle) + a2[1] * cos (angle);
                ABvec[0] = (alpha1[0] + alpha2[0])/3;
                ABvec[1] = (alpha1[1] + alpha2[1])/3;
                
                                
                if(diagt_in ==1)
                {
                     delx -= alpha1[0];
                     dely -= alpha1[1];
                }
                if(diagt_in ==2)
                {
                    delx += alpha1[0];
                    dely += alpha1[0];
                }
                
                m1 = (int) round ( (delx/alpha2[0] - dely/alpha2[1]) / (alpha1[0]/alpha2[0] - alpha1[1]/alpha2[1]) ) ;
                n1 = (int) round ( (delx/alpha1[0] - dely/alpha1[1]) / (alpha2[0]/alpha1[0] - alpha2[1]/alpha1[1]) ) ;
//                 printf("->Z4\t (%d, %d, %d)\t", m1 , n1 , diagt );
                diagt_in = diagt;
                delx = m1*a1[0] + n1*a2[0] ;
                dely = m1*a1[1] + n1*a2[1] ;
            }
            
            
             //rotate from zone 5 back to zone 0
            if( (m1>0) && (n1<0) && (abs(n1) <= m1))
            {
                angle=-5*M_PI/3;
                alpha1[0] = a1[0] * cos(angle) - a1[1] * sin (angle);
                alpha1[1] = a1[0] * sin(angle) + a1[1] * cos (angle);
                alpha2[0] = a2[0] * cos(angle) - a2[1] * sin (angle);
                alpha2[1] = a2[0] * sin(angle) + a2[1] * cos (angle);
                ABvec[0] = (alpha1[0] + alpha2[0])/3;
                ABvec[1] = (alpha1[1] + alpha2[1])/3;
                
//                 printf("delx: %lf, dely: %lf\n", delx, dely);
                if(diagt_in ==1)
                {
                     diagt=2;
                     delx += alpha2[0];
                     dely += alpha2[1];
                }
                if(diagt_in ==2)
                {
                    diagt=1;
                    delx -= alpha2[0];
                    dely -= alpha2[0];
                }
//                 printf("delx: %lf, dely: %lf\n", delx, dely);

                
               /* 
                if(diagt_in ==1)
                {
                    diagt = 2;
                    delx += ABvec[0];
                    dely += ABvec[1];
                }
                if(diagt_in ==2)
                {
                    delx -= ABvec[0];
                    dely -= ABvec[1];
                    diagt = 1;
                }*/
                
                m1 = (int) round ( (delx/alpha2[0] - dely/alpha2[1]) / (alpha1[0]/alpha2[0] - alpha1[1]/alpha2[1]) ) ;
                n1 = (int) round ( (delx/alpha1[0] - dely/alpha1[1]) / (alpha2[0]/alpha1[0] - alpha2[1]/alpha1[1]) ) ;
//                 printf("->Z5\t (%d, %d, %d)\t", m1 , n1 , diagt );
                diagt_in = diagt;
                delx = m1*a1[0] + n1*a2[0] ;
                dely = m1*a1[1] + n1*a2[1] ;
// printf("delx: %lf, dely: %lf\n", delx, dely);

            }

//            printf("\n" );  
            loopcount++;
        }
        //printf("->final \t (%d, %d, %d)\n", m1 , n1 , diagt );
        
        


    }
    ABvec[0] = (a1[0] + a2[0])/3;
    ABvec[1] = (a1[1] + a2[1])/3;
    
    int swap;
    
     if(n1 > m1)
     {
         swap = n1;
         n1 = m1;
         m1=swap;
     }
    
    dx2 = m1*a1[0] + n1*a2[0] - (1.5*diagt*diagt - 2.5*diagt)*ABvec[0];
    dy2 = m1*a1[1] + n1*a2[1] - (1.5*diagt*diagt - 2.5*diagt)*ABvec[1];
    final_sep = sqrt(dx2*dx2 + dy2*dy2);
//     printf("\t \n", m1 , n1, diagt, init_sep, final_sep );

    //will print out if something is dodgy!
    if(init_sep - final_sep > 0.0001)
            printf("orig (%d, %d, %d)\t->final (%d, %d, %d)  %lf %lf\n", m , n, dorig, m1 , n1, diagt, init_sep, final_sep  );

     if (m1< 0 || n1 < 0)
     {
          printf("orig (%d, %d, %d)\t->final (%d, %d, %d)  %lf %lf\n", m , n, dorig, m1 , n1, diagt, init_sep, final_sep  );
     }
     
//      if(m==-1 && n == 2 && dorig ==2)
//     {
//         printf("orig (%d, %d, %d)\t->final (%d, %d, %d)  %lf %lf\n", m , n, dorig, m1 , n1, diagt, init_sep, final_sep  );
//     }

    
    //printf ("%d %d %d\n", m1, n1, diagt);
    //printf("## (%d %d %d) del %lf %lf    %lf %lf\n", m1, n1, diagt_in, input[1][0] - input[0][0], input[1][1] - input[0][1], (m1*a1[0] + n1*a2[0]), (m1*a1[1] + n1*a2[1]));
    
    output[0] = m1; output[1] = n1; output[2] = diagt;
    
//     exit(1);
    
}



void HallBarify (RectRedux *System, RectRedux **Leads, hallbpara *hall_para, lead_para *params, int struc_out, char *filename)
{
  int geo = (System->geo);
  int length = (System->length);
  int length2 = (System->length2);
  int **siteinfo = (System->siteinfo);		
  double **pos = (System->pos);
  double *site_pots = (System->site_pots);
//   int **chaininfo = (System->chaininfo);
  int *Nrem = (System->Nrem);
  int *Ntot = (System->Ntot);
  
  
  int ntop = hall_para->num_top_probes;
  int nbot = hall_para->num_bot_probes;
  int *toppx = hall_para->toppx;
  int *toppw = hall_para->toppw;
  int *toppc = hall_para->toppc;
  int *botpx = hall_para->botpx;
  int *botpw = hall_para->botpw;
  int *botpc = hall_para->botpc;
  
  
  int Nnew, i, j, k;
  int tbgeo;
  double ybot, ytop, y_cell_diff;
  double xleft, xright, x_cell_diff;
  
  FILE *out;
	  
  if(struc_out != 0)
  {
    out = fopen(filename, "w");
  }
  
  
  //will be ropey for odd integer ZGNRs
  if(geo == 0)
  {
    tbgeo=1;
    ybot = 1/(2*sqrt(3)) - sqrt(3);
    ytop = 1/(2*sqrt(3)) + length * sqrt(3)/2;
    xleft = 0.5;
    xright = length2*1.0 - 0.5;
    x_cell_diff = 1.0;
    y_cell_diff = sqrt(3);
  }
  
  if(geo == 1)
  {
    tbgeo=0;
    
    ybot = -1.0;
    ytop = 0.5*length;
    xleft = -1/(2*sqrt(3)) + sqrt(3)/2;
    xright = -1/(2*sqrt(3))+ sqrt(3)/2 + (length2-1)*sqrt(3);
    x_cell_diff = sqrt(3);
    y_cell_diff = 1.0;
    
  }
  
//   printf("#ntop %d, nbot %d, num y %d\n", ntop, nbot, toppw[0]);
  
      RectRedux *CellSecs[ntop + nbot];
      for(i=0; i < ntop; i++)
      {
	CellSecs[i] = (RectRedux *)malloc(sizeof(RectRedux));
	(CellSecs[i]->geo) = tbgeo;
	(CellSecs[i]->length) = toppw[i];
	(CellSecs[i]->length2) = 1;
	(CellSecs[i]->Nrem) = (int *)malloc(sizeof(int));
	(CellSecs[i]->Ntot) = (int *)malloc(sizeof(int));
	simpleRibbonGeo (CellSecs[i], NULL, 0, NULL);
	
	swapxy((CellSecs[i]->pos), *(CellSecs[i]->Ntot));
	
	
	//shift probe positions
	for(j=0; j< *(CellSecs[i]->Ntot); j++)
	{
	 (CellSecs[i]->pos)[j][0] += (xleft + toppx[i]*x_cell_diff);
	 (CellSecs[i]->pos)[j][1] += (ytop);
//   	 printf("%lf	%lf\n", (CellSecs[i]->pos)[j][0], (CellSecs[i]->pos)[j][1]);
	}
	
      }
      
      for(i=0; i < nbot; i++)
      {
	CellSecs[ntop+i] = (RectRedux *)malloc(sizeof(RectRedux));
	(CellSecs[ntop+i]->geo) = tbgeo;
	(CellSecs[ntop+i]->length) = botpw[i];
	(CellSecs[ntop+i]->length2) = 1;
	(CellSecs[ntop+i]->Nrem) = (int *)malloc(sizeof(int));
	(CellSecs[ntop+i]->Ntot) = (int *)malloc(sizeof(int));
	simpleRibbonGeo (CellSecs[ntop+i], NULL, 0, NULL);
	
	swapxy((CellSecs[ntop+i]->pos), *(CellSecs[ntop+i]->Ntot));

	
	//shift probe positions
	for(j=0; j< *(CellSecs[ntop+i]->Ntot); j++)
	{
	 (CellSecs[ntop+i]->pos)[j][0] += (xleft + botpx[i]*x_cell_diff);
	 (CellSecs[ntop+i]->pos)[j][1] += (ybot);
	 
//   	 printf("%lf	%lf\n", (CellSecs[ntop+i]->pos)[j][0], (CellSecs[ntop+i]->pos)[j][1]);
	}
	
      }
      
      int *cell_dims = createIntArray(ntop + nbot);
      int tot_new_dim=0;
      int rem_new_dim=0;
      for(i=0; i<ntop; i++)
      {
	cell_dims[i] = (*(CellSecs[i]->Ntot))*toppc[i];
	tot_new_dim += cell_dims[i];
	rem_new_dim += (*(CellSecs[i]->Nrem))*toppc[i];
      }
      for(i=0; i<nbot; i++)
      {
	cell_dims[ntop+i] = (*(CellSecs[ntop+i]->Ntot))*botpc[i];
	tot_new_dim += cell_dims[ntop+i];
	rem_new_dim += (*(CellSecs[ntop+i]->Nrem))*botpc[i];
      }
      
      
//       printf("#tot new sites: %d\n", tot_new_dim);
      int Nremnew = *Nrem + tot_new_dim;
      int Ntotnew = *Ntot + tot_new_dim;
      
      int **newsiteinfo = createNonSquareIntMatrix(Ntotnew, 2);
      double **newpos = createNonSquareDoubleMatrix(Ntotnew, 3);
      double *newpots = createDoubleArray(Ntotnew);

      int tempcount=0;
      
      for(i=0; i<*Ntot; i++)
      {
	newsiteinfo[i][0] = siteinfo[i][0];
	newsiteinfo[i][1] = siteinfo[i][1];
	newpos[i][0] = pos[i][0];
	newpos[i][1] = pos[i][1];
	newpos[i][2] = pos[i][2];
	newpots[i] = site_pots[i];
      }
      tempcount=*Ntot;
      for(i=0; i<ntop; i++)
      {
	for(j=0; j<toppc[i]; j++)
	{
	  for(k=0; k<*(CellSecs[i]->Ntot); k++)
	  {
	    newsiteinfo[tempcount][0] = (CellSecs[i]->siteinfo)[k][0];
	    newsiteinfo[tempcount][1] = (CellSecs[i]->siteinfo)[k][1];
	    
	    newpos[tempcount][0] = (CellSecs[i]->pos)[k][0];
	    newpos[tempcount][1] = (CellSecs[i]->pos)[k][1] + j *y_cell_diff;
	    newpos[tempcount][2] = (CellSecs[i]->pos)[k][2];
	    
	    newpots[tempcount] = (CellSecs[i]->site_pots)[k];

	    tempcount++;
	  }
	}
      }
      
      for(i=0; i<nbot; i++)
      {
	for(j=0; j<botpc[i]; j++)
	{
	  for(k=0; k<*(CellSecs[ntop+i]->Ntot); k++)
	  {
	    newsiteinfo[tempcount][0] = (CellSecs[ntop+i]->siteinfo)[k][0];
	    newsiteinfo[tempcount][1] = (CellSecs[ntop+i]->siteinfo)[k][1];
	    
	    newpos[tempcount][0] = (CellSecs[ntop+i]->pos)[k][0];
	    newpos[tempcount][1] = (CellSecs[ntop+i]->pos)[k][1] - j *y_cell_diff;
	    newpos[tempcount][2] = (CellSecs[ntop+i]->pos)[k][2];
	    
	    newpots[tempcount] = (CellSecs[ntop+i]->site_pots)[k];

	    tempcount++;
	  }
	}
      }
      
      //printf("# Nrem  %d    Nremnew %d    tempcount %d\n", *Nrem, Nremnew, tempcount);
	if(struc_out != 0)
	{
	  for(i=0; i<Ntotnew; i++)
	  {
	    if(newsiteinfo[i][0] == 0 && newsiteinfo[i][1] == 0)
	    {
	      fprintf(out, "%lf	%lf\n", newpos[i][0], newpos[i][1]);
	    }
	  }
	  fprintf(out, "\n");
	  for(i=0; i<Ntotnew; i++)
	  {
	    if(newsiteinfo[i][0] == 0 && newsiteinfo[i][1] == 1)
	    {
	      fprintf(out, "%lf	%lf\n", newpos[i][0], newpos[i][1]);
	    }
	  }
	}

      free(siteinfo[0]); free(siteinfo);
      free(pos[0]); free(pos);
      free(site_pots);
      
      *Nrem = Nremnew;
      *Ntot = Ntotnew;
      
      (System->pos) = newpos;
      (System->site_pots) = newpots;
      (System->siteinfo) = newsiteinfo;
      
      
      //LEAD STRUCTURE!
      int num_leads = 2 + ntop + nbot;    
      
	      //left and right leads
	      for(i=0; i < 2; i++)
	      {
		Leads[i] = (RectRedux *)malloc(sizeof(RectRedux));
		(Leads[i]->geo) = (System->geo);
		(Leads[i]->length) = (System->length);
		(Leads[i]->length2) = 1;
		(Leads[i]->Nrem) = (int *)malloc(sizeof(int));
		(Leads[i]->Ntot) = (int *)malloc(sizeof(int));
		simpleRibbonGeo (Leads[i], NULL, 0, NULL);
	      }
	      
	      //top probes
	      for(i=2; i < 2 + ntop; i++)
	      {
		Leads[i] = (RectRedux *)malloc(sizeof(RectRedux));
		(Leads[i]->geo) = tbgeo;
		(Leads[i]->length) = toppw[i-2];
		(Leads[i]->length2) = 1;
		(Leads[i]->Nrem) = (int *)malloc(sizeof(int));
		(Leads[i]->Ntot) = (int *)malloc(sizeof(int));
		simpleRibbonGeo (Leads[i], NULL, 0, NULL);
		
		swapxy((Leads[i]->pos), *(Leads[i]->Ntot));
		
		
		for(j=0; j< *(Leads[i]->Ntot); j++)
		{
		(Leads[i]->pos)[j][0] += (xleft + toppx[i-2]*x_cell_diff);
		(Leads[i]->pos)[j][1] += (ytop + toppc[i-2]*y_cell_diff);
		}
		
		
	      }
	      
	      //bottom probes
	      for(i=2 + ntop; i < 2 + ntop + nbot; i++)
	      {
		Leads[i] = (RectRedux *)malloc(sizeof(RectRedux));
		(Leads[i]->geo) = tbgeo;
		(Leads[i]->length) = botpw[i-2-ntop];
		(Leads[i]->length2) = 1;
		(Leads[i]->Nrem) = (int *)malloc(sizeof(int));
		(Leads[i]->Ntot) = (int *)malloc(sizeof(int));
		simpleRibbonGeo (Leads[i], NULL, 0, NULL);
		
		swapxy((Leads[i]->pos), *(Leads[i]->Ntot));
		
		for(j=0; j< *(Leads[i]->Ntot); j++)
		{
		  (Leads[i]->pos)[j][0] += (xleft + botpx[i-2-ntop]*x_cell_diff);
		  (Leads[i]->pos)[j][1] += (ybot - botpc[i-2-ntop]*y_cell_diff);
		}
		
	      }
	      
	      //shift vectors for left and right leads
	      (params->shift_vecs) = createNonSquareDoubleMatrix(num_leads, 3);
	      (params->shift_vecs)[0][0] = -((System->pos)[2*(System->length)][0] - (System->pos)[0][0]);
	      (params->shift_vecs)[1][0] = (System->pos)[2*(System->length)][0] - (System->pos)[0][0];
	      (params->shift_vecs)[0][1] = 0;
	      (params->shift_vecs)[1][1] = 0;
	      
	      //shift vectors for top and bottom leads
	      for(i=2; i < 2 + ntop; i++)
	      {
		(params->shift_vecs)[i][0] = 0;
		(params->shift_vecs)[i][1] = y_cell_diff;
	      }
	      
	      for(i=2 + ntop; i < 2 + ntop + nbot; i++)
	      {
		(params->shift_vecs)[i][0] = 0;
		(params->shift_vecs)[i][1] = -y_cell_diff;
	      }
	      
	      //shift left and right leads accordingly to desired positions
	      for(i=0; i<*(Leads[0]->Ntot); i++)
	      {
		(Leads[0]->pos)[i][0] += (params->shift_vecs)[0][0];      
	      }
	      
	      for(i=0; i<*(Leads[1]->Ntot); i++)
	      {
		(Leads[1]->pos)[i][0] += (System->length2)*(params->shift_vecs)[1][0];      
	      }
    
	if(struc_out != 0)
	{
	  	  fprintf(out, "\n");
		  fprintf(out, "\n");

		      for(i=0; i<num_leads; i++)
		      {
			for(j=0; j<*(Leads[i]->Ntot); j++)
			{
			  fprintf(out, "%lf	%lf\n", (Leads[i]->pos)[j][0], (Leads[i]->pos)[j][1]);
			}
			fprintf(out, "\n");
		      }
// 		      	  	  fprintf(out, "\n");
// 
// 		       for(i=0; i<num_leads; i++)
// 		      {
// 			for(j=0; j<*(Leads[i]->Ntot); j++)
// 			{
// 			  fprintf(out, "%lf	%lf\n", (Leads[i]->pos)[j][0] + (params->shift_vecs)[i][0], (Leads[i]->pos)[j][1]  + (params->shift_vecs)[i][1]);
// 			}
// 			fprintf(out, "\n");
// 		      }
		      
	}
  
  
  if(struc_out != 0)
  {
    fclose(out);
  }
  
}









//BILAYER ROUTINES

void simpleBilayerGeo (RectRedux *SiteArray, void *p, int struc_out, char *filename)
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
	  
	bilayer_para *params = (bilayer_para *)p;
	  
	  //atomic coordinates and the atoms that are in and out, similar to previous cases
	  int tot_sites = 4*length*length2;
	  double **site_coords = createNonSquareDoubleMatrix(tot_sites, 3);
	  double *site_pots = createDoubleArray(tot_sites);
	  int **siteinfo = createNonSquareIntMatrix(tot_sites, 2);
	  double *subpots = params->subpots;
          
	  
         

          
	  double xstart, ystart;
	  int isclean, l, m, i;
	  int *Nrem = (SiteArray->Nrem);
	  int *Ntot = (SiteArray->Ntot);
	  

          
	  //in AB stackd, (type_shift=2), dimers are siteinfo[1] = 1, 2 (L1-B, L2-A)
	  
	  *Ntot = tot_sites;

	  if(geo==0)
	  {
	    //layer 1
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
	    
		//layer 2
		if(params->type_shift==0)
		{
			params->shift_vec[0]=0.0;
			params->shift_vec[1]=0.0;
		}
		
		if(params->type_shift==1)
		{
			params->shift_vec[0]=0.5;
			params->shift_vec[1]=1/(2*sqrt(3));
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
			//layer 2
			if(params->type_shift==0)
			{
				params->shift_vec[0]=0.0;
				params->shift_vec[1]=0.0;
			}
			
			if(params->type_shift==1)
			{
				params->shift_vec[1]=0.5;
				params->shift_vec[0]=1/(2*sqrt(3));
			}
			if(params->type_shift==12)
			{
                                params->shift_vec[1]=0.0;
				params->shift_vec[0]=-1/(sqrt(3));
			}
		
	  }
	  
	  
	//layer 2
		     
	    for(m=0; m<2*length*length2; m++)
	    {      
		    site_coords[2*length*length2 + m][0] = site_coords[m][0] + params->shift_vec[0];
		    site_coords[2*length*length2 + m][1] = site_coords[m][1] + params->shift_vec[1] ;
		    site_coords[2*length*length2 + m][2] = params->zsep;
		    siteinfo[2*length*length2 + m][1] = siteinfo[m][1] + 2;
	     
	    }
	     
	     
	      	      
	int tempint=0, tempint2=0;

	
	for(m=0; m<*Ntot; m++)
	{
		if(siteinfo[m][0] == 0)
		{
			tempint2++;
		}
	}
	
	*Nrem = tempint2;
	
	
	//site potentials by sublattice
	for(m=0; m<*Ntot; m++)
	{
		if(siteinfo[m][1] >=0 && siteinfo[m][1] < 4)
		{
// 			printf("# %d\n", siteinfo[m][1]);
			site_pots[m] = subpots[siteinfo[m][1]];
// 			printf("%lf\n", site_pots[m]);
		}
	}

	
	if(struc_out == 1)
	{
		for(l=0; l<*Ntot; l++)
		{
			if(siteinfo[l][0] == 0 && siteinfo[l][1] == 0)
			{  
				fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			}
		}
		fprintf(out, "\n");
		for(l=0; l<*Ntot; l++)
		{
			if(siteinfo[l][0] == 0 && siteinfo[l][1] == 1)
			{  
				fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			}
		}
		fprintf(out, "\n");
		for(l=0; l<*Ntot; l++)
		{
			if(siteinfo[l][0] == 0 && siteinfo[l][1] == 2)
			{  
				fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			}
		}
		fprintf(out, "\n");
		for(l=0; l<*Ntot; l++)
		{
			if(siteinfo[l][0] == 0 && siteinfo[l][1] == 3)
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



void BLGPotentials (RectRedux *SiteArray, void *p, int struc_out, char *filename)
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


            
            blgpots_para *params = (blgpots_para *)p;
            bilayer_para *blgparams = (bilayer_para *)(params->BLGpara) ;
            int seed = (params->seed);
            
            srand(time(NULL) + seed);
	  
	  //atomic coordinates and the atoms that are in and out, similar to previous cases
	  int tot_sites = 4*length*length2;
	  double **site_coords = createNonSquareDoubleMatrix(tot_sites, 3);
	  double *site_pots = createDoubleArray(tot_sites);
	  int **siteinfo = createNonSquareIntMatrix(tot_sites, 2);
	  double *subpots = params->subpots;
          double *subconcs = params->subconcs;
          int buffer_rows = params->buffer_rows;

          double *nosubpots = blgparams->subpots;
	  
	  double xstart, ystart;
	  int isclean, l, m, i;
	  int *Nrem = (SiteArray->Nrem);
	  int *Ntot = (SiteArray->Ntot);
	  
	  //in AB stackd, (type_shift=2), dimers are siteinfo[1] = 1, 2 (L1-B, L2-A)
	  
	  *Ntot = tot_sites;
          
          


	  if(geo==0)
	  {
	    //layer 1
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
	    
		//layer 2
		if(blgparams->type_shift==0)
		{
			blgparams->shift_vec[0]=0.0;
			blgparams->shift_vec[1]=0.0;
		}
		
		if(blgparams->type_shift==1)
		{
			blgparams->shift_vec[0]=0.5;
			blgparams->shift_vec[1]=1/(2*sqrt(3));
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
			//layer 2
			if(blgparams->type_shift==0)
			{
				blgparams->shift_vec[0]=0.0;
				blgparams->shift_vec[1]=0.0;
			}
			
			if(blgparams->type_shift==1)
			{
				blgparams->shift_vec[1]=0.5;
				blgparams->shift_vec[0]=1/(2*sqrt(3));
			}
			
			if(blgparams->type_shift==12)
			{
                                blgparams->shift_vec[1]=0.0;
				blgparams->shift_vec[0]=-1/(sqrt(3));
			}
		
	  }
	  
	  
	//layer 2
		     
	    for(m=0; m<2*length*length2; m++)
	    {      
		    site_coords[2*length*length2 + m][0] = site_coords[m][0] + blgparams->shift_vec[0];
		    site_coords[2*length*length2 + m][1] = site_coords[m][1] + blgparams->shift_vec[1] ;
		    site_coords[2*length*length2 + m][2] = blgparams->zsep;
		    siteinfo[2*length*length2 + m][1] = siteinfo[m][1] + 2;
	     
	    }
	     
	     
	      	      
	int tempint=0, tempint2=0;

	
	for(m=0; m<*Ntot; m++)
	{
		if(siteinfo[m][0] == 0)
		{
			tempint2++;
		}
	}
	
	*Nrem = tempint2;
	
        
        double temprandnum;
        for(l=0; l<buffer_rows*2*length; l++)
        {
            if(siteinfo[l][1] >=0 && siteinfo[l][1] < 4)
            {
                    site_pots[l] = nosubpots[siteinfo[l][1]];
            }
        }
        for(l=2*length*length2; l< 2*length*length2 + buffer_rows*2*length; l++)
        {
            if(siteinfo[l][1] >=0 && siteinfo[l][1] < 4)
            {
                    site_pots[l] = nosubpots[siteinfo[l][1]];
            }
        }
        
        for(l=buffer_rows*2*length; l<2*length*length2 - 2*buffer_rows*length; l++)
        {
            site_pots[l] = 0.0;
            temprandnum = myRandNum(0.0, 1.0);
            if(temprandnum < subconcs[ siteinfo[l][1]])
		    site_pots[l] = subpots[ siteinfo[l][1]];
            temprandnum = myRandNum(0.0, 1.0);
            if(temprandnum < subconcs[ siteinfo[l][1]])
		    site_pots[l] = subpots[ siteinfo[l][1]];
            
        }
        for(l=2*length*length2 + buffer_rows*2*length; l< 4*length*length2 - 2*buffer_rows*length; l++)
        {
            site_pots[l] = 0.0;
            temprandnum = myRandNum(0.0, 1.0);
            if(temprandnum < subconcs[ siteinfo[l][1]])
		    site_pots[l] = subpots[ siteinfo[l][1]];
        }
        
        for(l=2*length*length2 - 2*buffer_rows*length; l<2*length*length2 ; l++)
        {
            if(siteinfo[l][1] >=0 && siteinfo[l][1] < 4)
            {
                    site_pots[l] = nosubpots[siteinfo[l][1]];
            }
        }
        for(l=4*length*length2 - 2*buffer_rows*length; l<4*length*length2; l++)
        {
            if(siteinfo[l][1] >=0 && siteinfo[l][1] < 4)
            {
                    site_pots[l] = nosubpots[siteinfo[l][1]];
            }
        }
        
        


	
	if(struc_out == 1)
	{
		for(l=0; l<*Ntot; l++)
		{
			if(siteinfo[l][0] == 0 && siteinfo[l][1] == 0)
			{  
				fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			}
		}
		fprintf(out, "\n");
		for(l=0; l<*Ntot; l++)
		{
			if(siteinfo[l][0] == 0 && siteinfo[l][1] == 1)
			{  
				fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			}
		}
		fprintf(out, "\n");
		for(l=0; l<*Ntot; l++)
		{
			if(siteinfo[l][0] == 0 && siteinfo[l][1] == 2)
			{  
				fprintf(out, "%lf	%lf\n", site_coords[l][0], site_coords[l][1]);
			}
		}
		fprintf(out, "\n");
		for(l=0; l<*Ntot; l++)
		{
			if(siteinfo[l][0] == 0 && siteinfo[l][1] == 3)
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






//Generation of generic multilayer devices
//Each layer is generated independently from different single layer systems
//Combined with specified z coordinates, rescaling, shifts and rotations into one RectRedux

void customMultilayer (RectRedux *SiteArray, void *p, int struc_out, char *filename)
{
    
    multilayer_para *params = (multilayer_para *)p;
    int num_layers = (params->num_layers);
    
    int i, j, k, l;
    int ntot=0, nrem=0, ntc, nrc;
    int *Nrem = (SiteArray->Nrem);
    int *Ntot = (SiteArray->Ntot);
    int need_pert=0;
    FILE *out;
    
    RectRedux *LayerCells[num_layers];
    
    double origin[2], origintemp[2];
    origin[0]=0.0; origin[1]=0.0;
    origintemp[0]=0.0; origintemp[1]=0.0;
    
    
    if ( (params->origin) == 1)
    {
        
        if ( (SiteArray->geo) ==0)
        {
            origin[0] = (int) (SiteArray->length2/2) - 0.5* (  (SiteArray->length % 2) ); 
            origin[1] = (sqrt(3)/2.0) * (int) (SiteArray->length/2)   ;
        }
        
    }
    
    
    //Individual layer structures
    for(i=0; i< num_layers; i++)
    {
        LayerCells[i] = (RectRedux *)malloc(sizeof(RectRedux));
        (LayerCells[i]->geo) = (params->geo)[i];
        (LayerCells[i]->length) = (params->length)[i];
        (LayerCells[i]->length2) = (params->length2)[i];
        (LayerCells[i]->Nrem) = (int *)malloc(sizeof(int));
        (LayerCells[i]->Ntot) = (int *)malloc(sizeof(int));
        
        //generate individual layers
        ((params->layerfn)[i]) (LayerCells[i], (params->layerpara)[i], 0, NULL) ;
        
        //total numbers
        nrem += *(LayerCells[i]->Nrem);
        ntot += *(LayerCells[i]->Ntot);
        
        if((LayerCells[i]->pert_pos) != NULL)
            need_pert=1;
        
    }
    
    *Nrem = nrem;
    *Ntot = ntot;
    
    double **site_coords = createNonSquareDoubleMatrix(ntot, 3);
    double *site_pots = createDoubleArray(ntot);
    int **siteinfo = createNonSquareIntMatrix(ntot, 2);
    double **pert_pos = NULL;
    double x1, y1, z1, cost, sint;
    
    if(need_pert ==1)
    {
        pert_pos = createNonSquareDoubleMatrix(ntot, 3);
    }
    
    ntc=0; 
    
    //new positions of sites in layers
    for(i=0; i< num_layers; i++)
    {
        cost = cos((params->theta)[i]);
        sint = sin((params->theta)[i]);
        
        if ( (params->origin) == 1)
        {
            
            if ( (params->geo)[i] ==0)
            {
                origintemp[0] = (int) ((params->length2)[i]/2) - 0.5* (  ((params->length)[i] % 2) ); 
                origintemp[1] = (sqrt(3)/2.0) * (int) ((params->length)[i]/2)   ;
            }
            
        }
        
        
        
        for(j=0; j< *(LayerCells[i]->Ntot); j++)
        {
            x1 = (LayerCells[i]->pos)[j][0] + (params->delta)[i][0] - origintemp[0];
            y1 = (LayerCells[i]->pos)[j][1] + (params->delta)[i][1] - origintemp[1];
            z1 = (LayerCells[i]->pos)[j][2] + (params->delta)[i][2] ;
            
            site_coords[ntc + j][0] =  origin[0] + (params->epsilon)[i] * (x1*cost -y1*sint);
            site_coords[ntc + j][1] =  origin[1] + (params->epsilon)[i] * (x1*sint +y1*cost);
            site_coords[ntc + j][2] =  z1;
                        
            site_pots[ntc + j] = (LayerCells[i]->site_pots)[j] ;
            siteinfo[ntc + j][0] = (LayerCells[i]->siteinfo)[j][0] ;
            siteinfo[ntc + j][1] = (LayerCells[i]->siteinfo)[j][1] ;
            
            
            if(need_pert==1)
            {
                if((LayerCells[i]->pert_pos) != NULL)
                {
                    x1 = (LayerCells[i]->pert_pos)[j][0] + (params->delta)[i][0] ;
                    y1 = (LayerCells[i]->pert_pos)[j][1] + (params->delta)[i][1] ;
                    z1 = (LayerCells[i]->pert_pos)[j][2] + (params->delta)[i][2] ;
                }
                
                pert_pos[ntc + j][0] =  (params->epsilon)[i] * (x1*cost -y1*sint);
                pert_pos[ntc + j][1] =  (params->epsilon)[i] * (x1*sint +y1*cost);
                pert_pos[ntc + j][2] =  z1;
            }
            
        }
        
        ntc+=*(LayerCells[i]->Ntot);
    }
 
    
    //Free memory associated with individual layers
    
    for(i=0; i< num_layers; i++)
    {
        free((LayerCells[i]->siteinfo)[0]); free((LayerCells[i]->siteinfo));
        free((LayerCells[i]->pos)[0]); free((LayerCells[i]->pos));
        free((LayerCells[i]->site_pots));
        
        if((LayerCells[i]->pert_pos) != NULL)
        {
            free((LayerCells[i]->pert_pos)[0]); free((LayerCells[i]->pert_pos));
        }
    }

    
    (SiteArray->pos) = site_coords;
    (SiteArray->site_pots) = site_pots;
    (SiteArray->siteinfo) = siteinfo;
    (SiteArray->pert_pos) = pert_pos;
    
    
    if(struc_out != 0)
    {
        out = fopen(filename, "w");
        
        if(struc_out == 1)
        {
            ntc=0; 
            for(i=0; i<num_layers; i++)
            {
                for(j=0; j<*(LayerCells[i]->Ntot); j++)
                {
                    
                    if(need_pert == 0)
                    {
                        if(siteinfo[ntc+j][0] == 0 )
                        {  
                            fprintf(out, "%lf	%lf   %lf\n", site_coords[ntc+j][0], site_coords[ntc+j][1], site_coords[ntc+j][2]);
                        }
                    }
                    if(need_pert == 1)
                    {
                        if(siteinfo[ntc+j][0] == 0 )
                        {  
                            fprintf(out, "%lf	%lf   %lf\n", pert_pos[ntc+j][0], pert_pos[ntc+j][1], pert_pos[ntc+j][2]);
                        }
                    }
                }
                ntc+=*(LayerCells[i]->Ntot);
                fprintf(out, "\n");
            }
            
        }
        fclose(out);
    }
}






















