#include "devices.h"
#include "mapping.h"
#include <stdio.h>


//LDOS out, can be split by sublattice and impurity sites
void LDOSMapOut(RectRedux *DeviceCell, double *ldoses, char *out)
{
    int Nrem = *(DeviceCell -> Nrem);
    int length1 = (DeviceCell->length);
    int length2 = (DeviceCell->length2);
    double **pos = (DeviceCell->pos);
    int **siteinfo = (DeviceCell->siteinfo);
    double *site_pots = (DeviceCell->site_pots);

   
    int num=0, i, j, k, numtodate=0;
    FILE *output;
    char filename[160];
    sprintf(filename, "%s.dmap", out);
    output = fopen(filename, "w");
    
    
    		
	      for(i=0; i<2*length1*length2; i++)
	      {
		if(siteinfo[i][0]==0 )
		{
		  //if(siteinfo[i][1]==0 && site_pots[i] == 0.0)
		  //{
		    fprintf(output, "%lf	%lf	%lf\n", pos[i][0], pos[i][1], ldoses[numtodate] );
		  //}
		  numtodate++;
		}
		
	      }
// 	      fprintf(output, "\n");
	      
// 	      numtodate=0;
// 	      for(i=0; i<2*length1*length2; i++)
// 	      {
// 		if(siteinfo[i][0]==0)
// 		{
// 		  if(siteinfo[i][1]==1  && site_pots[i] == 0.0)
// 		  {
// 		    fprintf(output, "%lf	%lf	%lf\n", pos[i][0], pos[i][1], ldoses[numtodate] );
// 		  }
// 		  numtodate++;
// 		}
// 		
// 	      }
// 	      	      fprintf(output, "\n");
// 		      	      numtodate=0;
// 
// 	       for(i=0; i<2*length1*length2; i++)
// 	      {
// 		if(siteinfo[i][0]==0 )
// 		{
// 		  if(siteinfo[i][1]==0 && site_pots[i] != 0.0 )
// 		  {
// 		    fprintf(output, "%lf	%lf	%lf\n", pos[i][0], pos[i][1], ldoses[numtodate] );
// 		  }
// 		  numtodate++;
// 		}
// 		
// 	      }
// 	      fprintf(output, "\n");
// 	      	      numtodate=0;
// 
// 	      
// 	      for(i=0; i<2*length1*length2; i++)
// 	      {
// 		if(siteinfo[i][0]==0 )
// 		{
// 		  if(siteinfo[i][1]==1 && site_pots[i] != 0.0)
// 		  {
// 		    fprintf(output, "%lf	%lf	%lf\n", pos[i][0], pos[i][1], ldoses[numtodate] );
// 		  }
// 		  numtodate++;
// 		}
// 		
// 	      }
// 	      fprintf(output, "\n"); 	      
		      

		fclose(output);
    
    
}


void CurrentMapOut(RectRedux *DeviceCell, double **conds, int *indices, int **neigh, char *out)
{
    int Nrem = *(DeviceCell -> Nrem);
    int length1 = (DeviceCell->length);
    int length2 = (DeviceCell->length2);
    double **pos = (DeviceCell ->pos);
    
    double xstart, ystart, x2, y2;
    int num=0, i, j, k, numtodate=0;
    
    int chain, cell4, atom, type;
    double condxavg, condyavg;
    
    FILE *output;
    char filename[160];
    sprintf(filename, "%s.cmap", out);
    output = fopen(filename, "w");
    
    for(i=0; i<2*length1*length2; i++)
    {
      
	  condxavg=0; condyavg=0;

	  for(j=0; j<3; j++)
	  {
	    if(indices[neigh[i][j]] > -1 && neigh[i][j] != -1)
	    {

	      
	      condxavg += (pos[neigh[i][j]][0]-pos[i][0])*conds[indices[i]][j];
	      condyavg += (pos[neigh[i][j]][1]-pos[i][1])*conds[indices[i]][j];
	      
 
	      
	      //specific bond currents for awk codes later
	      fprintf(output, "##	%lf	%lf	%lf	%lf %e %d	%d\n", pos[i][0], pos[i][1], pos[neigh[i][j]][0], pos[neigh[i][j]][1], conds[indices[i]][j], indices[i], indices[neigh[i][j]]);
	      
	      
	    }
	  }
	  //total current vector from each site
	  fprintf(output, "%e	%e %e %e\n", pos[i][0], pos[i][1], condxavg, condyavg);

      
	if(indices[i] < 0)
	{	  
	  fprintf(output, "%e	%e %e %e\n", pos[i][0], pos[i][1], 0.0, 0.0);
	}
      
      
    }
    fclose(output);
    
  
    
    
}
