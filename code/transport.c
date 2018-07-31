#include "transport.h"
#include "test.h"


//mode=0 - just transmissions!
//mode=1 - transmissions and ldos/current maps - can be generalised further later
void genTransmissions(double _Complex En, RectRedux *DeviceCell, RectRedux **Leads, cnxProfile *cnxp, 
		      cellDivision *cellinfo, hoppingfunc *hoppingfn, void *hoppingparams,
		      lead_para *leadsparams, trans_params *tpara, int mode, double *ldoses, double ***currents)
{
    //important definitions and quantities
//     int length1 = (DeviceCell->length);
//     int length2 = (DeviceCell->length2);

//     printf("#started transmission calc\n");
    int geo = (DeviceCell->geo);
    
    int cell1dim = (cellinfo->cell1dim);
    int num_leads = tpara->num_leads;
    
    int this_cell, dim, dim_old=0, dim_new;
    int Nrem = *(DeviceCell->Nrem);
    int Ntot = *(DeviceCell->Ntot);
    double _Complex **g_old, g_new, g00inv;
    double _Complex **g_sys_r, **g_sys_a, **SigmaR, **SigmaA, **Gamma, **temp1, **temp2, *gii, **gi1, **g1i;
    double _Complex **Gamma1, **Gamma2, **gsr, **gsa;
    int d1, d2, s1, s2;
    
    int i, j, k, l, m, n;
    
    int devicemode=0, devicemode2=0; 
    
    if(mode==0)
    {
      devicemode = 0; 
      devicemode2 = 0;
    }
    
    if(mode==1)
    {
     devicemode = 1;
     devicemode2 = 1;
    }
    
    gii=NULL; gi1=NULL; g1i=NULL;
    if(mode ==1)
    {
      gii = createCompArray(Ntot);
      gi1 = createNonSquareMatrix(Ntot, cell1dim);
      g1i = createNonSquareMatrix(cell1dim, Ntot);
    }
          
    
    SigmaR = createSquareMatrix(cell1dim);
    SigmaA = createSquareMatrix(cell1dim);
    g_sys_r = createSquareMatrix(cell1dim);
    g_sys_a = createSquareMatrix(cell1dim);
    Gamma = createSquareMatrix(cell1dim);

      leadfunction *leadfn = (leadfunction *)(leadsparams->leadsfn);

      //lead sigmas 
	 (leadfn)(En, DeviceCell, Leads, cellinfo, leadsparams, SigmaR);
	 
	 

      //calculate retarded GF of system 
         genDeviceGF(En, DeviceCell, cnxp, cellinfo, hoppingfn, hoppingparams, devicemode, devicemode2, g_sys_r, gii, gi1, SigmaR);
// 	 printEMatrix(g_sys_r, cell1dim);
// 	 listNonZero(g_sys_r, cell1dim, cell1dim);
	 
      //calculate advanced quantities
// 	 if( (tpara->TRsym) == 0)
// 	 {
	   for(i=0; i<cell1dim; i++)
	   {
	     for(j=0; j<cell1dim; j++)
	     {
		SigmaA[i][j] = conj(SigmaR[j][i]);
		g_sys_a[i][j] = conj(g_sys_r[j][i]);		
	     }
	   }
	   
	   if(mode==1)
	   {
	      for(i=0; i<cell1dim; i++)
	      {
		  for(j=0; j<Ntot; j++)
		  {
		    g1i[i][j] = conj(gi1[j][i]);
		  }
	      }
		
	    }
	   
	/*   
	 }
	 else if( (tpara->TRsym) == 1)
	 {
	   
	   if(mode==1)
	   {
	      devicemode2 = 2;
	   }
	   
	    (leadfn)(creal(En) - I*cimag(En), DeviceCell, Leads, cellinfo, leadsparams, SigmaA);
	    genDeviceGF(creal(En) - I*cimag(En), DeviceCell, cnxp, cellinfo, hoppingfn, hoppingparams, devicemode, devicemode2, g_sys_a, NULL, g1i, SigmaA);
	 }*/
	 
// 	 printEMatrix(g_sys_a, cell1dim);
	 
	 //Gamma
	  for(i=0; i<cell1dim; i++)
	  {
	    for(j=0; j<cell1dim; j++)
	    {
	      Gamma[i][j] = I*(SigmaR[i][j] - SigmaA[i][j]);
	    }
	  }
// 	  printEMatrix(Gamma, cell1dim);
//   	  listNonZero(SigmaR, cell1dim, cell1dim);
      FreeMatrix(SigmaR); FreeMatrix(SigmaA); 
      
      s1=0; 
      for(i=0; i<num_leads; i++)
      {
	d1=(cellinfo->lead_dims)[i];
	s2=0;
	Gamma1 = createSquareMatrix(d1);
	MatrixCopyPart(Gamma, Gamma1, s1, s1, 0, 0, d1, d1);
	
	for(j=0; j< num_leads; j++)
	{
	  d2=(cellinfo->lead_dims)[j];
	  gsr = createNonSquareMatrix(d1, d2);
	  MatrixCopyPart(g_sys_r, gsr, s1, s2, 0, 0, d1, d2);
	  Gamma2 = createSquareMatrix(d2);
	  MatrixCopyPart(Gamma, Gamma2, s2, s2, 0, 0, d2, d2);
	  gsa = createNonSquareMatrix(d2, d1);
	  MatrixCopyPart(g_sys_a, gsa, s2, s1, 0, 0, d2, d1);
	  
	  temp1=createNonSquareMatrix(d1, d2);
	  MatrixMultNS(Gamma1, gsr, temp1, d1, d1, d2);
	  temp2=createNonSquareMatrix(d1, d2);
	  MatrixMultNS(temp1, Gamma2, temp2, d1, d2, d2);
	  FreeMatrix(temp1);
	  temp1=createSquareMatrix(d1);
	  MatrixMultNS(temp2, gsa, temp1, d1, d2, d1);
	  FreeMatrix(temp2);
	  
	  (tpara->transmissions)[i][j] = creal(MatrixTrace(temp1, d1));
	  FreeMatrix(temp1);
	  s2+=d2;
	  FreeMatrix(gsr); FreeMatrix(Gamma2); FreeMatrix(gsa);
	  
// 	  printf("#%d	%d	%d	%d	%lf\n", i, j, d1, d2, (tpara->transmissions)[i][j]);
	}
	
	s1+=d1;
	FreeMatrix(Gamma1);
      }
        double *bshifts = createDoubleArray(3);
	double current_mag;
	double _Complex hopmag;
	
	FILE *bigdump;
	char bigfile[200];
	


      if(mode==1)
      {
	
	// a way to bulk print out all the bond currents.
	//for debugging only. see loop below also.
// 	for(i=0; i<num_leads; i++)
// 	{
// 	  sprintf(bigfile, "current_dump_l%d", i);
// 	  bigdump = fopen(bigfile, "w");
// 	  fclose(bigdump);
// 	}
	
	for(i=0; i<Ntot; i++)
	{
	  ldoses[i] = -cimag(gii[i])/M_PI;
	
	  s1=0; d1=0;
	  for(k=0; k<num_leads; k++)
	  {
	      d1=(cellinfo->lead_dims)[k];
 	      currents[k][i][0] = 0.0;
	      currents[k][i][1] = 0.0;
	      
	      for(l=0; l<(cnxp->site_cnxnum)[i]; l++)
	      {
		j = (cnxp->site_cnx)[i][l];
		hopmag = hoppingfn(DeviceCell, DeviceCell, j, i, bshifts, hoppingparams);
		current_mag=0;
		
		for(m=0; m<d1; m++)
		{
		  for(n=0; n<d1; n++)
		  {
		    current_mag += cimag(hopmag * gi1[i][s1+m]  * Gamma[s1+m][s1+n]  *g1i[s1+n][j]);
		  }
		}
// 		sprintf(bigfile, "current_dump_l%d", k);
// 		bigdump = fopen(bigfile, "a");
//  		if((DeviceCell->pos)[j][0]>=(DeviceCell->pos)[i][0])
// 		  fprintf(bigdump, "%lf	%lf	%lf	%lf	%.15e\n", (DeviceCell->pos)[i][0], (DeviceCell->pos)[i][1], (DeviceCell->pos)[j][0],(DeviceCell->pos)[j][1], current_mag);
// 		fclose(bigdump);
		
	
		
		if(current_mag>0)
		{
		  currents[k][i][0] += current_mag * ( (DeviceCell->pos)[j][0] -(DeviceCell->pos)[i][0]  );
		  currents[k][i][1] += current_mag * ( (DeviceCell->pos)[j][1] -(DeviceCell->pos)[i][1]  );
		}
		
	      }
	      
	      s1+=d1;
	      
	  }
	      
	    
	  
	  
	}
	
	
	
	
      }
      
      
      
//       for(i=0; i<Ntot; i++)
//       {
// 	printf("%lf	%lf	%e\n", (DeviceCell->pos)[i][0], (DeviceCell->pos)[i][1], -cimag(gii[i]));
//       }
       
	
       if(mode ==1)
	{
	  free(gii); 
	  FreeMatrix(g1i); 
	  FreeMatrix(gi1); 
	}
         
      free(bshifts);
      FreeMatrix(Gamma);
      FreeMatrix(g_sys_r); FreeMatrix(g_sys_a);
        
  

}

//mode indicates either single sweep (mode=0) or double sweep (mode=1)
//(mode=2) for double sweep with reversed sweeping directions (for cases yet to be determined, patched GF?)
//mode2 determines what is calculated in the reverse sweep (if mode>0) :
//0 - Gii only (e.g. for DOS maps only)
//1 - Gii and Gi1 (e.g. for LDOS maps and current maps with time reversal symmetry)
//2 - G1i only (for advanced GFs in second current map sweep where TRS is broken)

void genDeviceGF(double _Complex En, RectRedux *DeviceCell, cnxProfile *cnxp, 
		      cellDivision *cellinfo, hoppingfunc *hoppingfn, void *hoppingparams, int mode, int mode2, 
		      double _Complex **Gon, double _Complex *Gdiags, double _Complex **Goff, double _Complex **Sigma)
{

  int geo = (DeviceCell->geo);
    
  int num_cells = (cellinfo->num_cells);
  
  int this_cell, dim, dim_old=0, dim_new;
  double _Complex **g_old,  **g00inv, **V12, **V21, **smallSigma, **temp1, **temp2;
  int i, j, k, l, index1, index2, this0, last0;
  
  int cell_start, cell_end, cell_iter, it_count;
  int dim1  = (cellinfo->cell_dims)[0];
  int Nrem = *(DeviceCell->Nrem);
    int Ntot = *(DeviceCell->Ntot);

  double *bshifts = createDoubleArray(3);
  
  double _Complex ***allgs, **gtemp, **off_old, **off_new, **t1, **t2, **gprev;
  
  
  if(mode == 0 || mode == 1)
  {
    cell_start = num_cells-1;
    cell_end = 0;
    cell_iter = -1;
  }
  
  if(mode == 2)
  {
    cell_start = 0 ;
    cell_end = num_cells-1;
    cell_iter = 1;
  }
  
  
  if(mode>0)
  {
    //requests memory for an array of pointers to matrices
    allgs = (double _Complex ***)malloc(num_cells * sizeof(double _Complex **));
    
  }


  //calculate system GF	
  //backwards recursive sweep!
  for(it_count=0; it_count<(cellinfo->num_cells); it_count ++)
  {
	this_cell = cell_start + it_count*cell_iter;
  
//    	printf("########\n#cell %d\n", this_cell);
	
// 	printf("### backwards sweep, cell %3d of %3d\n", this_cell, (cellinfo->num_cells));
    //generate disconnected cell GF
	dim = (cellinfo->cell_dims)[this_cell];
	g00inv = createSquareMatrix(dim);
	
	this0=(cellinfo->starting_index)[this_cell];
	if(it_count>0)
	{
	  last0=(cellinfo->starting_index)[this_cell-cell_iter];
	  
	  V21 = createNonSquareMatrix(dim, dim_old);
	  V12 = createNonSquareMatrix(dim_old, dim);
	}
	
	
	//onsites
	for(i=0; i<dim; i++)
	{
	  index1 = (cellinfo->cells_site_order)[this0 +i];
	  
	  g00inv[i][i] = En - (DeviceCell->site_pots)[index1] - hoppingfn(DeviceCell, DeviceCell, index1, index1, bshifts, hoppingparams);
          
          
          //magnitude is negative, hence += to add to the inverse GF
          if(DeviceCell->cap_pots != NULL)
          {
              g00inv[i][i] += I* (DeviceCell->cap_pots)[index1];
          }
	  
	  
	  
	}

	//internal (& external?) hoppings
	for(i=0; i<dim; i++)
	{
	  
	  index1 = (cellinfo->cells_site_order)[this0 +i];
	  for(j=0; j<(cnxp->site_cnxnum)[index1]; j++)
	  {
	    index2 = (cnxp->site_cnx)[index1][j];
	    
	    //is this neighbour in the same cell?
	    if( (cellinfo->sites_by_cell)[index2] == this_cell)
	    {
		//what is its index within the cell?
		k=0;
		for(l=0; l<dim; l++)
		{
		  if( index2 == (cellinfo->cells_site_order)[this0 +l])
		  {
		    k=l;
		  }
		}
		//printf("neighbours	%d	%d\n", i, k);
 		g00inv[i][k] -= hoppingfn(DeviceCell, DeviceCell, index1, index2, bshifts, hoppingparams);
	      
	    }
	    
	    
	    //is this neighbour in the previous cell?
	    if( (cellinfo->sites_by_cell)[index2] == (this_cell-cell_iter))
	    {
		//what is its index within the cell?
		k=0;
		for(l=0; l<dim_old; l++)
		{
		  if( index2 == (cellinfo->cells_site_order)[last0 +l])
		  {
		    k=l;
		  }
		}
 		V21[i][k] = hoppingfn(DeviceCell, DeviceCell, index1, index2, bshifts, hoppingparams);
		V12[k][i] = hoppingfn(DeviceCell, DeviceCell, index2, index1, bshifts, hoppingparams);
	      
	    }
	    
	    
	  }
	}

	  	//check non-zero elements at the cell stage
//  	  	printf("#g00inv\n");
//  	  	listNonZero(g00inv, dim, dim);
// // 	  
// 	  	if(it_count>0)
// 	  	{
// 	  	  printf("#V21\n");
// 	  	  listNonZero(V21, dim, dim_old);
// 	  	  
// 	  	  printf("#V12\n");
// 	  	  listNonZero(V12, dim_old, dim);
// 	  	}
      
      
      	  

    //generate self energy term from previous cells
	smallSigma=createSquareMatrix(dim);
	
	if(it_count>0)
	{
	  //previous cells self energy: smallSigma = V21 g_old V12
	    temp1 = createNonSquareMatrix(dim, dim_old);
	    MatrixMultNS(V21, g_old, temp1, dim, dim_old, dim_old);
	    MatrixMultNS(temp1, V12, smallSigma, dim, dim_old, dim);
	    FreeMatrix(temp1);
	  
	  
	}
      
	
	
	
    //account for lead self energies if this is cell 0
	if(this_cell==0)
	{
	  temp1 = createSquareMatrix(dim);
	  MatrixAdd(smallSigma, Sigma, temp1, dim);
	  MatrixCopy(temp1, smallSigma, dim);
// 	  	listNonZero(Sigma, dim, dim);

	  FreeMatrix(temp1);
	}
	
	temp1 = createSquareMatrix(dim);
	MatrixSubtract(g00inv, smallSigma, temp1, dim);
	if(it_count>0)
	{
	  FreeMatrix(g_old);
	}
	g_old = createSquareMatrix(dim);
	InvertMatrixGSL(temp1, g_old, dim);
// 	printf("# it: %d, %e %e\n", it_count, creal(g_old[1][1]), cimag(g_old[1][1]));
	FreeMatrix(temp1); FreeMatrix(g00inv); FreeMatrix(smallSigma);
	
	if(it_count>0)
	{
	  FreeMatrix(V12); FreeMatrix(V21);
	}
	
	//save the temporary edges if needed for the dual sweep
	if(mode>0)
	{
	  allgs[this_cell] = createSquareMatrix(dim);
	  MatrixCopy(g_old, allgs[this_cell], dim);
	}
	
	dim_old=dim;
// 		listNonZero(g_old, dim, dim);
		
	
	
      
  }//end of first sweep 

  MatrixCopy(g_old, Gon, dim1);
  FreeMatrix(g_old);

  
  //perform forward sweep
  if(mode>0)
  {

    //setups and initialisations before running loop
    
    //copy info about cell 0
    off_old = createNonSquareMatrix(dim1, dim1);
    MatrixCopyPart(Gon, off_old, 0, 0, 0, 0, dim1, dim1);

    if(mode2<2)
    {
      gprev = createSquareMatrix(dim1);
      MatrixCopyPart(Gon, gprev, 0, 0, 0, 0, dim1, dim1);
    }
//           	        	printf("ok-tz %d\n", this0);

    if(mode2==1)
    {
	for(i=0; i<dim1; i++)
	{
// 	  printf("ok %d / %d\n", i, dim1);
	  Gdiags[(cellinfo->cells_site_order)[i]] = Gon[i][i];
// 		printf("ok %d / %d	Nr %d\n", i, dim1, Nrem);

	
	  for(j=0; j<dim1; j++)
	  {
// 	    printf("okj %d / %d siteorder %d\n", j, dim1, (cellinfo->cells_site_order)[i]);
	    Goff[(cellinfo->cells_site_order)[i]][j] = Gon[i][j];
	  }
// 	  printf("ok %d / %d\n", i, dim1);
	}
            

    }
		
    if(mode2==2)
    {
	for(i=0; i<dim1; i++)
	{
	  for(j=0; j<dim1; j++)
	  {           	        	
	    Goff[i][(cellinfo->cells_site_order)[j]] = Gon[i][j];
	  }
	}
    }

//     printf ("\n"); 
    //insert second sweep here
    for(it_count=1; it_count<(cellinfo->num_cells); it_count ++)
    {
	this_cell = cell_end - it_count*cell_iter;
	
// 	printf("### forwards sweep, cell %3d of %3d\n", this_cell, (cellinfo->num_cells));

    
	dim = (cellinfo->cell_dims)[this_cell];
	dim_old = (cellinfo->cell_dims)[this_cell + cell_iter];
	
	this0=(cellinfo->starting_index)[this_cell];
	last0=(cellinfo->starting_index)[this_cell+cell_iter];
	
	
	//g_old is the current cell but before update
	g_old = createSquareMatrix(dim);
	MatrixCopyPart(allgs[this_cell], g_old, 0,0,0,0, dim, dim);
	
	//gtemp is the fully connected and updated version of the current cell
	if(mode2<2)
	{
	  gtemp = createSquareMatrix(dim);
	}
	
	//gprev is the fully connected version of the previous cell

	
	//off_old 
	    //is filled at the end of the iteration by copying off_new
	    //has dimension (dim_old, dim1) (or vice versa) and is the connected version
	
	//off_new
	if(mode2 == 1)
	{
	    off_new = createNonSquareMatrix(dim, dim1);
	}
	if(mode2 == 2)
	{
	    off_new = createNonSquareMatrix(dim1, dim);
	}
	    
	    
	//generate V matrices as above 
		V12 = createNonSquareMatrix(dim_old, dim);
		V21 = createNonSquareMatrix(dim, dim_old);
	    
		for(i=0; i<dim; i++)
		{
		  
		  index1 = (cellinfo->cells_site_order)[this0 +i];
		  for(j=0; j<(cnxp->site_cnxnum)[index1]; j++)
		  {
		    index2 = (cnxp->site_cnx)[index1][j];
		    
		    //is this neighbour in the previous cell?
		    if( (cellinfo->sites_by_cell)[index2] == (this_cell + cell_iter))
		    {
			//what is its index within the cell?
			k=0;
			for(l=0; l<dim_old; l++)
			{
			  if( index2 == (cellinfo->cells_site_order)[last0 +l])
			  {
			    k=l;
			  }
			}
			
			V21[i][k] = hoppingfn(DeviceCell, DeviceCell, index1, index2, bshifts, hoppingparams);
			V12[k][i] = hoppingfn(DeviceCell, DeviceCell, index2, index1, bshifts, hoppingparams);
		      
		    }
		    
		    
		  }
		}
		
	
	
	if(mode2<2)
	{
	
	
		t1 = createNonSquareMatrix(dim, dim_old);
		MatrixMultNS(g_old, V21, t1, dim, dim, dim_old);
		
		//update off-diagonals
		if(mode2==1)
		{
		  MatrixMultNS(t1, off_old, off_new, dim, dim_old, dim1);
		}
		
		
	//update diagonals
		t2 = createNonSquareMatrix(dim, dim_old);
		MatrixMultNS(t1, gprev, t2, dim, dim_old, dim_old);
		FreeMatrix(t1);
		t1 = createNonSquareMatrix(dim, dim);
		MatrixMultNS(t2, V12, t1, dim, dim_old, dim);
		FreeMatrix(t2); 
		t2= createSquareMatrix(dim);
		MatrixMult(t1, g_old, t2, dim);
		FreeMatrix(t1);
		MatrixAdd(g_old, t2, gtemp, dim);
		FreeMatrix(t2);
      
	    
	//copying and freeing
		FreeMatrix(g_old);
		FreeMatrix(V12);
		FreeMatrix(V21);
		FreeMatrix(gprev);
		
		gprev = createSquareMatrix(dim);
		MatrixCopy(gtemp, gprev, dim);
		
		for(i=0; i<dim; i++)
		{
		  Gdiags[(cellinfo->cells_site_order)[this0 +i]] = gtemp[i][i];
		  
		  for(j=0; j<dim1; j++)
		  {
		    Goff[(cellinfo->cells_site_order)[this0 +i]][j] = off_new[i][j];
		  }
		}
		
		FreeMatrix(gtemp);
		FreeMatrix(off_old);
		
		off_old=createNonSquareMatrix(dim, dim1);
		MatrixCopyPart(off_new, off_old, 0, 0, 0, 0, dim, dim1);
		FreeMatrix(off_new);
		
	}
	
	if(mode2 == 2)
	{
		t1 = createNonSquareMatrix(dim1, dim);
		MatrixMultNS(off_old, V12, t1, dim1, dim_old, dim);
		MatrixMultNS(t1, g_old, off_new, dim1, dim, dim);
		
		FreeMatrix(t1);
		FreeMatrix(g_old);
		FreeMatrix(V12);
		FreeMatrix(V21);
	
		for(i=0; i<dim; i++)
		{
		  for(j=0; j<dim1; j++)
		  {
		    Goff[j][(cellinfo->cells_site_order)[this0 +i]] = off_new[j][i];
		  }
		}
		FreeMatrix(off_old);
		off_old=createNonSquareMatrix(dim1, dim);
		MatrixCopyPart(off_new, off_old, 0, 0, 0, 0, dim1, dim);
		FreeMatrix(off_new);
	}
    }
	

    
  }

   free(bshifts);
   
   if(mode>0)
   {
     for(i=0; i<num_cells; i++)
     {
       FreeMatrix(allgs[i]);
     }
     free(allgs);
   }
  
}



//generates bands/projections/weighted band struc assuming a system periodic in x
//mode indicates either bands only (mode=0), projections at given k (mode=1), or an attempt at an unfolded weighted bandstructure (mode=3)
//(mode=2) for double sweep with reversed sweeping directions (for broken time-reversal symmetry cases)
void genKXbandproj(RectRedux *DeviceCell,  hoppingfunc *hoppingfn, void *hoppingparams, int mode,
		      double kx, double *bands, double **projs, double **weights)
{
  
  int geo = (DeviceCell->geo);
  int Nrem = *(DeviceCell->Nrem);
    int Ntot = *(DeviceCell->Ntot);
      int length2 = (DeviceCell->length2);
      int length1 = (DeviceCell->length);


  double _Complex **g_old,  **g00inv, **V12, **V21, **smallSigma, **temp1, **temp2;
  int i, j, k, l, index1, index2, this0, last0;
  
  double *bshifts = createDoubleArray(3);
  
  int *rem_tot_mapping = createIntArray(Nrem);
  int **unit_cell_map = createNonSquareIntMatrix(Nrem, 2);
  double _Complex **Ham = createSquareMatrix(Nrem);
  
  double cell_sep;
  if(geo==0)
  {
    cell_sep = 1.0 * length2;
  }
  if(geo==1)
  {
    cell_sep = sqrt(3) * length2;
  }
  
  
  if(mode == 0 || mode == 1)
  {
  
  }
  
  //create mapping of indices
  j=0;
  for(i=0; i <Ntot; i++)
  {
    if((DeviceCell->siteinfo)[i][0] == 0)
    {
      rem_tot_mapping[j] = i;
      unit_cell_map[j][0] = (i / (2*length1));
      unit_cell_map[j][1] = (i % (2*length1));
      // printf("%d	%d	%d	%d\n", i, j, unit_cell_map[j][0], unit_cell_map[j][1]);
      j++;
    }
  }
  
  
    
  
  
  for(i=0; i<Nrem; i++)
  {
    //onsites (& hopping corrections to these)
    k=rem_tot_mapping[i];
    Ham[i][i] = (DeviceCell->site_pots)[k];   
		  //check that the onsite correction here works properly at some stage
    
     if(DeviceCell->cap_pots != NULL)
          {
              Ham[i][i] -= I* (DeviceCell->cap_pots)[k];
          }
	  
    
    //intra and inter -cell hoppings - this should also fix onsite hopping corrections
    for(j=0; j< Nrem; j++)
    {
      l = rem_tot_mapping[j];
      bshifts[0] = 0;
      //Ham[i][j] += hoppingfn(DeviceCell, DeviceCell, k, l, bshifts, hoppingparams)* cexp(I*kx* (DeviceCell->pos[l][0] + bshifts[0] -DeviceCell->pos[k][0] ) );
      Ham[i][j] += hoppingfn(DeviceCell, DeviceCell, k, l, bshifts, hoppingparams);
      
      //connections to other cells
      bshifts[0] = -cell_sep;
      Ham[i][j] += hoppingfn(DeviceCell, DeviceCell, k, l, bshifts, hoppingparams) * cexp(I*kx*length2); 
       // Ham[i][j] += hoppingfn(DeviceCell, DeviceCell, k, l, bshifts, hoppingparams) * cexp(I*kx* (DeviceCell->pos[l][0] + bshifts[0] -DeviceCell->pos[k][0] ) );
      
      
      bshifts[0] = cell_sep;
      //Ham[i][j] += hoppingfn(DeviceCell, DeviceCell, k, l, bshifts, hoppingparams) * cexp(I*kx* (DeviceCell->pos[l][0] + bshifts[0] -DeviceCell->pos[k][0] ) ); 
      Ham[i][j] += hoppingfn(DeviceCell, DeviceCell, k, l, bshifts, hoppingparams) * cexp(-I*kx*length2); 
    }
  }
  
//    printEMatrix(Ham, Nrem);
    
   double _Complex *evalues, **evecs, **weight, **projects;
    
    evalues = createCompArray(Nrem);
    if(mode ==0)
    {
	EigenvaluesGSL(Ham, evalues, Nrem);
	
 	for(i=0; i< Nrem; i++)
	{
	  bands[i] = evalues[i];
	}
    }
    
    if(mode==1 || mode ==2)
    {
      evecs = createSquareMatrix(Nrem);
      EigenvectorsGSL(Ham, evalues, evecs, Nrem);
      
      for(i=0; i< Nrem; i++)
	{
	  bands[i] = evalues[i];
	}
    }
    
    double kprime;
    
    if(mode==1 || mode ==2)
    {
	projects = createSquareMatrix(Nrem);
	//band index
	for(j=0; j< Nrem; j++)
	{
	  //site
	  for(k=0; k<Nrem; k++)
	  {
	    projects[j][k] = evecs[j][k]  * conj(evecs[j][k]);
	    projs[j][k] = creal(projects[j][k]);
	  }
	}
	
    }
	 
    
    if(mode == 2)
    {
      weight = createNonSquareMatrix(Nrem, length2);
      
      //loop over G vectors to other BZs
      for(i=0; i<length2; i++)
      {
	kprime = kx + i* (2*M_PI/(length2));
	
	//loop over band index
	for(j=0; j< Nrem; j++)
	{
	  weight[j][i] = 0.0;
	  for(k=0; k<Nrem; k++)
	  {
	      for(l=0; l<Nrem; l++)
	      {
		if(unit_cell_map[k][1] == unit_cell_map[l][1])
		{
		    weight[j][i] += conj(evecs[j][k]) * evecs[j][l] * cexp(- I * kprime * 1.0 *(unit_cell_map[k][0] - unit_cell_map[l][0])) / length2;
		  
		}
		
		
	      }
	    
	  }
	  weights[j][i] = creal(weight[j][i]);
	  
	  
	}
	
      }
	
      
      
    }
   

    
    
  free(rem_tot_mapping);  
  FreeMatrix(Ham);  
  free(bshifts);
  free(evalues);
  free(unit_cell_map[0]);
  free(unit_cell_map);
  
  if(mode==1 || mode ==2)
  {
    FreeMatrix(evecs);
  }
  if(mode ==2)
  {
    FreeMatrix(weight);
  }
  if(mode ==1 || mode ==2)
  {
    FreeMatrix(projects);
  }
  
// 	dim = (cellinfo->cell_dims)[this_cell];
// 	g00inv = createSquareMatrix(dim);
// 	
// 	this0=(cellinfo->starting_index)[this_cell];
// 	if(it_count>0)
// 	{
// 	  last0=(cellinfo->starting_index)[this_cell-cell_iter];
// 	  
// 	  V21 = createNonSquareMatrix(dim, dim_old);
// 	  V12 = createNonSquareMatrix(dim_old, dim);
// 	}
// 	
// 	
// 	//onsites
// 	for(i=0; i<dim; i++)
// 	{
// 	  index1 = (cellinfo->cells_site_order)[this0 +i];
// 	  
// 	  g00inv[i][i] = En - (DeviceCell->site_pots)[index1];
// 	  
// 	}
// 
// 	//internal (& external?) hoppings
// 	for(i=0; i<dim; i++)
// 	{
// 	  
// 	  index1 = (cellinfo->cells_site_order)[this0 +i];
// 	  for(j=0; j<(cnxp->site_cnxnum)[index1]; j++)
// 	  {
// 	    index2 = (cnxp->site_cnx)[index1][j];
// 	    
// 	    //is this neighbour in the same cell?
// 	    if( (cellinfo->sites_by_cell)[index2] == this_cell)
// 	    {
// 		//what is its index within the cell?
// 		k=0;
// 		for(l=0; l<dim; l++)
// 		{
// 		  if( index2 == (cellinfo->cells_site_order)[this0 +l])
// 		  {
// 		    k=l;
// 		  }
// 		}
// 		//printf("neighbours	%d	%d\n", i, k);
//  		g00inv[i][k] -= hoppingfn(DeviceCell, DeviceCell, index1, index2, bshifts, hoppingparams);
// 	      
// 	    }
// 	    
// 	    
// 	    //is this neighbour in the previous cell?
// 	    if( (cellinfo->sites_by_cell)[index2] == (this_cell-cell_iter))
// 	    {
// 		//what is its index within the cell?
// 		k=0;
// 		for(l=0; l<dim_old; l++)
// 		{
// 		  if( index2 == (cellinfo->cells_site_order)[last0 +l])
// 		  {
// 		    k=l;
// 		  }
// 		}
//  		V21[i][k] = hoppingfn(DeviceCell, DeviceCell, index1, index2, bshifts, hoppingparams);
// 		V12[k][i] = hoppingfn(DeviceCell, DeviceCell, index2, index1, bshifts, hoppingparams);
// 	      
// 	    }
// 	    
// 	    
// 	  }
// 	}
	

      	  

      
}










//has been edited to accept 2 device cell structures, so can calculate hopping between leads and device also
//bshifts allows rigid shifting of one cell, e.g. for neighbouring cell calculations
double _Complex simpleTB(RectRedux *aDeviceCell, RectRedux *bDeviceCell, int a, int b, double *bshifts, void *hoppingparams)
{
  gen_hop_params *para = (gen_hop_params *)hoppingparams; 
  int num_neigh = para->num_neigh;
  double _Complex *hops = para->hops;
  double  *NN_lowdis = para->NN_lowdis;
  double  *NN_highdis = para->NN_highdis;
  double  *NN_shifts = para->NN_shifts;
  double  *NN_zmin= para->NN_zmin;
  double  *NN_zmax = para->NN_zmax;
  double t0;
  double kpar = para->kpar;
  double _Complex ans=0.0;
  double x1 = (aDeviceCell->pos)[a][0], y1 = (aDeviceCell->pos)[a][1];
  double x2 = (bDeviceCell->pos)[b][0] + bshifts[0], y2 = (bDeviceCell->pos)[b][1] + bshifts[1];
  double y2p, distp;
  
  double dist = sqrt(pow(x2-x1, 2.0) + pow(y2-y1, 2.0));
  double zdiff = fabs((bDeviceCell->pos)[b][2] - (aDeviceCell->pos)[a][2]);
  double ycelldist;
  int i;
//   printf("#%lf	%lf\n", dist, zdiff);
  
  //if((para->isperiodic)==0)
  //{
  
    //which coupling to use
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(dist >= (para->NN_lowdis[i]) && dist < (para->NN_highdis[i]))
      {
		if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
		{ 
			t0 = hops[i];
// 			printf("#hopping type %d of %d used\n", i, num_neigh);
// 			
// 			printf("#type 0: %lf - %lf, z: %lf - %lf\n", (para->NN_lowdis[0]) , (para->NN_highdis[0]), NN_zmin[0], NN_zmax[0]);
// 			printf("#type 1: %lf - %lf, z: %lf - %lf\n", (para->NN_lowdis[1]) , (para->NN_highdis[1]), NN_zmin[1], NN_zmax[1]);
			
		}
      }
      
      if(dist == 0.0)
      {
	t0 += NN_shifts[i];
      }
    }
    ans=t0;
 
  

      //note the += to allow more than one connection between the same atoms (or their images)
  if((para->isperiodic)==1)
  {
    //set separation to up and down cells
    //this is only sensible for even-indexed ribbons, but will run with odd results for odd indices
    //this
    if((aDeviceCell->geo)==0)
    {
	ycelldist = (aDeviceCell->length)*sqrt(3)/2;
    }
    if((aDeviceCell->geo)==1)
    {
	ycelldist = (aDeviceCell->length)*0.5;
    }
    
    //check if b is in cell above
    y2p = y2 + ycelldist;
    distp= sqrt(pow(x2-x1, 2.0) + pow(y2p-y1, 2.0));
    
//     if(distp > (para->NN_lowdis) && distp < (para->NN_highdis))
//     {
//      ans+=t0*cexp(-I*kpar);    
//     }
    
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(distp > (para->NN_lowdis[i]) && distp < (para->NN_highdis[i]))
      {
	      if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
	      {
		t0 = hops[i];
	      }
      }
    }
    ans+=t0*cexp(-I*kpar); 
    
    
    //check if b is in cell below
    y2p = y2 - ycelldist;
    distp= sqrt(pow(x2-x1, 2.0) + pow(y2p-y1, 2.0));
    
//     if(distp > (para->NN_lowdis) && distp < (para->NN_highdis))
//     {
//      ans+=t0*cexp(I*kpar); 
//     }
    
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(distp > (para->NN_lowdis[i]) && distp < (para->NN_highdis[i]))
      {
	      if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
	      {
			t0 = hops[i];
	      }
      }
    }
    ans+=t0*cexp(I*kpar); 
    
  }
  //printf("# hopping %d	%d: %lf	%lf\n", a, b, ans, dist);
  return ans;
  
}


//A simple implementation of hopping parameters for a strained graphene system
double _Complex strainedTB(RectRedux *aDeviceCell, RectRedux *bDeviceCell, int a, int b, double *bshifts, void *hoppingparams)
{
  gen_hop_params *para = (gen_hop_params *)hoppingparams; 
  int num_neigh = para->num_neigh;
  double _Complex *hops = para->hops;
  double  *NN_lowdis = para->NN_lowdis;
  double  *NN_highdis = para->NN_highdis;
  double  *NN_shifts = para->NN_shifts;
  double  *NN_zmin= para->NN_zmin;
  double  *NN_zmax = para->NN_zmax;
  double t0;
  double kpar = para->kpar;
  double _Complex ans=0.0;
  double x1 = (aDeviceCell->pos)[a][0], y1 = (aDeviceCell->pos)[a][1];
  double x2 = (bDeviceCell->pos)[b][0] + bshifts[0], y2 = (bDeviceCell->pos)[b][1] + bshifts[1];
  double y2p, distp;
  double basedist[] = {1/sqrt(3), 1, 2/sqrt(3)};
  double beta = 3.37;
  
  double dist = sqrt(pow(x2-x1, 2.0) + pow(y2-y1, 2.0));
  double zdiff = fabs((bDeviceCell->pos)[b][2] - (aDeviceCell->pos)[a][2]);
  double ycelldist;
  int i;
  double x1r, x2r, y1r, y2r, z1r, z2r, distr, distpr, y2pr;
  
    if(aDeviceCell->pert_pos != NULL)
    {
        x1r = (aDeviceCell->pert_pos)[a][0];
        y1r = (aDeviceCell->pert_pos)[a][1];
        z1r = (aDeviceCell->pert_pos)[a][2];
    }
    else
    {
        x1r=x1;
        y1r=y1;
        z1r=0.0;
    }
    if(bDeviceCell->pert_pos != NULL)
    {
        x2r = (bDeviceCell->pert_pos)[b][0] + bshifts[0];
        y2r = (bDeviceCell->pert_pos)[b][1] + bshifts[1];
        z2r = (bDeviceCell->pert_pos)[b][2];
    }
    else
    {
        x2r=x2;
        y2r=y2;
        z2r=0.0;
    }
    distr = sqrt(pow(x2r-x1r, 2.0) + pow(y2r-y1r, 2.0) + pow(z2r-z1r, 2.0));
  

  
    //which coupling to use
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(dist >= (para->NN_lowdis[i]) && dist < (para->NN_highdis[i]))
      {
		if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
		{ 
			t0 = hops[i];
                        
                        //strain
                        if(aDeviceCell->pert_pos != NULL || bDeviceCell-> pert_pos != NULL)
                        {
                            t0=t0*exp(-beta * ( (distr/basedist[i]) - 1.0) );
                        }
		}
      }
      
      if(dist == 0.0)
      {
	t0 += NN_shifts[i];
      }
    }
    ans=t0;
    
        
 
  

      //note the += to allow more than one connection between the same atoms (or their images)
  if((para->isperiodic)==1)
  {
    //set separation to up and down cells
    //this is only sensible for even-indexed ribbons, but will run with odd results for odd indices
    if((aDeviceCell->geo)==0)
    {
	ycelldist = (aDeviceCell->length)*sqrt(3)/2;
    }
    if((aDeviceCell->geo)==1)
    {
	ycelldist = (aDeviceCell->length)*0.5;
    }
    
    //check if b is in cell above
    y2p = y2 + ycelldist;
    y2pr = y2r + ycelldist;

    distp= sqrt(pow(x2-x1, 2.0) + pow(y2p-y1, 2.0));
    distpr = sqrt(pow(x2r-x1r, 2.0) + pow(y2pr-y1r, 2.0) + pow(z2r-z1r, 2.0));


    
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(distp > (para->NN_lowdis[i]) && distp < (para->NN_highdis[i]))
      {
	      if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
	      {
		t0 = hops[i];
                
                //strain
                if(aDeviceCell->pert_pos != NULL || bDeviceCell-> pert_pos != NULL)
                {
                    t0=t0*exp(-beta * ( (distpr/basedist[i]) - 1.0) );
                }
	      }
      }
    }
    ans+=t0*cexp(-I*kpar); 
    
    
    //check if b is in cell below
    y2p = y2 - ycelldist;
    y2pr = y2r - ycelldist;
    distp= sqrt(pow(x2-x1, 2.0) + pow(y2p-y1, 2.0));
    distpr = sqrt(pow(x2r-x1r, 2.0) + pow(y2pr-y1r, 2.0) + pow(z2r-z1r, 2.0));

    
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(distp > (para->NN_lowdis[i]) && distp < (para->NN_highdis[i]))
      {
	      if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
	      {
			t0 = hops[i];
                        
                        //strain
                        if(aDeviceCell->pert_pos != NULL || bDeviceCell-> pert_pos != NULL)
                        {
                            t0=t0*exp(-beta * ( (distpr/basedist[i]) - 1.0) );
                        }
	      }
      }
    }
    ans+=t0*cexp(I*kpar); 
    
  }
  return ans;
  
}



double _Complex peierlsTB(RectRedux *aDeviceCell, RectRedux *bDeviceCell, int a, int b, double *bshifts, void *hoppingparams)
{
  gen_hop_params *para = (gen_hop_params *)hoppingparams; 
  int num_neigh = para->num_neigh;
  double _Complex *hops = para->hops;
  double  *NN_lowdis = para->NN_lowdis;
  double  *NN_highdis = para->NN_highdis;
  double  *NN_shifts = para->NN_shifts;
  double  *NN_zmin= para->NN_zmin;
  double  *NN_zmax = para->NN_zmax;
  double t0;
  double kpar = para->kpar;
  double _Complex ans=0.0;
  double x1 = (aDeviceCell->pos)[a][0], y1 = (aDeviceCell->pos)[a][1];
  double x2 = (bDeviceCell->pos)[b][0] + bshifts[0], y2 = (bDeviceCell->pos)[b][1] + bshifts[1];
  
  double zdiff = fabs((bDeviceCell->pos)[b][2] - (aDeviceCell->pos)[a][2]);
  
  double y2p, distp;
  double dist = sqrt(pow(x2-x1, 2.0) + pow(y2-y1, 2.0));
  double ycelldist;
  int i;
  
  
    //which coupling to use
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(dist > (para->NN_lowdis[i]) && dist < (para->NN_highdis[i]))
      {
	      if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
	      {
		t0 = hops[i];
	      }
      }
      
      //a==b condition has been added to distinguish between different layer cases
      if(dist == 0.0 && a==b)
      {
	t0 += NN_shifts[i];
      }
      
      
    }
    ans=t0*graphenePeierlsPhase(x1, y1, x2, y2, para->gauge, para->Btes, para->restrics, para->limits);     

    
      //note the += to allow more than one connection between the same atoms (or their images)
  if((para->isperiodic)==1)
  {
    //set separation to up and down cells
    //this is only sensible for even-indexed ribbons, but will run with odd results for odd indices
    if((aDeviceCell->geo)==0)
    {
	ycelldist = (aDeviceCell->length)*sqrt(3)/2;
    }
    if((aDeviceCell->geo)==1)
    {
	ycelldist = (aDeviceCell->length)*0.5;
    }
    
    //check if b is in cell above
    y2p = y2 + ycelldist;
    distp= sqrt(pow(x2-x1, 2.0) + pow(y2p-y1, 2.0));
    

    
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(distp > (para->NN_lowdis[i]) && distp < (para->NN_highdis[i]))
      {
	      if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
	      {
			t0 = hops[i];
	      }
      }
    }
    ans+=t0*graphenePeierlsPhase(x1, y1, x2, y2p, para->gauge, para->Btes, para->restrics, para->limits)*cexp(-I*kpar);    
    
    
    
    //check if b is in cell below
    y2p = y2 - ycelldist;
    distp= sqrt(pow(x2-x1, 2.0) + pow(y2p-y1, 2.0));
    

    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(distp > (para->NN_lowdis[i]) && distp < (para->NN_highdis[i]))
      {
	      if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
	      {
			t0 = hops[i];
	      }
      }
    }
    ans+=t0*graphenePeierlsPhase(x1, y1, x2, y2p, para->gauge, para->Btes, para->restrics, para->limits)*cexp(I*kpar);    
    
    
  }
  return ans;
}






//returns peierls phase factor for two sites in graphene lattice
//note assumes the graphene lattice constant, so should be generalised for nongraphene
//gauge = 0 -> phase along x
//gauge = 1 -> phase along y

//gauges 2 and 3 returns peierls phase factor for two sites in graphene lattice in a Hall bar setup
//there is a twist in the gauge to allow periodic leads in both directions
//gauge = 2 -> fields in leads - phase along y initially
//gauge = 3 -> no fields in leads - phase along x initially

double _Complex graphenePeierlsPhase(double x1, double y1, double x2, double y2, int gauge, double BTesla, int *res, double **limits)
{
 
  double xb, yb, delx, dely, midx, midy;
  double beta = BTesla * 1.46262E-5;
  double phase, weight;
  dely=y2-y1; delx=x2-x1;
  midx=x1 + delx/2; midy = y1+dely/2;
  double xrenorm=0;
  

  //limits for gauge=0 and gauge=1 if finite cutoff
  if (gauge == 0 || gauge == 1)
  {
	if(res[0] == 0)
	{
	  xb=midx;
	}
	else if(res[0] == 1)
	{
	  xb=midx;
	  if(xb<limits[0][0])
	  {
	    xb=limits[0][0]  ;
	  }
	  if(xb>limits[0][1])
	  {
	    xb=limits[0][1]  ;
	  }
	  
	}
	//printf("#x %lf	%lf\n", x, x1);
	
	if(res[1] == 0)
	{
	  yb=midy;
	}
	else if(res[0] == 1)
	{
	  yb=midy;
	  if(yb<limits[1][0])
	  {
	    yb=limits[1][0];
	  }
	  if(yb>limits[1][1])
	  {
	    yb=limits[1][1];
	  }
	}
  }
  
  if(gauge == 0)
  {
  //  phase = beta*(x*dely + 0.5*delx*dely);
  // phase = beta*(x*dely );
    phase = beta*xb*dely;
  }
  
  else if(gauge == 1)
  {
   // phase = beta*(y*delx + 0.5*delx*dely);
   //phase = beta*(y*delx );
        phase = -beta*yb*delx;
    
  }
  
  double phi=0, phii=0, phij=0, weight2=0;
  double l6, m1, m2, m3;
  if(gauge == 2)
  {
//     xrenorm = limits[0][1] + (limits[0][2] - limits[0][1])/2;
       // xrenorm = 0;

//     l6 = limits[0][4] + (limits[0][5] - limits[0][4])/2;
//     
//     m1 = limits[0][1] + (limits[0][4] - limits[0][1])/2;
//     m2 = limits[0][5] + (limits[0][2] - limits[0][5])/2;
//     m3 = limits[0][4] + (limits[0][5] - limits[0][4])/2;
//     
//     m1= limits[0][1];
//     m2 = limits[0][2];
//     
//     
//     if(midx <= limits[0][0])
//     {
//       weight = 0.0;
//       xrenorm = m1;
//     }
//     if(midx > limits[0][0] && midx < limits[0][1])
//     {
//       weight = (midx - limits[0][0])/(limits[0][1] - limits[0][0]);
//       xrenorm = m1;
//     }
//     if(midx >= limits[0][1] && midx <= limits[0][4])
//     {
//       weight = 1.0;
//       xrenorm = m1;
//     }
//     if(midx > limits[0][4] && midx <= l6)
//     {
//       weight = 1.0 - (midx - limits[0][4])/(l6 - limits[0][4]);
//       xrenorm = m1;
//     }
//     if(midx > l6 && midx < limits[0][5])
//     {
//       weight = (midx - l6)/(limits[0][5] - l6);
//       xrenorm = m2;
//     }
//     if(midx >= limits[0][5] && midx <= limits[0][2])
//     {
//       weight = 1.0;
//       xrenorm = m2;
//     }
//     if(midx > limits[0][2] && midx < limits[0][3])
//     {
//       weight = 1.0 - (midx - limits[0][2])/(limits[0][3] - limits[0][2]);
//       xrenorm = m2;
//     }
//     if(midx >= limits[0][3])
//     {
//       weight = 0.0;
//       xrenorm = m2;
//     }
//     phase = beta*( (weight -1)*midy*delx + weight*(midx-xrenorm)*dely);
    
    
    xrenorm = limits[0][0]; 
    //xrenorm=0;
    
//     if(midx <= limits[0][0])
//     {
//       weight = 0.0;
//       phase = beta*( (weight -1)*midy*delx + weight*(midx-xrenorm)*dely);
//     }
//     if(midx > limits[0][0] && midx < limits[0][1])
//     {
//       weight = (midx - limits[0][0])/(limits[0][1] - limits[0][0]);
//       phase = beta*( (weight -1)*midy*delx + weight*(midx-xrenorm)*dely);
//     }
//     if(midx >= limits[0][1] && midx <= limits[0][2])
//     {
//       weight = 1.0;
//       phase = beta*( (weight -1.0)*midy*delx + weight*(midx-xrenorm)*dely);
//     }
//     if(midx > limits[0][2] && midx < limits[0][3])
//     {
//       weight = 1.0 - (midx - limits[0][2])/(limits[0][3] - limits[0][2]);
//       phase = beta*( (weight -1.0)*(midy*delx - (limits[0][3]-xrenorm)*dely)  + weight*(midx-xrenorm)*dely);
//     }
//     if(midx >= limits[0][3])
//     {
//       weight = 0.0;
//       phase = beta*( (weight -1.0)*(midy*delx - (limits[0][3]-xrenorm)*dely)  + weight*(midx-xrenorm)*dely);
//     }
    
    
        xrenorm = 0; 

     if(midx <= limits[0][0])
    {
      weight = 0.0;
      phase = beta*( (weight -1)*midy*delx + weight*(midx-xrenorm)*dely);
    }
    
    double len;
    
    if(midx > limits[0][0] && midx < limits[0][1])
    {
      len = limits[0][1] - limits[0][0];
      
      //smooth transition phase
      phase = (beta / (2*len*len*len)) * (2 *dely *x1 * pow(x1 - limits[0][0], 2) * (3*len - 2*x1 + 2* limits[0][0]) 
				      -delx*dely*(len*len*len + 4 * pow(x1 - limits[0][0], 2) * (4*x1 - limits[0][0]) - 6*len * (3*x1*x1 - 4*x1*limits[0][0] +limits[0][0]* limits[0][0]))
				      -2*delx*(len-x1+limits[0][0])*(len*len + len*(x1 - limits[0][0]) - 2*(4*x1*x1 -5*x1*limits[0][0] + limits[0][0] * limits[0][0]))*y1
				      -4*pow(delx,4) * (y1 + dely)
				      +2*pow(delx, 3) * (3*len - 8*x1 + 6*limits[0][0]) * (y1 + dely)
				      + 6*pow(delx, 2) * (3*len*x1 - 4*x1*x1 - 2*len*limits[0][0] + 6*x1*limits[0][0] - 2*pow(limits[0][0], 2))* (y1 + dely) );
      
      //linear transition phase
//       phase = (beta / (2*len)) * (2*x1*dely*(x1 - limits[0][0]) + 2*delx*delx*(y1+dely)
// 				      + delx*( dely*(4*x1 - limits[0][0] - limits[0][1]) + 2*y1*(2*x1 - limits[0][1])));
//       
      
    }
    
    
    if(midx >= limits[0][1] && midx <= limits[0][2])
    {
      weight = 1.0;
      phase = beta*( (weight -1.0)*midy*delx + weight*(midx-xrenorm)*dely);
    }
    
    
    if(midx > limits[0][2] && midx < limits[0][3])
    {
      //weight = 1.0 - (midx - limits[0][2])/(limits[0][3] - limits[0][2]);
     // phase = beta*( (weight -1.0)*(midy*delx - (limits[0][3]-xrenorm)*dely)  + weight*(midx-xrenorm)*dely);
            len = limits[0][3] - limits[0][2];

	    
	//smooth transition phase
      phase = (beta / (2*len*len*len)) * (2 *dely *x1 * pow(len - x1 + limits[0][2], 2) * (len + 2*x1 - 2* limits[0][2]) 
				      +delx*dely*(len*len*len + 4 * pow(x1 - limits[0][2], 2) * (4*x1 - limits[0][2]) - 6*len * (3*x1*x1 - 4*x1*limits[0][2] +limits[0][2]* limits[0][2]))
				      +2*delx*(x1-limits[0][2])*( -9*len*x1 + 8*x1*x1 + 3*len*limits[0][2] -10*x1*limits[0][2] + 2*limits[0][2]*limits[0][2]  )*y1
				      +4*pow(delx,4) * (y1 + dely)
				      -2*pow(delx, 3) * (3*len - 8*x1 + 6*limits[0][2]) * (y1 + dely)
				      - 6*pow(delx, 2) * (3*len*x1 - 4*x1*x1 - 2*len*limits[0][2] + 6*x1*limits[0][2] - 2*pow(limits[0][2], 2))* (y1 + dely) );
      
     //linear transition phase
//       phase = (beta / (2*len)) * (2*x1*dely*( limits[0][3] - x1) - 2*delx*delx*(y1+dely)
// 				      + delx*( dely*(-4*x1 + limits[0][2] + limits[0][3]) + 2*y1*(limits[0][2] -2*x1)));
//       
     
      
    }
    if(midx >= limits[0][3])
    {
      weight = 0.0;
      //phase = beta*( (weight -1.0)*(midy*delx - (limits[0][3]-xrenorm)*dely)  + weight*(midx-xrenorm)*dely);
      phase = beta*( (weight -1.0)*(midy*delx)  + weight*(midx-xrenorm)*dely);

    }
    

	

    
  }
    
  if(gauge == 3)  //this is not running correctly yet
  {
      
      //x-limits
      xb=midx;
      if(xb<limits[0][0])
      {
	xb=limits[0][0]  ;
      }
      if(xb>limits[0][3])
      {
	xb=limits[0][3]  ;
      }
    
 
  
      //ylimits  -limits[1][1]? and [1][2] not used in this code, but indices consistent with x limits
      yb=midy;
      if(yb<limits[1][0])
      {
	yb=limits[1][0];
      }
      if(yb>limits[1][3])
      {
	yb=limits[1][3];
      }
      
      
      if(midx <= limits[0][0])
      {
	weight = 0.0;
      }
      if(midx > limits[0][0] && midx < limits[0][1])
      {
	weight = (midx - limits[0][0])/(limits[0][1] - limits[0][0]);
      }
      if(midx >= limits[0][1] && midx <= limits[0][2])
      {
	weight = 1.0;
      }
      if(midx > limits[0][2] && midx < limits[0][3])
      {
	weight = 1.0 - (midx - limits[0][2])/(limits[0][3] - limits[0][2]);
      }
      if(midx >= limits[0][3])
      {
	weight = 0.0;
      }
      phase = beta*( (weight -1)*midx*dely + weight*midy*delx);
      
      
  
  }
  
  
  
  
  
  
//  printf("# %lf (%lf,  %lf) phase %+e\n", x, delx, dely, phase);
//   if(delx > 0)
   //printf("%lf	%lf	%e	\n", midx, midy, phase);
  
  return cexp(2 * M_PI * I * phase);  
  
  
}



//returns peierls phase factor for two sites in graphene lattice in a Hall bar setup
//note assumes the graphene lattice constant, so should be generalised for nongraphene
//there is a twist in the gauge to allow periodic leads in both directions
//gauge = 3 -> fields in leads
//gauge = 4 -> no fields in leads
double _Complex grapheneHallPhase(double x1, double y1, double x2, double y2, int gauge, double BTesla, int *res, double **limits)
{
 
  double xb, yb, delx, dely, midx, midy;
  double beta = BTesla * 1.46262E-5;
  double phase, weight;
  dely=y2-y1; delx=x2-x1;
  midx=x1 + delx/2; midy = y1+dely/2;
  

  if(gauge == 3)
  {
    if(midx <= limits[0][0])
    {
      weight = 0.0;
    }
    if(midx > limits[0][0] && midx < limits[0][1])
    {
      weight = (midx - limits[0][0])/(limits[0][1] - limits[0][0]);
    }
    if(midx >= limits[0][1] && midx <= limits[0][2])
    {
      weight = 1.0;
    }
    if(midx > limits[0][2] && midx < limits[0][3])
    {
      weight = 1.0 - (midx - limits[0][2])/(limits[0][3] - limits[0][2]);
    }
    if(midx >= limits[0][3])
    {
      weight = 0.0;
    }
    phase = beta*( (weight -1)*midy*delx + weight*midx*dely);
  }
    
  if(gauge == 4)
  {
      
      //x-limits
      xb=midx;
      if(xb<limits[0][0])
      {
	xb=limits[0][0]  ;
      }
      if(xb>limits[0][3])
      {
	xb=limits[0][3]  ;
      }
    
 
  
      //ylimits  -limits[1][3] and [1][2] not used in this code, but indices consistent with x limits
      yb=midy;
      if(yb<limits[1][0])
      {
	yb=limits[1][0];
      }
      if(yb>limits[1][3])
      {
	yb=limits[1][3];
      }
      
      
      if(midx <= limits[0][0])
      {
	weight = 0.0;
      }
      if(midx > limits[0][0] && midx < limits[0][1])
      {
	weight = (midx - limits[0][0])/(limits[0][1] - limits[0][0]);
      }
      if(midx >= limits[0][1] && midx <= limits[0][2])
      {
	weight = 1.0;
      }
      if(midx > limits[0][2] && midx < limits[0][3])
      {
	weight = 1.0 - (midx - limits[0][2])/(limits[0][3] - limits[0][2]);
      }
      if(midx >= limits[0][3])
      {
	weight = 0.0;
      }
      phase = beta*( (weight -1)*midx*dely + weight*midy*delx);
      
      
  
  }
    
  
//  printf("# %lf (%lf,  %lf) phase %+e\n", x, delx, dely, phase);
  
  return cexp(2 * M_PI * I * phase);  
  
  
}


//generalised lead Sigmas (based on simple2leads)
void multipleLeads (double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, lead_para *params, double _Complex **Sigma)
{
  
    int leadloop, dim1, dim1a,  dimcounta=0, lcount;
    double _Complex **ginv, **V12, **V21, **g00, **SL, **SR, **VLD, **VDL, **smallSigma, **temp1;
    double elemerr=1.0e-15;

    int num_leads = (cellinfo->num_leads);
    int i, j, k;
  
    double *bshifts0 = createDoubleArray(3);
    hoppingfunc *hopfn = (hoppingfunc *)(params->hopfn);
    
    
    for(leadloop=0; leadloop < num_leads; leadloop++)
    {
  
      //generate leads SGFs
  
	  dim1 = *(LeadCells[leadloop]->Nrem);
	  ginv = createSquareMatrix(dim1);
	  V12 = createSquareMatrix(dim1);
	  V21 = createSquareMatrix(dim1);
	  g00 = createSquareMatrix(dim1);

	  SL = createSquareMatrix(dim1);

	  //generate the info required for Rubio method
	  lead_prep(En, LeadCells[leadloop], leadloop, params, ginv, V12, V21);

	  InvertMatrixGSL(ginv, g00, dim1);
	  RubioSGF(SL, g00, V12, V21, dim1, &lcount, elemerr*dim1*dim1);
	  
	  
	    if(leadloop==0)
	  {
// 	    printf("DIM %d\n", dim1);
// 	    listNonZero(ginv, dim1, dim1);
// 	    listNonZero(V12, dim1, dim1);
// 	    listNonZero(V21, dim1, dim1);
	  }
	  FreeMatrix(ginv); FreeMatrix(V12); FreeMatrix(V21); FreeMatrix (g00);
	

  
  //connections of leads to device
	  //leads are connected to the device using the hopping rules of the lead region, 
	  //not the device region (if different)
	  
	  dim1a = (cellinfo->lead_dims)[leadloop];
// 	            printf("#lead %d    dim %d %d\n",leadloop, dim1, dim1a);

	  
	  VLD = createNonSquareMatrix(dim1, dim1a);
	  VDL = createNonSquareMatrix(dim1a, dim1);
	  
	  for(i=0; i <dim1; i++)
	  {
	    for(j=0; j<dim1a; j++)
	    {
		k=(cellinfo->lead_sites)[dimcounta + j];
		VLD[i][j] =  (hopfn)(LeadCells[leadloop], DeviceCell, i, k, bshifts0, (params->hoppara) );
	      	VDL[j][i] =  (hopfn)(DeviceCell, LeadCells[leadloop], k, i, bshifts0, (params->hoppara) );
	      	//VDL[j][i] =  conj(VLD[i][j] );

	    }
		  
	  }
//  	  listNonZero(VLD, dim1, dim1a);
//    	  listNonZero(VDL, dim1a, dim1);

	  temp1 = createNonSquareMatrix(dim1a, dim1);
	  MatrixMultNS(VDL, SL, temp1, dim1a, dim1, dim1);
	  smallSigma = createSquareMatrix(dim1a);

	  MatrixMultNS(temp1, VLD, smallSigma, dim1a, dim1, dim1a);

	  FreeMatrix(temp1); FreeMatrix(VLD); FreeMatrix(VDL);
	  
	  MatrixCopyPart(smallSigma, Sigma, 0, 0, dimcounta, dimcounta, dim1a, dim1a);
	  FreeMatrix(smallSigma); FreeMatrix(SL);
	  
	  dimcounta += dim1a;
	  
    }
// 	      	  listNonZero(Sigma, dimcounta, dimcounta);

    free(bshifts0);
}


//very simple lead Sigmas 
//assumes a constant density of states in metal leads
//diagonal terms in Sigma constant
//off diagonal terms with optional, variable distance decay
//Does not use LeadCells (no repeated unit cells or SGF to calculate...) or hoppingfunc
//model params are stored within params->(hoppara->hops)
// [0] is 'Sig' - the magnitude of the imaginary part of the diagonal SGF term 
// [1] is 'alpha', [2] is 'beta' in the off diagonal expression i alpha (diagonal value) / separation^beta
// [3] is the hopping between lead and device (this should allow is to be made spin dependent easier?)
void multipleSimplestMetalLeads (double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, lead_para *params, double _Complex **Sigma)
{
  
    int leadloop, dim1, dim1a,  dimcounta=0, lcount;
    double _Complex **smallSigma, **temp1;

    int num_leads = (cellinfo->num_leads);
    gen_hop_params *hopp = (params->hoppara);
    double sep;
    double _Complex *hops = (hopp->hops);
    int i, j, k;
    int iprime, jprime;
  
    
    
    for(leadloop=0; leadloop < num_leads; leadloop++)
    {
  
  
            dim1a = (cellinfo->lead_dims)[leadloop];

            
            smallSigma = createSquareMatrix(dim1a);
            
            for(i=0; i< dim1a; i++)
            {
                iprime = (cellinfo->lead_sites)[dimcounta + i];
                for(j=0; j<dim1a; j++)
                {
                    jprime = (cellinfo->lead_sites)[dimcounta + j];
                    sep = sqrt( pow((DeviceCell->pos)[jprime][0] - (DeviceCell->pos)[iprime][0], 2) + pow((DeviceCell->pos)[jprime][1] - (DeviceCell->pos)[iprime][1], 2));
                    
                    smallSigma[i][j] = I * hops[0] * hops[3] * conj(hops[3]);
                    
                    if(i!=j)
                    {
                        smallSigma[i][j] = smallSigma[i][j] * hops[1] / (pow(sep, hops[2]));
                    }
                    
                }
            }
            
            
            MatrixCopyPart(smallSigma, Sigma, 0, 0, dimcounta, dimcounta, dim1a, dim1a);
            FreeMatrix(smallSigma); 
            
            dimcounta += dim1a;
	  
    }

}




//simplest way to generate sigmas
void simple2leads (double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, lead_para *params, double _Complex **Sigma)
{
  

	int geo = (DeviceCell->geo);
	int i, j, k;
	int dim1, dim2; 
  	double elemerr=1.0e-15;
	int lcount, rcount;
	double _Complex **ginv, **V12, **V21, **g00, **SL, **SR, **VLD, **VDL, **smallSigma, **temp1;
	int dim1a, dim2a;
	double *bshifts0 = createDoubleArray(3);
	    hoppingfunc *hopfn = (hoppingfunc *)(params->hopfn);

	
	    //a cleverer, more general version of this routine could simply loop over an arbitrary # of leads
	    //not much would need to be changed for this to work
	    
//    //generate lead SGFs
	
      //lead 1
	  dim1 = *(LeadCells[0]->Nrem);
	  ginv = createSquareMatrix(dim1);
	  V12 = createSquareMatrix(dim1);
	  V21 = createSquareMatrix(dim1);
	  g00 = createSquareMatrix(dim1);

	  SL = createSquareMatrix(dim1);
      
	  //generate the info required for Rubio method
	  lead_prep(En, LeadCells[0], 0, params, ginv, V12, V21);
	  
// 	  listNonZero(ginv, dim1, dim1);
// 	  listNonZero(V12, dim1, dim1);
// 	  listNonZero(V21, dim1, dim1);

	  InvertMatrixGSL(ginv, g00, dim1);
	  RubioSGF(SL, g00, V12, V21, dim1, &lcount, elemerr*dim1*dim1);
	  FreeMatrix(ginv); FreeMatrix(V12); FreeMatrix(V21); FreeMatrix (g00);
	  
      //lead 2
	  dim2 = *(LeadCells[1]->Nrem);
	  ginv = createSquareMatrix(dim2);
	  V12 = createSquareMatrix(dim2);
	  V21 = createSquareMatrix(dim2);
	  g00 = createSquareMatrix(dim2);

	  SR = createSquareMatrix(dim2);
      
	  //generate the info required for Rubio method
	  lead_prep(En, LeadCells[1], 1, params, ginv, V12, V21);
	  
// 	  listNonZero(ginv, dim2, dim2);
// 	  listNonZero(V12, dim2, dim2);
// 	  listNonZero(V21, dim2, dim2);

	  InvertMatrixGSL(ginv, g00, dim2);
	  RubioSGF(SR, g00, V12, V21, dim2, &lcount, elemerr*dim2*dim2);
	  FreeMatrix(ginv); FreeMatrix(V12); FreeMatrix(V21); FreeMatrix (g00);

 	
	 //connections of left and right leads to device
	  //leads are connected to the device using the hopping rules of the lead region, not the device region
	  
	  //lead1
	  dim1a = (cellinfo->lead_dims)[0];
	  VLD = createNonSquareMatrix(dim1, dim1a);
	  VDL = createNonSquareMatrix(dim1a, dim1);
	  
	  for(i=0; i <dim1; i++)
	  {
	    for(j=0; j<dim1a; j++)
	    {
		k=(cellinfo->lead_sites)[j];
		VLD[i][j] =  (hopfn)(LeadCells[0], DeviceCell, i, k, bshifts0, (params->hoppara) );
	      	VDL[j][i] =  (hopfn)(DeviceCell, LeadCells[0], k, i, bshifts0, (params->hoppara) );
	      	//VDL[j][i] =  conj(VLD[i][j] );

	    }
		  
	  }
//  	  listNonZero(VLD, dim1, dim1a);
//  	  listNonZero(VDL, dim1a, dim1);

	  temp1 = createNonSquareMatrix(dim1a, dim1);
	  MatrixMultNS(VDL, SL, temp1, dim1a, dim1, dim1);
	  smallSigma = createSquareMatrix(dim1a);
	  MatrixMultNS(temp1, VLD, smallSigma, dim1a, dim1, dim1a);
	  FreeMatrix(temp1); FreeMatrix(VLD); FreeMatrix(VDL);
	  
	  MatrixCopyPart(smallSigma, Sigma, 0, 0, 0, 0, dim1a, dim1a);
	  FreeMatrix(smallSigma); FreeMatrix(SL);


	  //lead2
	  dim2a = (cellinfo->lead_dims)[1];
	  VLD = createNonSquareMatrix(dim2, dim2a);
	  VDL = createNonSquareMatrix(dim2a, dim2);
	  
	  for(i=0; i <dim2; i++)
	  {
	    for(j=0; j<dim2a; j++)
	    {
		k=(cellinfo->lead_sites)[dim1a+j];
		VLD[i][j] =  (hopfn)(LeadCells[1], DeviceCell, i, k, bshifts0, (params->hoppara) );
	      	VDL[j][i] =  (hopfn)(DeviceCell, LeadCells[1], k, i, bshifts0, (params->hoppara) );
// 	      	VDL[j][i] =  conj(VLD[i][j] );

	    }
		  
	  }
// 	  listNonZero(VLD, dim2, dim2a);
//   	  listNonZero(VDL, dim2a, dim2);
	  
	  temp1 = createNonSquareMatrix(dim2a, dim2);
	  MatrixMultNS(VDL, SR, temp1, dim2a, dim2, dim2);
	  smallSigma = createSquareMatrix(dim2a);
	  MatrixMultNS(temp1, VLD, smallSigma, dim2a, dim2, dim2a);
	  FreeMatrix(temp1); FreeMatrix(VLD); FreeMatrix(VDL);
	  
	  MatrixCopyPart(smallSigma, Sigma, 0, 0, dim1a, dim1a, dim2a, dim2a);
	  FreeMatrix(smallSigma); FreeMatrix(SR);
	  free(bshifts0);

}


//general routine to prepare lead cells for Rubio-esque routine
//takes positions, cell sep. vector and hopping rule
//returns unit cell ginv, V12 and V21
void lead_prep(double _Complex En, RectRedux *LeadCell, int leadindex, lead_para *params, double _Complex **ginv, double _Complex **V12, double _Complex **V21)
{
    int i, j, k;
    
    int dim = *(LeadCell->Nrem);
    
    hoppingfunc *hopfn = (hoppingfunc *)(params->hopfn);
    double *bshifts = createDoubleArray(3);
    
  
    //g00
    for(i=0; i <dim; i++)
    {
      ginv[i][i] = En - (LeadCell->site_pots[i]) - (hopfn)(LeadCell, LeadCell, i, i, bshifts, (params->hoppara) ) ;
      
        //magnitude is negative, hence += to add to the inverse GF
          if(LeadCell->cap_pots != NULL)
          {
              ginv[i][i] += I* (LeadCell->cap_pots)[i];
          }
      
      
      for(j=0; j<dim; j++)
      {
	if(j!=i)
	{
	  ginv[i][j] = - (hopfn)(LeadCell, LeadCell, i, j, bshifts, (params->hoppara) );
	}
	
      }
      
            
    }
   
//   if(leadindex ==0)
//   {
//      printf("GINV %d\n", dim);
// 	listNonZero(ginv, dim, dim);
//   }

      //Vs
    for(i=0; i<3; i++)
      bshifts[i] = (params->shift_vecs)[leadindex][i];
    
    
      for(i=0; i <dim; i++)
      {
	for(j=0; j<dim; j++)
	{
	  
	    V12[i][j] =  (hopfn)(LeadCell, LeadCell, i, j, bshifts, (params->hoppara) );
	  
	  
	}
	      
      }

      
    for(i=0; i<3; i++)
      bshifts[i] = -(params->shift_vecs)[leadindex][i];
    
    
      for(i=0; i <dim; i++)
      {
	for(j=0; j<dim; j++)
	{
	  
	    V21[i][j] =  (hopfn)(LeadCell, LeadCell, i, j, bshifts, (params->hoppara) );
// 	    V21[i][j] =  conj(V12[j][i]);

	}
	      
      }
    free(bshifts);

}


//general routine to prepare lead cells for Rubio-esque routine
//takes positions, cell sep. vector and hopping rule
//returns unit cell ginv, V12 and V21
void lead_prep2(double _Complex En, RectRedux *LeadCell, int leadindex, rib_lead_para *params, double _Complex **ginv, double _Complex **V12, double _Complex **V21)
{
    int i, j, k;
    
    int dim = *(LeadCell->Nrem);
    
    hoppingfunc *hopfn = (hoppingfunc *)(params->hopfn);
    double *bshifts = createDoubleArray(3);
    
  
    //g00
    for(i=0; i <dim; i++)
    {
      ginv[i][i] = En - (LeadCell->site_pots[i]) - (hopfn)(LeadCell, LeadCell, i, i, bshifts, (params->hoppara) ) ;
      
        //magnitude is negative, hence += to add to the inverse GF
          if(LeadCell->cap_pots != NULL)
          {
              ginv[i][i] += I* (LeadCell->cap_pots)[i];
          }
      
      for(j=0; j<dim; j++)
      {
	if(j!=i)
	{
	  ginv[i][j] = - (hopfn)(LeadCell, LeadCell, i, j, bshifts, (params->hoppara) );
	}
	
      }
      
            
    }
   
//   if(leadindex ==0)
//   {
//      printf("GINV %d\n", dim);
// 	listNonZero(ginv, dim, dim);
//   }

      //Vs
    for(i=0; i<3; i++)
      bshifts[i] = (params->shift_vec)[i];
    
    
      for(i=0; i <dim; i++)
      {
	for(j=0; j<dim; j++)
	{
	  
	    V12[i][j] =  (hopfn)(LeadCell, LeadCell, i, j, bshifts, (params->hoppara) );
	  
	  
	}
	      
      }

      
    for(i=0; i<3; i++)
      bshifts[i] = -(params->shift_vec)[i];
    
    
      for(i=0; i <dim; i++)
      {
	for(j=0; j<dim; j++)
	{
	  
	    V21[i][j] =  (hopfn)(LeadCell, LeadCell, i, j, bshifts, (params->hoppara) );
// 	    V21[i][j] =  conj(V12[j][i]);

	}
	      
      }
    free(bshifts);

}
  
  
//gate induced potential - variations on the efetov model.
//places a potential on sites depending on gate and geometry parameters

void gate_induced_pot ( int vgtype, RectRedux *DeviceCell, double *engdeppots, double gate_voltage, double edge_cut_off, double subs_thick, double subs_epsr)
{
  int length1 = (DeviceCell->length);
  int length2 = (DeviceCell->length2);
  int geo = (DeviceCell->geo);
  int Nrem = *(DeviceCell->Nrem);
  double Aconst = 1.579, Bconst = 0.0001887, Cconst=0.173;
  double **pos = (DeviceCell->pos);
  
  double top_edge, bottom_edge, rib_center, ypos, ny, nmax, n0, zval, zmax;
  int i;
  
  //top and bottom edges of ribbon for potential calculations
      //zigzag naturally offset by definition of atomic positions (none along y=0)
	if(geo==0)
	{
	    top_edge=length1*sqrt(3)/2;
	    bottom_edge=0.0;
	    rib_center = top_edge/2;
	}
	
	  //armchair offset slightly to avoid division by zero for very edge sites
	if(geo==1 )
	{
	    top_edge=(length1)*0.5;
	    bottom_edge = -0.5;
	    rib_center = top_edge/2 -0.25;
	}
	
	
	//density induced by this setup in an infinite width graphene sheet (SI units, m^-2)
	//note the 2 in the denominator - this assumes a spinless calculation.
	n0= (subs_epsr)*(8.85E-12)*fabs(gate_voltage)/(2*subs_thick * 1.6E-19);
	nmax = n0 * ((top_edge-bottom_edge)/2) / sqrt( pow( (top_edge-bottom_edge)/2, 2)  - pow((top_edge-rib_center) - edge_cut_off, 2) );
  
	for(i=0; i<Nrem; i++)
	{
	  ypos = pos[i][1] - rib_center;
	  
	  if (vgtype==0)
	  {
	    ny = n0;
	  }
	  
	  if (vgtype==1)
	  {
	    ny = n0 * ((top_edge-bottom_edge)/2)  / sqrt( pow((top_edge-bottom_edge)/2, 2)  - pow(ypos, 2) );
	
	    if(ny>nmax)
	      ny=nmax;
	  }
	  
	  if (vgtype==3)
	  {
		zval = 2*ypos / (top_edge-bottom_edge);
		zmax = 2*((top_edge-rib_center) - edge_cut_off) / (top_edge-bottom_edge);
		ny = n0 * (Aconst * zval*zval + Bconst*zval + Cconst) / sqrt(1 - zval*zval) ;
		nmax = n0 * (Aconst * zmax*zmax + Bconst*zmax + Cconst) / sqrt(1 - zmax*zmax) ;
		
		
		
	   // ny = n0 * ((top_edge-bottom_edge)/2)  / sqrt( pow((top_edge-bottom_edge)/2, 2)  - pow(ypos, 2) );
	
	    if(ny>nmax)
	      ny=nmax;
	  }
	  
	  
	  engdeppots[i] = (- gate_voltage/ fabs(gate_voltage)) * (1.05E-28) * (ny/fabs(ny))* sqrt(M_PI * fabs(ny) ) / (2.7 * 1.6E-19) ;
	  
	  if(gate_voltage==0.0)
	  {
	    engdeppots[i] = 0.0;
	  }
	}
  
}



void multipleCustomLeads (double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, lead_para *params, double _Complex **Sigma)
{
  
    int leadloop, dim1, dim1a,  dimcounta=0, lcount;
    double _Complex **smallSigma, **temp1;

    int num_leads = (cellinfo->num_leads);
    gen_hop_params *hopp = (params->hoppara);
    double sep;
    double _Complex *hops = (hopp->hops);
    int i, j, k;
    int iprime, jprime;
    
    multiple_para *multiple;
	singleleadfunction *singlefn;
    void * singleparams;
    
    for(leadloop=0; leadloop < num_leads; leadloop++)
    {
  
		multiple = (multiple_para *)(params-> multiple)[leadloop];
		singlefn = (singleleadfunction *)(multiple -> indiv_lead_fn);
		singleparams = (multiple -> indiv_lead_para);
		
            dim1a = (cellinfo->lead_dims)[leadloop];

            
            smallSigma = createSquareMatrix(dim1a);
            
           //external function call of multiple.indiv_lead_fn to calculate sigma
	    (singlefn)(leadloop, En, DeviceCell, LeadCells, cellinfo, singleparams, smallSigma);
            
            
            MatrixCopyPart(smallSigma, Sigma, 0, 0, dimcounta, dimcounta, dim1a, dim1a);
            FreeMatrix(smallSigma); 
            
            dimcounta += dim1a;
	  
    }

}


//generate Sigma for a single Ribbon type lead
void singleRibbonLead (int leadnum, double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, void *params, double _Complex **Sigma)
{
    int leadloop, dim1, dim1a,  lcount, dimcounta=0;
    double _Complex **ginv, **V12, **V21, **g00, **SL, **SR, **VLD, **VDL, **smallSigma, **temp1;
    double elemerr=1.0e-15;
	
    rib_lead_para *ribpara = (rib_lead_para *)params;
    
    int i, j, k;
  
    double *bshifts0 = createDoubleArray(3);
    hoppingfunc *hopfn = (hoppingfunc *)(ribpara->hopfn);
    
	for(leadloop=0; leadloop < leadnum; leadloop++)
	{
		dim1a = (cellinfo->lead_dims)[leadloop];
		dimcounta += dim1a;
	}
      //generate leads SGFs
  
	  dim1 = *(LeadCells[leadnum]->Nrem);
	  ginv = createSquareMatrix(dim1);
	  V12 = createSquareMatrix(dim1);
	  V21 = createSquareMatrix(dim1);
	  g00 = createSquareMatrix(dim1);

	  SL = createSquareMatrix(dim1);

	  //generate the info required for Rubio method
	  lead_prep2(En, LeadCells[leadnum], leadnum, ribpara, ginv, V12, V21);

	  InvertMatrixGSL(ginv, g00, dim1);
	  RubioSGF(SL, g00, V12, V21, dim1, &lcount, elemerr*dim1*dim1);
	  
	  

	  FreeMatrix(ginv); FreeMatrix(V12); FreeMatrix(V21); FreeMatrix (g00);
	

    //connections of leads to device
	  //leads are connected to the device using the hopping rules of the lead region, 
	  //not the device region (if different)
	  
	  dim1a = (cellinfo->lead_dims)[leadnum];

	  VLD = createNonSquareMatrix(dim1, dim1a);
	  VDL = createNonSquareMatrix(dim1a, dim1);
	  
	  for(i=0; i <dim1; i++)
	  {
	    for(j=0; j<dim1a; j++)
	    {
		k=(cellinfo->lead_sites)[dimcounta + j];
		VLD[i][j] =  (hopfn)(LeadCells[leadnum], DeviceCell, i, k, bshifts0, (ribpara->hoppara) );
	      	VDL[j][i] =  (hopfn)(DeviceCell, LeadCells[leadnum], k, i, bshifts0, (ribpara->hoppara) );
	    }
		  
	  }


	  temp1 = createNonSquareMatrix(dim1a, dim1);
	  MatrixMultNS(VDL, SL, temp1, dim1a, dim1, dim1);
	  MatrixMultNS(temp1, VLD, Sigma, dim1a, dim1, dim1a);
	  FreeMatrix(temp1); FreeMatrix(VLD); FreeMatrix(VDL);
	  FreeMatrix(SL);
	
    free(bshifts0);
}




//very simple metal lead Sigma 
//assumes a constant density of states in metal leads
//diagonal terms in Sigma constant
//off diagonal terms with optional, variable distance decay
//Does not use LeadCells (no repeated unit cells or SGF to calculate...) or hoppingfunc
//model params are stored within params->(hoppara->hops)
// [0] is 'Sig' - the magnitude of the imaginary part of the diagonal SGF term 
// [1] is 'alpha', [2] is 'beta' in the off diagonal expression i alpha (diagonal value) / separation^beta
// [3] is the hopping between lead and device (this should allow is to be made spin dependent easier?)
void singleSimplestMetalLead (int leadnum, double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, void *params, double _Complex **Sigma)
{
  
	metal_lead_para *metalpara = (metal_lead_para *)params;
	
	
	
    int leadloop, dim1, dim1a,  dimcounta=0, lcount;
    double _Complex **temp1;

    int num_leads = (cellinfo->num_leads);
    gen_hop_params *hopp = (metalpara->hoppara);
    double sep;
    double _Complex *hops = (hopp->hops);
    int i, j, k;
    int iprime, jprime;
  
    
    
    for(leadloop=0; leadloop < leadnum; leadloop++)
	{
		dim1a = (cellinfo->lead_dims)[leadloop];
		dimcounta += dim1a;
	}
  
  //generate leads SGFs
  
	  dim1 = (cellinfo->lead_dims)[leadnum];
  
            
            for(i=0; i< dim1; i++)
            {
                iprime = (cellinfo->lead_sites)[dimcounta + i];
                for(j=0; j<dim1; j++)
                {
                    jprime = (cellinfo->lead_sites)[dimcounta + j];
                    sep = sqrt( pow((DeviceCell->pos)[jprime][0] - (DeviceCell->pos)[iprime][0], 2) + pow((DeviceCell->pos)[jprime][1] - (DeviceCell->pos)[iprime][1], 2));
                    
                    Sigma[i][j] = I * hops[0] * hops[3] * conj(hops[3]);
		    if(cimag(En)<0)
		    {
			    Sigma[i][j] = - Sigma[i][j];
		    }
		    
                    
                    if(i!=j && sep> 0)
                    {
                        Sigma[i][j] = Sigma[i][j] * hops[1] / (pow(sep, hops[2]));
                    }
                    if(i!=j && sep== 0)
                    {
                        Sigma[i][j] = 0.0;
                    }
                    
                    
                }
            }
            
            
            
	  
    
	    

}





