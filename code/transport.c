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


      //calculate advanced quantities
	 if( (tpara->TRsym) == 0)
	 {
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
	   
	   
	 }
	 else if( (tpara->TRsym) == 1)
	 {
	   
	   if(mode==1)
	   {
	      devicemode2 = 2;
	   }
	   
	    (leadfn)(creal(En) - I*cimag(En), DeviceCell, Leads, cellinfo, leadsparams, SigmaA);
	    genDeviceGF(creal(En) - I*cimag(En), DeviceCell, cnxp, cellinfo, hoppingfn, hoppingparams, devicemode, devicemode2, g_sys_a, NULL, g1i, SigmaA);
	 }
	 
	 
	 //Gamma
	  for(i=0; i<cell1dim; i++)
	  {
	    for(j=0; j<cell1dim; j++)
	    {
	      Gamma[i][j] = I*(SigmaR[i][j] - SigmaA[i][j]);
	    }
	  }
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
	  
	  g00inv[i][i] = En - (DeviceCell->site_pots)[index1] - hoppingfn(DeviceCell, DeviceCell, index1, index1, bshifts, hoppingparams);;
	  
	  
	  
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
// 	  	printf("#g00inv\n");
// 	  	listNonZero(g00inv, dim, dim);
// 	  
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







//this has been generalised relative to the antidot code version (RectRedux param set)
double genConduc5(double _Complex En, RectRedux *DeviceCell, double hopping)
{
	  //important definitions and quantities
	    int length1 = (DeviceCell->length);
	    int length2 = (DeviceCell->length2);
	    int geo = (DeviceCell->geo);
//	    int vcells = botcells + midcells + topcells;

	    int Nrem = *(DeviceCell->Nrem);
	    int **chaininfo = (DeviceCell->chaininfo);
	    int **siteinfo = (DeviceCell->siteinfo);
	    double *site_pots = (DeviceCell->site_pots);

	    int i, j, k, l, foundcleanyet, theclean, lead_cell_dim; //dimV1=0, dimV2=0;
	    int cell, kstart, lstart;
	    
	   //generate unit cell for the leads and a clean cell GF
		double _Complex **G_lead_cell;
		double _Complex **VL1temp, **VL1, **V1L, **SL, **SR, **VLR, **VRL, **ginv, **g00, **ginv2;
		double _Complex **SLred, **SRred;
		

		    
	  //generate left and right lead SGFs
	      double elemerr=1.0e-15;
	      int lcount, rcount;
	      

	      SLred = createSquareMatrix(2*length1);
	      SRred = createSquareMatrix(2*length1);
	      VL1temp = createSquareMatrix(2*length1);
 	      VRL = createSquareMatrix(2*length1);
	      ginv = createSquareMatrix(2*length1);
	      g00 = createSquareMatrix(2*length1);
	      
	      
	    if(geo==0)
	    {
	      ZGNR (ginv, VL1temp, VRL, length1, hopping, En);
	    }
	    if(geo==1)
	    {
	      AGNR (ginv, VL1temp, VRL, length1, hopping, En);
	    }
	    InvertMatrixGSL(ginv, g00, 2*length1);

	      RubioSGF(SLred, g00, VRL, VL1temp, 2*length1, &lcount, elemerr*length1*length1);
	      RubioSGF(SRred, g00, VL1temp, VRL, 2*length1, &rcount, elemerr*length1*length1);
	
		  FreeMatrix(VRL); 
		  //FreeMatrix(ginv);
		  FreeMatrix(g00);
		  
		 // printf("%lf	%e	%e\n", creal(En), creal(SRred[0][0]), cimag(SRred[0][0]));

		  
		int dim_old, dim_new, l_count, chain;  
		double _Complex **old_SGF, **new_SGF, **Sigma_L, **Sigma_R, **t1, **t2, **t3;

		
	    //left_lead Sigma
		dim_old = chaininfo[0][0];
		kstart=0; lstart=0; l_count=0;
		VL1 = createNonSquareMatrix(2*length1, dim_old);
		    
		Sigma_L = createSquareMatrix(dim_old);
		MatrixCopy(VL1temp, VL1, 2*length1);


		V1L = createNonSquareMatrix(dim_old, 2*length1);
		for(i=0; i<dim_old; i++)
		{
		  for(j=0; j<2*length1; j++)
		  {
		    V1L[i][j] = VL1[j][i];
		  }
		}
		
		t1=createNonSquareMatrix(dim_old, 2*length1);
		MatrixMultNS(V1L, SLred, t1, dim_old, 2*length1, 2*length1);
		MatrixMultNS(t1, VL1, Sigma_L, dim_old, 2*length1, dim_old);
		FreeMatrix(t1); FreeMatrix(V1L); FreeMatrix(VL1); FreeMatrix(SLred);    
		
		//gamma_L
		    double _Complex **Gamma_L = createSquareMatrix(dim_old);
		    
		    for(i=0; i<dim_old;i++)
		    {
		      for(j=0; j<dim_old; j++)
		      {
			Gamma_L[i][j] = -2*cimag(Sigma_L[i][j]);
		      }
		    }
		    
		double _Complex **Gamma_R = createSquareMatrix(dim_old);

		  
	      //right_lead Sigma
		dim_old = chaininfo[length2-1][0];
		kstart=0; lstart=0; l_count=0;
		VL1 = createNonSquareMatrix(2*length1, dim_old);
		Sigma_R = createSquareMatrix(dim_old);
		    
		kstart=0; lstart=0; l_count=0;


		for(i=0; i<2*length1; i++)
		{
		  for(j=0; j<dim_old; j++)
		  {
		    VL1[i][j] = VL1temp[j][i];
		  }
		}

		V1L = createNonSquareMatrix(dim_old, 2*length1);
		for(i=0; i<dim_old; i++)
		{
		  for(j=0; j<2*length1; j++)
		  {
		    V1L[i][j] = VL1[j][i];
		  }
		}
		
		t1=createNonSquareMatrix(dim_old, 2*length1);
		MatrixMultNS(V1L, SRred, t1, dim_old, 2*length1, 2*length1);
		MatrixMultNS(t1, VL1, Sigma_R, dim_old, 2*length1, dim_old);
		FreeMatrix(t1); FreeMatrix(V1L); FreeMatrix(VL1); FreeMatrix(SRred);    


		
		//folding back from the right to g11^R
		  double _Complex **Sigma  = createSquareMatrix(dim_old);

		  double _Complex **VLR_temp, **gtemp;
		  dim_old = 2*length1; //old dimension is dimension of RHS lead chain
		  double _Complex **S_temp = createSquareMatrix(dim_old);

		 /* 
		  VLR_temp = createSquareMatrix(2*length1);
		    //unit cell left -> right connections
		      for(i=0; i < length1; i++)
		      {
			  VLR_temp[4*i][4*i + 1] = hopping;
			  VLR_temp[4*i+3][4*i+2] = hopping;
		      }*/
		  
		  for(chain=length2-1; chain >=0; chain--)
		  {

		    dim_new = chaininfo[chain][0];
		    ginv2 = createSquareMatrix(dim_new);
		    FreeMatrix(Sigma);
		    Sigma = createSquareMatrix(dim_new);
		    
		    //bare bones H
			k=0;
			for(i=0; i<2*length1; i++)
			{
			    if(siteinfo[(chain)*2*length1 + i][0] == 0)
			    {
				ginv2[k][k] = En - site_pots[(chain)*2*length1 + i];
				
				for(j=0; j<2*length1; j++)
				{
				  if(i!=j)
				  {
				    if(siteinfo[(chain)*2*length1 + j][0] == 0)
				    {
				      ginv2[k][j] = ginv[k][j];
				    }
				  }
				}
				
			
				k++;
			    }
			}
			
		    //self energy terms
		    
			  //right-most chain has RHS lead self energy term
			  if(chain == length2 - 1)
			    MatrixCopy(Sigma_R, Sigma, dim_new);
			  
			  //other chains need their self energy terms constructed from existing GFs and generated V matrices
			  if(chain < length2 - 1)
			  {
			      
			      VLR = createNonSquareMatrix(dim_new, dim_old);
			      VRL = createNonSquareMatrix(dim_old, dim_new);
			      
			      k=0;
			      for(i=0; i<2*length1; i++)
			      {
				  if(siteinfo[(chain)*2*length1 + i][0] == 0)
				  {
				      l=0;
				      for(j=0; j< 2*length1; j++)
				      {
					  if(siteinfo[(chain+1)*2*length1 + j][0] == 0)
					  {
					      
					    VLR[k][l] =  VL1temp[i][j];
					    VRL[l][k] =  VL1temp[i][j];
					    l++;
					  }
				      }
				    k++;
				  }
			      }
			      
			      		      //printf("use Stemp: dim %d\n", dim_old);

			      t1=createNonSquareMatrix(dim_new, dim_old);
			      MatrixMultNS(VLR, S_temp, t1, dim_new, dim_old, dim_old);
			      MatrixMultNS(t1, VRL, Sigma, dim_new, dim_old, dim_new);	
			      FreeMatrix(t1); FreeMatrix(VLR); FreeMatrix(VRL);
			    
			  }
			  

			  //left most chain also has SE contribution from left lead
			  if(chain == 0)
			  {
			    t1 = createSquareMatrix(dim_new);
			    
			    //Gamma_R
			    for(i=0; i<dim_new;i++)
			    {
			      for(j=0; j<dim_new; j++)
			      {
				Gamma_R[i][j] = -2*cimag(Sigma[i][j]);
			      }
			    }	
			    
			    
			    MatrixAdd(Sigma, Sigma_L, t1, dim_new);
			    MatrixCopy(t1, Sigma, dim_new);
			    FreeMatrix(t1);
			  }
			  
			  
		  //Add self energies to bare hamiltonian
		      t1=createSquareMatrix(dim_new);
		      gtemp = createSquareMatrix(dim_new);
		      //printf("ok1\n");

		      MatrixSubtract(ginv2, Sigma, t1, dim_new);
		     // printEMatrix(t1, dim_new);
		      InvertMatrixGSL(t1, gtemp, dim_new);
		      FreeMatrix(t1); FreeMatrix(ginv2); 
		      				    		    //printf("ok2\n");

		      FreeMatrix(S_temp);
		      S_temp=createSquareMatrix(dim_new);
		      MatrixCopyPart(gtemp, S_temp, 0, 0, 0, 0, chaininfo[chain][0], chaininfo[chain][0]);
		      FreeMatrix(gtemp);
		      
		      dim_old = dim_new;

		  }  
		  
		  
		  double _Complex **Sadv = createSquareMatrix(dim_new);
		  t1=createSquareMatrix(dim_new);
		  t2=createSquareMatrix(dim_new);

		  for(i=0; i<dim_new; i++)
		  {
		    for(j=0; j<dim_new; j++)
		    {
		      Sadv[i][j] = conj(S_temp[j][i]);
		    }
		  }
		  MatrixMult(Gamma_R, S_temp, t1, dim_new);
		  MatrixMult(t1, Gamma_L, t2, dim_new);
		  MatrixMult(t2, Sadv, t1, dim_new);
		  
		  
		  double cond = creal(MatrixTrace(t1, dim_new));
		  FreeMatrix(S_temp); FreeMatrix(VL1temp);  
		  FreeMatrix(Sigma); FreeMatrix(ginv);
		  FreeMatrix(t1); FreeMatrix(t2);
		  FreeMatrix(Sadv);
		  FreeMatrix(Sigma_L); FreeMatrix(Sigma_R);
		 // FreeMatrix(VLR_temp);
		  FreeMatrix(Gamma_L);
		  FreeMatrix(Gamma_R);
		  printf("%lf	%e\n", creal(En), cond);

		  return cond;

		  
}




double genConduc4(double _Complex En, 
		 RectRedux *DeviceCell,								//device params
		 double *ldoses, double **conds, int *indices, int **neigh, double hopping)								//extra output params
{
  
	//important definitions and quantities
	    int length1 = (DeviceCell->length);
	    int length2 = (DeviceCell->length2);
	    int geo = (DeviceCell->geo);
	    //int vcells = botcells + midcells + topcells;
	    
	    int Nrem = *(DeviceCell->Nrem);
	    int **chaininfo = (DeviceCell->chaininfo);
	    int **siteinfo = (DeviceCell->siteinfo);
	    double *site_pots = (DeviceCell->site_pots);
	    double **pos = (DeviceCell->pos);


	    int i, j, k, l, foundcleanyet, theclean, lead_cell_dim; //dimV1=0, dimV2=0;
	    int cell, kstart, lstart;
	    
	   
// 	   for(i=0; i<Nrem;i++)
// 	   {
// 	     printf("%lf\n", site_pots[i]);
// 	   }
// 	   
	   
	//generate unit cell for the leads and a clean cell GF
		double _Complex **G_lead_cell;
		double _Complex **VL1temp, **VL1, **V1L, **SL, **SR, **VLR, **VRL, **ginv, **g00, **ginv2;
		double _Complex **SLred, **SRred;
		
	
		    
		    
	//generate left and right lead SGFs
	      double elemerr=1.0e-15;
	      int lcount, rcount;
	      
	
	      SLred = createSquareMatrix(2*length1);
	      SRred = createSquareMatrix(2*length1);
	      VL1temp = createSquareMatrix(2*length1);
 	      VRL = createSquareMatrix(2*length1);
	      ginv = createSquareMatrix(2*length1);
	      g00 = createSquareMatrix(2*length1);
	      
	      
	    if(geo==0)
	    {
	      ZGNR (ginv, VL1temp, VRL, length1, hopping, En);
	    }
	    if(geo==1)
	    {
	      AGNR (ginv, VL1temp, VRL, length1, hopping, En);
	    }
	    InvertMatrixGSL(ginv, g00, 2*length1);

	      RubioSGF(SLred, g00, VRL, VL1temp, 2*length1, &lcount, elemerr*length1*length1);
	      RubioSGF(SRred, g00, VL1temp, VRL, 2*length1, &rcount, elemerr*length1*length1);
	
		  FreeMatrix(VRL); 
		  //FreeMatrix(ginv);
		  FreeMatrix(g00);
		  
		  
		int dim_old, dim_new, l_count, chain;  
		double _Complex **old_SGF, **new_SGF, **Sigma_L, **Sigma_R, **t1, **t2, **t3;

		
	    //left_lead Sigma
		dim_old = chaininfo[0][0];
		kstart=0; lstart=0; l_count=0;
		VL1 = createNonSquareMatrix(2*length1, dim_old);
		    
		Sigma_L = createSquareMatrix(dim_old);
		MatrixCopy(VL1temp, VL1, 2*length1);


		V1L = createNonSquareMatrix(dim_old, 2*length1);
		for(i=0; i<dim_old; i++)
		{
		  for(j=0; j<2*length1; j++)
		  {
		    V1L[i][j] = VL1[j][i];
		  }
		}
		
		t1=createNonSquareMatrix(dim_old, 2*length1);
		MatrixMultNS(V1L, SLred, t1, dim_old, 2*length1, 2*length1);
		MatrixMultNS(t1, VL1, Sigma_L, dim_old, 2*length1, dim_old);
		FreeMatrix(t1); FreeMatrix(V1L); FreeMatrix(VL1); FreeMatrix(SLred);    
		
		//gamma_L
		    double _Complex **Gamma_L = createSquareMatrix(dim_old);
		    
		    for(i=0; i<dim_old;i++)
		    {
		      for(j=0; j<dim_old; j++)
		      {
			Gamma_L[i][j] = -2*cimag(Sigma_L[i][j]);
		      }
		    }
		    
		double _Complex **Gamma_R = createSquareMatrix(dim_old);

		  
	      //right_lead Sigma
		dim_old = chaininfo[length2-1][0];
		kstart=0; lstart=0; l_count=0;
		VL1 = createNonSquareMatrix(2*length1, dim_old);
		Sigma_R = createSquareMatrix(dim_old);
		    
		kstart=0; lstart=0; l_count=0;


		for(i=0; i<2*length1; i++)
		{
		  for(j=0; j<dim_old; j++)
		  {
		    VL1[i][j] = VL1temp[j][i];
		  }
		}

		V1L = createNonSquareMatrix(dim_old, 2*length1);
		for(i=0; i<dim_old; i++)
		{
		  for(j=0; j<2*length1; j++)
		  {
		    V1L[i][j] = VL1[j][i];
		  }
		}
		
		t1=createNonSquareMatrix(dim_old, 2*length1);
		MatrixMultNS(V1L, SRred, t1, dim_old, 2*length1, 2*length1);
		MatrixMultNS(t1, VL1, Sigma_R, dim_old, 2*length1, dim_old);
		FreeMatrix(t1); FreeMatrix(V1L); FreeMatrix(VL1); FreeMatrix(SRred);  
		
		
		
		
		  
			  
	//folding back from the right to g11^R
	    double _Complex **allgs = createNonSquareMatrix(Nrem, 2*length1);
	    double _Complex **Sigma, **S_temp, **VLR_temp, **gtemp;
	    dim_old = 2*length1; //old dimension is dimension of RHS lead chain
	    
	  
	    for(chain=length2-1; chain >=0; chain--)
	    {
	      //printf("# folding back %d of %d\n", chain, length2);
	      dim_new = chaininfo[chain][0];
	      ginv2 = createSquareMatrix(dim_new);
	      Sigma = createSquareMatrix(dim_new);
	      
	      //bare bones H
		  k=0;
			for(i=0; i<2*length1; i++)
			{
			    if(siteinfo[(chain)*2*length1 + i][0] == 0)
			    {
				ginv2[k][k] = En - site_pots[(chain)*2*length1 + i];
				//printf("%lf\n", site_pots[(chain)*2*length1 + i]);
				for(j=0; j<2*length1; j++)
				{
				  if(i!=j)
				  {
				    if(siteinfo[(chain)*2*length1 + j][0] == 0)
				    {
				      ginv2[k][j] = ginv[k][j];
				    }
				  }
				}
				
			
				k++;
			    }
			}
		  
	      //self energy terms
	      
		    //right-most chain has RHS lead self energy term
		    if(chain == length2 - 1)
		      MatrixCopy(Sigma_R, Sigma, dim_new);
		    
		    
		  // S_temp = createSquareMatrix(dim_old);

		    //other chains need their self energy terms constructed from existing GFs and generated V matrices
		    if(chain < length2 - 1)
		    {
			//FreeMatrix(S_temp);
			S_temp = createSquareMatrix(dim_old);
			MatrixCopyPart(allgs, S_temp, chaininfo[chain+1][1], 0, 0, 0,  chaininfo[chain+1][0], chaininfo[chain+1][0]);
			
			VLR = createNonSquareMatrix(dim_new, dim_old);
			VRL = createNonSquareMatrix(dim_old, dim_new);
			
			k=0;
			      for(i=0; i<2*length1; i++)
			      {
				  if(siteinfo[(chain)*2*length1 + i][0] == 0)
				  {
				      l=0;
				      for(j=0; j< 2*length1; j++)
				      {
					  if(siteinfo[(chain+1)*2*length1 + j][0] == 0)
					  {
					      
					    VLR[k][l] =  VL1temp[i][j];
					    VRL[l][k] =  VL1temp[i][j];
					    l++;
					  }
				      }
				    k++;
				  }
			      }
			
			
			t1=createNonSquareMatrix(dim_new, dim_old);
			MatrixMultNS(VLR, S_temp, t1, dim_new, dim_old, dim_old);
			MatrixMultNS(t1, VRL, Sigma, dim_new, dim_old, dim_new);	
			FreeMatrix(t1);  FreeMatrix(VLR); FreeMatrix(VRL); FreeMatrix(S_temp);
		      
		    }
		    
		    
		    //left most chain also has SE contribution from left lead
		    if(chain == 0)
		    {
		      t1 = createSquareMatrix(dim_new);
		      
		      //Gamma_R
			    for(i=0; i<dim_new;i++)
			    {
			      for(j=0; j<dim_new; j++)
			      {
				Gamma_R[i][j] = -2*cimag(Sigma[i][j]);
			      }
			    }	
		      
		      
		      MatrixAdd(Sigma, Sigma_L, t1, dim_new);
		      MatrixCopy(t1, Sigma, dim_new);
		      FreeMatrix(t1);
		    }
		    
		     
	    //Add self energies to bare hamiltonian
		t1=createSquareMatrix(dim_new);
		gtemp = createSquareMatrix(dim_new);

		MatrixSubtract(ginv2, Sigma, t1, dim_new);
		InvertMatrixGSL(t1, gtemp, dim_new);
		FreeMatrix(t1); FreeMatrix(ginv2); FreeMatrix(Sigma); 
		
		MatrixCopyPart(gtemp, allgs, 0, 0, chaininfo[chain][1], 0, chaininfo[chain][0], chaininfo[chain][0]);
		FreeMatrix(gtemp);
		
		dim_old = dim_new;
	    }  
	    
		FreeMatrix(Sigma_L); FreeMatrix(Sigma_R); 
		//printf("%lf	%lf	%lf\n", creal(En), creal(allgs[0][0]), cimag(allgs[0][0]));
		    double _Complex **Sadv = createSquareMatrix(dim_new);
		    S_temp = createSquareMatrix(dim_old);
		    MatrixCopyPart(allgs, S_temp, 0, 0, 0, 0, chaininfo[0][0], chaininfo[0][0]);
		    
		     t1=createSquareMatrix(dim_new);
		     t2=createSquareMatrix(dim_new);

		      for(i=0; i<dim_new; i++)
		      {
			for(j=0; j<dim_new; j++)
			{
			  Sadv[i][j] = conj(S_temp[j][i]);
			}
		      }
		      MatrixMult(Gamma_R, S_temp, t1, dim_new);
		      MatrixMult(t1, Gamma_L, t2, dim_new);
		      MatrixMult(t2, Sadv, t1, dim_new);
		      
		      
		      double cond = creal(MatrixTrace(t1, dim_new));
		      FreeMatrix(t1); FreeMatrix(t2); FreeMatrix(Sadv); FreeMatrix(S_temp);
		      FreeMatrix(Gamma_R);
		
		
		
	    //Now build left-to-right, updating diagonals and off-diagonals
		
		
		double _Complex *diagGs = createCompArray(Nrem);
		double _Complex **offdiags = createNonSquareMatrix(Nrem, chaininfo[0][0]);
		
		//copy of first chain info
		for(i=0; i<chaininfo[0][0]; i++)
		{
		  diagGs[i] = allgs[i][i];
		  
		  for(j=0; j<chaininfo[0][0]; j++)
		  {
		    offdiags[i][j] = allgs[i][j];
		  }
		}
		
		
		double _Complex **Gnew, **gold, **off_old, **off_new, **Gprev;
		Gprev=createSquareMatrix(chaininfo[0][0]);
		MatrixCopyPart(allgs, Gprev, 0, 0, 0, 0, chaininfo[0][0], chaininfo[0][0]);
		int dim0 = chaininfo[0][0];
		off_old = createNonSquareMatrix(dim0, dim0);
		MatrixCopy(offdiags, off_old, dim0);

		
		for(chain = 1; chain < length2; chain ++)
		{
		    //printf("# folding forward %d of %d\n", chain, length2);
		    dim_old = chaininfo[chain-1][0];
		    dim_new = chaininfo[chain][0];
		    
		    //gold is the current chain connected from the right 
		      gold = createSquareMatrix(dim_new);
		      MatrixCopyPart(allgs, gold, chaininfo[chain][1], 0, 0, 0, dim_new, dim_new);
		      
		    //Gprev is the completely connected previous cell
		      //is filled at the end of the iteration by copying Gnew
		      //has dimension dim_old
		      
		    //Gnew is the current chain completely connected
		      Gnew = createSquareMatrix(dim_new);
		      
		    //off_old 
		      //is filled at the end of the iteration by copying Gnew
		      //has dimension (dim_old, dim0) and is the connected version
		      
		    //off_new
		      off_new = createNonSquareMatrix(dim_new, dim0);
		    
		    
		    //generate Vs similar to before!
		      VLR = createNonSquareMatrix(dim_old, dim_new);
		      VRL = createNonSquareMatrix(dim_new, dim_old);
			  
			 k=0;
			      for(i=0; i<2*length1; i++)
			      {
				  if(siteinfo[(chain-1)*2*length1 + i][0] == 0)
				  {
				      l=0;
				      for(j=0; j< 2*length1; j++)
				      {
					  if(siteinfo[(chain)*2*length1 + j][0] == 0)
					  {
					      
					    VLR[k][l] =  VL1temp[i][j];
					    VRL[l][k] =  VL1temp[i][j];
					    l++;
					  }
				      }
				    k++;
				  }
			      }
			  
		    //update off-diagonals
		      t1 = createNonSquareMatrix(dim_new, dim_old);
		      MatrixMultNS(gold, VRL, t1, dim_new, dim_new, dim_old);
		      MatrixMultNS(t1, off_old, off_new, dim_new, dim_old, dim0);
		      
		    //update diagonals
		      t2 = createNonSquareMatrix(dim_new, dim_old);
		      MatrixMultNS(t1, Gprev, t2, dim_new, dim_old, dim_old);
		      FreeMatrix(t1);
		      t1 = createNonSquareMatrix(dim_new, dim_new);
		      MatrixMultNS(t2, VLR, t1, dim_new, dim_old, dim_new);
		      FreeMatrix(t2); 
		      t2= createSquareMatrix(dim_new);
		      MatrixMult(t1, gold, t2, dim_new);
		      FreeMatrix(t1);
		      MatrixAdd(gold, t2, Gnew, dim_new);
		      FreeMatrix(t2);
		  
		  
		    //copying and freeing
		      FreeMatrix(gold);
		      FreeMatrix(VLR);
		      FreeMatrix(VRL);
		      
		      FreeMatrix(Gprev);
		      Gprev = createSquareMatrix(dim_new);
		      MatrixCopy(Gnew, Gprev, dim_new);
		      for(i=0; i<dim_new; i++)
		      {
			  diagGs[chaininfo[chain][1] + i] = Gnew[i][i];
		      }
		      FreeMatrix(Gnew);
		      
		      FreeMatrix(off_old);
		      off_old = createNonSquareMatrix(dim_new, dim0);
		      MatrixCopyPart(off_new, off_old, 0,0,0,0,dim_new, dim0);
		      MatrixCopyPart(off_new, offdiags, 0, 0, chaininfo[chain][1], 0, dim_new, dim0);
		      FreeMatrix(off_new);
		  
		  
		  
		}
		
		
		//printf("%lf	%e	%e\n", creal(En), cimag(allgs[1][1]), cimag(diagGs[chaininfo[1][1] +1]));
		
		FreeMatrix(allgs); 
		FreeMatrix(off_old); FreeMatrix(Gprev);

		  
	
			    for(i=0; i<Nrem; i++)
			    {
			      ldoses[i] = (-1/M_PI) * cimag(diagGs[i]) ;
			    }
			    
			    
			    //current map outputs
				//need to index neighbours, and associate full and AL indices first
				//(this is a little sloppy, simply an index for each remaining atom)
				//would this be better in device generation step?
				j=0;
				for(i=0; i<2*length1*length2; i++)
				{
				  if((DeviceCell->siteinfo)[i][0] == 0)
				  {
				    indices[i] = j;
				    j++;
				  }
				  else
				  {
				    indices[i] = -1;
				  }
				}
				
				int atom, cell4, type, neighnum;
				
				for(i=0; i <2*length1*length2; i++)
				{
				  neigh[i][0] = -1; neigh[i][1] = -1; neigh[i][2] = -1;
				  neighnum = 0;

				    if((DeviceCell->siteinfo)[i][0] == 0)
				    {
					chain = i / (2*length1);
					atom = i - chain * 2 * length1;
					
					if(chain>0)
					{
					  for(j=0; j<2*length1; j++)
					  {
					      if(VL1temp[j][atom] != 0.0)
					      {
						if((DeviceCell->siteinfo)[2*length1*(chain-1) + j][0] == 0)
						{
						  neigh[i][neighnum] = 2*length1*(chain-1) + j;
						  neighnum++;
						}
					      }
					  }
					}
					
					for(j=0; j<2*length1; j++)
					{
					  if(j!=atom)
					  {
					    if(ginv[atom][j] != 0.0)
					    {
					      if((DeviceCell->siteinfo)[2*length1*(chain) + j][0] == 0)
					      {
						neigh[i][neighnum] = 2*length1*(chain) + j;
						neighnum++;
					      }
					    }
					  }
					}
					
					if(chain < length2-1)
					{
					  for(j=0; j<2*length1; j++)
					  {
					      if(VL1temp[atom][j] != 0.0)
					      {
						if((DeviceCell->siteinfo)[2*length1*(chain+1) + j][0] == 0)
						{
						  neigh[i][neighnum] = 2*length1*(chain+1) + j;
						  neighnum++;
						}
					      }
					  }
					}
					      
					
				   
				    
				    
				    
				    
					for(j=0;j<3; j++)
					{
					  if(indices[neigh[i][j]] != -1)
					  {
					      conds[indices[i]][j] = 0.0;
					      
					      for(k=0; k<dim0; k++)
					      {	
						  for(l=0; l<dim0; l++)
						  {
						    conds[indices[i]][j] += hopping * cimag (offdiags[indices[i]][k] * Gamma_L[k][l] * conj(offdiags[indices[neigh[i][j]]][l]) ) ;
						  }
					      }
					      
					      
					  }
					}
				    
				    }
				    
				}
				


		
		    FreeMatrix(ginv);
		    FreeMatrix(VL1temp);
		    FreeMatrix(offdiags);
		    FreeMatrix(Gamma_L);
		    free(diagGs);
		    
		    return cond;
		    
		    
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

  double t0;
  double kpar = para->kpar;
  double _Complex ans=0.0;
  double x1 = (aDeviceCell->pos)[a][0], y1 = (aDeviceCell->pos)[a][1];
  double x2 = (bDeviceCell->pos)[b][0] + bshifts[0], y2 = (bDeviceCell->pos)[b][1] + bshifts[1];
  double y2p, distp;
  
  double dist = sqrt(pow(x2-x1, 2.0) + pow(y2-y1, 2.0));
  double ycelldist;
  int i;
  //if((para->isperiodic)==0)
  //{
  
    //which coupling to use
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(dist > (para->NN_lowdis[i]) && dist < (para->NN_highdis[i]))
      {
           t0 = hops[i];
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
           t0 = hops[i];
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
           t0 = hops[i];
      }
    }
    ans+=t0*cexp(I*kpar); 
    
  }
  //printf("# hopping %d	%d: %lf	%lf\n", a, b, ans, dist);
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

  double t0;
  double kpar = para->kpar;
  double _Complex ans=0.0;
  double x1 = (aDeviceCell->pos)[a][0], y1 = (aDeviceCell->pos)[a][1];
  double x2 = (bDeviceCell->pos)[b][0] + bshifts[0], y2 = (bDeviceCell->pos)[b][1] + bshifts[1];
  double y2p, distp;
  
  double dist = sqrt(pow(x2-x1, 2.0) + pow(y2-y1, 2.0));
  double ycelldist;
  int i;
  
  //if((para->isperiodic)==0)
  //{
//     if(dist > (para->NN_lowdis) && dist < (para->NN_highdis))
//     {
//      ans=t0*graphenePeierlsPhase(x1, y1, x2, y2, para->gauge, para->Btes, para->restrics, para->limits);     
//     }
//   //}
  
    //which coupling to use
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(dist > (para->NN_lowdis[i]) && dist < (para->NN_highdis[i]))
      {
           t0 = hops[i];
      }
      
      if(dist == 0.0)
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
//      ans+=t0*graphenePeierlsPhase(x1, y1, x2, y2p, para->gauge, para->Btes, para->restrics, para->limits)*cexp(-I*kpar);    
//     }
    
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(distp > (para->NN_lowdis[i]) && distp < (para->NN_highdis[i]))
      {
           t0 = hops[i];
      }
    }
    ans+=t0*graphenePeierlsPhase(x1, y1, x2, y2p, para->gauge, para->Btes, para->restrics, para->limits)*cexp(-I*kpar);    
    
    
    
    //check if b is in cell below
    y2p = y2 - ycelldist;
    distp= sqrt(pow(x2-x1, 2.0) + pow(y2p-y1, 2.0));
    
//     if(distp > (para->NN_lowdis) && distp < (para->NN_highdis))
//     {
//      ans+=t0*graphenePeierlsPhase(x1, y1, x2, y2p, para->gauge, para->Btes, para->restrics, para->limits)*cexp(I*kpar);    
//     }
//     
    t0=0.0;
    for(i=0; i< num_neigh; i++)
    {
      if(distp > (para->NN_lowdis[i]) && distp < (para->NN_highdis[i]))
      {
           t0 = hops[i];
      }
    }
    ans+=t0*graphenePeierlsPhase(x1, y1, x2, y2p, para->gauge, para->Btes, para->restrics, para->limits)*cexp(I*kpar);    
    
    
    
  }
  //printf("# hopping %d	%d: %lf	%lf\n", a, b, ans, dist);
    //printf("%lf	%lf	%e\n", x1 + (x2-x1)/2, y1 + (y2-y1)/2, fabs(creal(ans)));

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
	  
	  
// 	    if(leadloop==0)
// 	  {
// // 	    printf("DIM %d\n", dim1);
// 	    listNonZero(ginv, dim1, dim1);
// 	    listNonZero(V12, dim1, dim1);
// 	    listNonZero(V21, dim1, dim1);
// 	  }
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
  
  
//gate induced potential - variations on the efetov model.
//places a potential on sites depending on gate and geometry parameters

void gate_induced_pot ( int vgtype, RectRedux *DeviceCell, double *engdeppots, double gate_voltage, double edge_cut_off, double subs_thick, double subs_epsr)
{
  int length1 = (DeviceCell->length);
  int length2 = (DeviceCell->length2);
  int geo = (DeviceCell->geo);
  int Nrem = *(DeviceCell->Nrem);
  double **pos = (DeviceCell->pos);
  
  double top_edge, bottom_edge, rib_center, ypos, ny, nmax, n0;
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
	n0= (subs_epsr)*(8.85E-12)*gate_voltage/(2*subs_thick * 1.6E-19);
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
	  
	  engdeppots[i] = (- gate_voltage/ fabs(gate_voltage)) * (1.05E-28) * (ny/fabs(ny))* sqrt(M_PI * fabs(ny) ) / (2.7 * 1.6E-19) ;
	  
	  if(gate_voltage==0.0)
	  {
	    engdeppots[i] = 0.0;
	  }
	}
  
}
