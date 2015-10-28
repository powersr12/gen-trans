#include "transport.h"
#include "test.h"

void genTransmissions(double _Complex En, RectRedux *DeviceCell, RectRedux **Leads, cnxProfile *cnxp, 
		      cellDivision *cellinfo, hoppingfunc *hoppingfn, void *hoppingparams,
		      lead_para *leadsparams, trans_params *tpara)
{
    //important definitions and quantities
    int length1 = (DeviceCell->length);
    int length2 = (DeviceCell->length2);
    int geo = (DeviceCell->geo);
    
    int cell1dim = (cellinfo->cell1dim);
    int num_leads = tpara->num_leads;
    
    int this_cell, dim, dim_old=0, dim_new;
    double _Complex **g_old, g_new, g00inv;
    double _Complex **g_sys_r, **g_sys_a, **SigmaR, **SigmaA, **Gamma, **temp1, **temp2;
    double _Complex **Gamma1, **Gamma2, **gsr, **gsa;
    int d1, d2, s1, s2;
    
    int i, j, k;
    
    SigmaR = createSquareMatrix(cell1dim);
    SigmaA = createSquareMatrix(cell1dim);
    g_sys_r = createSquareMatrix(cell1dim);
    g_sys_a = createSquareMatrix(cell1dim);
    Gamma = createSquareMatrix(cell1dim);

      leadfunction *leadfn = (leadfunction *)(leadsparams->leadsfn);

      //lead sigmas 
	 (leadfn)(En, DeviceCell, Leads, cellinfo, leadsparams, SigmaR);
	 


      //calculate retarded GF of system 
         genDeviceGF(En, DeviceCell, cnxp, cellinfo, hoppingfn, hoppingparams, 0, g_sys_r, NULL, SigmaR);
      
      
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
	 }
	 else if( (tpara->TRsym) == 1)
	 {
	    (leadfn)(creal(En) - I*cimag(En), DeviceCell, Leads, cellinfo, leadsparams, SigmaA);
	    genDeviceGF(creal(En) - I*cimag(En), DeviceCell, cnxp, cellinfo, hoppingfn, hoppingparams, 0, g_sys_a, NULL, SigmaA);
	 }
	 
	 
	 //Gamma
	  for(i=0; i<cell1dim; i++)
	  {
	    for(j=0; j<cell1dim; j++)
	    {
	      Gamma[i][j] = I*(SigmaR[i][j] - SigmaA[i][j]);
	    }
	  }
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
      
      
      
      FreeMatrix(Gamma);
      FreeMatrix(g_sys_r); FreeMatrix(g_sys_a);
        
  
  
}

//mode indicates either single sweep (mode=0) or double sweep (mode=1)
//(mode=2) for double sweep with reversed sweeping directions (for broken time-reversal symmetry cases)
void genDeviceGF(double _Complex En, RectRedux *DeviceCell, cnxProfile *cnxp, 
		      cellDivision *cellinfo, hoppingfunc *hoppingfn, void *hoppingparams, int mode, 
		      double _Complex **Gon, double _Complex **Goff, double _Complex **Sigma)
{
  
  int geo = (DeviceCell->geo);
    
  int num_cells = (cellinfo->num_cells);
  
  int this_cell, dim, dim_old=0, dim_new;
  double _Complex **g_old,  **g00inv, **V12, **V21, **smallSigma, **temp1, **temp2;
  int i, j, k, l, index1, index2, this0, last0;
  
  int cell_start, cell_end, cell_iter, it_count;
  int dim1  = (cellinfo->cell_dims)[0];
  double *bshifts = createDoubleArray(3);
  
  if(mode == 0 || mode == 1)
  {
    cell_start = num_cells-1;
    cell_end = 0;
    cell_iter = -1;
  }
    
  
  //calculate system GF	
  //backwards recursive sweep!
  for(it_count=0; it_count<(cellinfo->num_cells); it_count ++)
  {
	this_cell = cell_start + it_count*cell_iter;
  
//    	printf("########\n#cell %d\n", this_cell);
	
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
	  
	  g00inv[i][i] = En - (DeviceCell->site_pots)[index1];
	  
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
	if(mode==1)
	{
	}
	
	dim_old=dim;
      
  }//end of first sweep 
  
  MatrixCopy(g_old, Gon, dim1);
  FreeMatrix(g_old);
  
  if(mode>0)
  {
    
    //insert second sweep here
    
  }
  
   free(bshifts);
  
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
	      double elemerr=1.0e-10;
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
	      double elemerr=1.0e-8;
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
  double t0 = para->t0;
  double kpar = para->kpar;
  double _Complex ans=0.0;
  double x1 = (aDeviceCell->pos)[a][0], y1 = (aDeviceCell->pos)[a][1];
  double x2 = (bDeviceCell->pos)[b][0] + bshifts[0], y2 = (bDeviceCell->pos)[b][1] + bshifts[1];
  double y2p, distp;
  
  double dist = sqrt(pow(x2-x1, 2.0) + pow(y2-y1, 2.0));
  double ycelldist;
  
  //if((para->isperiodic)==0)
  //{
    if(dist > (para->NN_lowdis) && dist < (para->NN_highdis))
    {
     ans=t0;     
    }
  //}
  

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
    
    if(distp > (para->NN_lowdis) && distp < (para->NN_highdis))
    {
     ans+=t0*cexp(-I*kpar);    
    }
    
    //check if b is in cell below
    y2p = y2 - ycelldist;
    distp= sqrt(pow(x2-x1, 2.0) + pow(y2p-y1, 2.0));
    
    if(distp > (para->NN_lowdis) && distp < (para->NN_highdis))
    {
     ans+=t0*cexp(I*kpar); 
    }
    
    
  }
  //printf("# hopping %d	%d: %lf	%lf\n", a, b, ans, dist);
  return ans;
  
}



double _Complex peierlsTB(RectRedux *aDeviceCell, RectRedux *bDeviceCell, int a, int b, double *bshifts, void *hoppingparams)
{
  gen_hop_params *para = (gen_hop_params *)hoppingparams; 
  double t0 = para->t0;
  double kpar = para->kpar;
  double _Complex ans=0.0;
  double x1 = (aDeviceCell->pos)[a][0], y1 = (aDeviceCell->pos)[a][1];
  double x2 = (bDeviceCell->pos)[b][0] + bshifts[0], y2 = (bDeviceCell->pos)[b][1] + bshifts[1];
  double y2p, distp;
  
  double dist = sqrt(pow(x2-x1, 2.0) + pow(y2-y1, 2.0));
  double ycelldist;
  
  //if((para->isperiodic)==0)
  //{
    if(dist > (para->NN_lowdis) && dist < (para->NN_highdis))
    {
     ans=t0*graphenePeierlsPhase(x1, y1, x2, y2, para->gauge, para->Btes, para->restrics, para->limits);     
    }
  //}
  

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
    
    if(distp > (para->NN_lowdis) && distp < (para->NN_highdis))
    {
     ans+=t0*graphenePeierlsPhase(x1, y1, x2, y2p, para->gauge, para->Btes, para->restrics, para->limits)*cexp(-I*kpar);    
    }
    
    //check if b is in cell below
    y2p = y2 - ycelldist;
    distp= sqrt(pow(x2-x1, 2.0) + pow(y2p-y1, 2.0));
    
    if(distp > (para->NN_lowdis) && distp < (para->NN_highdis))
    {
     ans+=t0*graphenePeierlsPhase(x1, y1, x2, y2p, para->gauge, para->Btes, para->restrics, para->limits)*cexp(I*kpar);    
    }
    
    
  }
  //printf("# hopping %d	%d: %lf	%lf\n", a, b, ans, dist);
  return ans;
  
}

//returns peierls phase factor for two sites in graphene lattice
//note assumes the graphene lattice constant, so should be generalised for nongraphene
double _Complex graphenePeierlsPhase(double x1, double y1, double x2, double y2, int gauge, double BTesla, int *res, double **limits)
{
 
  double xb, yb, delx, dely, midx, midy;
  double beta = BTesla * 1.46262E-5;
  double phase;
  dely=y2-y1; delx=x2-x1;
  midx=x1 + delx/2; midy = y1+dely/2;
  

  
  
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
        phase = beta*yb*delx;

    
  }
//  printf("# %lf (%lf,  %lf) phase %+e\n", x, delx, dely, phase);
  
  return cexp(2 * M_PI * I * phase);  
  
  
}


//simplest way to generate sigmas
void simple2leads (double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, lead_para *params, double _Complex **Sigma)
{
  

	int geo = (DeviceCell->geo);
	int i, j, k;
	int dim1, dim2; 
  	double elemerr=1.0e-10;
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

	  InvertMatrixGSL(ginv, g00, dim1);
	  RubioSGF(SR, g00, V12, V21, dim1, &lcount, elemerr*dim2*dim2);
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
      ginv[i][i] = En - (LeadCell->site_pots[i]);
      
      for(j=0; j<dim; j++)
      {
	if(j!=i)
	{
	  ginv[i][j] = - (hopfn)(LeadCell, LeadCell, i, j, bshifts, (params->hoppara) );
	}
	
      }
            
    }

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
  
  

