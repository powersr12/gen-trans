#include "transport.h"

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