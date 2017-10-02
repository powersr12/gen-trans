#include "transport.h"
#include "test.h"


//mode=0 - just transmissions!
//mode=1 - transmissions and ldos/current maps - can be generalised further later
void genTransmissions(double _Complex En, RectRedux *DeviceCell, RectRedux **Leads, cnxProfile *cnxp, cellDivision *cellinfo, hoppingfunc *hoppingfn, gen_hop_params *hoppingparams,
		      lead_para *leadsparams, trans_params *tpara, int mode, double *ldoses, double ***currents)
{


	//device parameters, modes, matrix declarations etc.
		int geo = (DeviceCell->geo);
		
		int tcell1dim = (cellinfo->cell1dim);
		int num_leads = tpara->num_leads;
		
		int this_cell, dim, dim_old=0, dim_new;
		int tNrem = *(DeviceCell->Nrem);
		int tNtot = *(DeviceCell->Ntot);
		int spindep = (DeviceCell->spindep);
		int cell1dim, Nrem, Ntot;
		
		if(spindep == 0 || spindep == 2)
		{
			cell1dim = tcell1dim;
			Nrem = tNrem;
			Ntot = tNtot;
		}
		if(spindep == 1)
		{
			cell1dim = 2*tcell1dim;
			Nrem = 2*tNrem;
			Ntot = 2*tNtot;
		}
			
		

		double _Complex **g_old, g_new, g00inv;
		double _Complex **g_sys_r, **g_sys_a, **SigmaR, **SigmaA, **Gamma, **temp1, **temp2, *gii, **gi1, **g1i;
		double _Complex **tSigmaR, **tg_sys_r, *tgii, **tgi1, **tg1i;
		double _Complex **Gamma1, **Gamma2, **gsr, **gsa;
		int d1, d2, s1, s2;
		
		int i, j, k, l, m, n;
		
		int devicemode=0, devicemode2=0; 
		int spinA = 0, spinB=0;
		
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
    
    
	//GFs for LDOS/currents
		gii=NULL; gi1=NULL; g1i=NULL;
		tgii=NULL; tgi1=NULL; tg1i=NULL;
		
          
    
	//retarded and advanced self energies and device GFs
		//matrix sizes are now dependent on the spin mode
		if(spindep == 0 || spindep == 1)
		{
			SigmaR = createSquareMatrix(cell1dim);
			SigmaA = createSquareMatrix(cell1dim);
			g_sys_r = createSquareMatrix(cell1dim);
			g_sys_a = createSquareMatrix(cell1dim);
			Gamma = createSquareMatrix(cell1dim);
			
			if(mode ==1)
			{
				gii = createCompArray(Ntot);
				gi1 = createNonSquareMatrix(Ntot, cell1dim);
				g1i = createNonSquareMatrix(cell1dim, Ntot);
			}
			
		}
		if(spindep == 2)
		{
			SigmaR = createNonSquareMatrix(2*cell1dim, cell1dim );
			SigmaA = createNonSquareMatrix(2*cell1dim, cell1dim);
			g_sys_r = createNonSquareMatrix(2*cell1dim, cell1dim);
			g_sys_a = createNonSquareMatrix(2*cell1dim, cell1dim);
			Gamma = createNonSquareMatrix(2*cell1dim, cell1dim);
			
			
			tSigmaR = createNonSquareMatrix(cell1dim, cell1dim );
			tg_sys_r = createNonSquareMatrix(cell1dim, cell1dim);
			
			if(mode ==1)
			{
				gii = createCompArray(2*Ntot);
				gi1 = createNonSquareMatrix(2*Ntot, cell1dim);
				g1i = createNonSquareMatrix(2*cell1dim, Ntot);
				
				tgii = createCompArray(Ntot);
				tgi1 = createNonSquareMatrix(Ntot, cell1dim);
			}
			
			
		}
		

		
			
	//specify lead generation function, and calculate self energies of the leads
		//lead functions take care of spin inside themselves! (multipleCustomLeads does..., need to fix others)
		leadfunction *leadfn = (leadfunction *)(leadsparams->leadsfn);
		(leadfn)(En, DeviceCell, Leads, cellinfo, leadsparams, SigmaR);
	 
		
	 
		//printEMatrix(SigmaR, cell1dim);
     
	//calculate retarded GF of system 
	if(spindep == 0)
		genDeviceGF(En, DeviceCell, cnxp, cellinfo, hoppingfn, hoppingparams, devicemode, devicemode2, 0, g_sys_r, gii, gi1, SigmaR);
	
	if(spindep == 1)
		genDeviceGF(En, DeviceCell, cnxp, cellinfo, hoppingfn, hoppingparams, devicemode, devicemode2, 2, g_sys_r, gii, gi1, SigmaR);
	
	//printEMatrix(SigmaR, cell1dim);
	
	if(spindep == 2)
	{
		//first spin
		MatrixCopyPart(SigmaR, tSigmaR, 0, 0, 0, 0, cell1dim, cell1dim);
		genDeviceGF(En, DeviceCell, cnxp, cellinfo, hoppingfn, hoppingparams, devicemode, devicemode2, 0, tg_sys_r, tgii, tgi1, tSigmaR);
		
		MatrixCopyPart(tg_sys_r, g_sys_r, 0, 0, 0, 0, cell1dim, cell1dim);
		
		if(mode ==1)
		{
			for(i=0;i<Ntot; i++)
				gii[i] = tgii[i];
			
			MatrixCopyPart(tgi1, gi1, 0, 0, 0, 0, Ntot, cell1dim);
		}
		
		//second spin
		MatrixCopyPart(SigmaR, tSigmaR, cell1dim, 0, 0, 0, cell1dim, cell1dim);
		genDeviceGF(En, DeviceCell, cnxp, cellinfo, hoppingfn, hoppingparams, devicemode, devicemode2, 1, tg_sys_r, tgii, tgi1, tSigmaR);
		
		MatrixCopyPart(tg_sys_r, g_sys_r, 0, 0, cell1dim, 0, cell1dim, cell1dim);
		
		if(mode ==1)
		{
			for(i=0;i<Ntot; i++)
				gii[Ntot+i] = tgii[i];
			
			MatrixCopyPart(tgi1, gi1, 0, 0, Ntot, 0, Ntot, cell1dim);
		}
		
		FreeMatrix(tg_sys_r);
		FreeMatrix(tgi1);
		FreeMatrix(tSigmaR);
		free(tgii);
		
	}
		


	//calculate advanced quantities from the retarded versions 
		for(i=0; i<cell1dim; i++)
		{
			for(j=0; j<cell1dim; j++)
			{
				SigmaA[i][j] = conj(SigmaR[j][i]);
				g_sys_a[i][j] = conj(g_sys_r[j][i]);	
				
				if(spindep==2)
				{
					SigmaA[cell1dim+i][j] = conj(SigmaR[cell1dim+j][i]);
					g_sys_a[cell1dim+i][j] = conj(g_sys_r[cell1dim+j][i]);	
				}
			}
		}
		
				
		if(mode==1)
		{
			for(i=0; i<cell1dim; i++)
			{
				for(j=0; j<Ntot; j++)
				{
					g1i[i][j] = conj(gi1[j][i]);
					
					if(spindep==2)
					{
						g1i[cell1dim+i][j] = conj(gi1[cell1dim+j][i]);
					}
				}
			}
		}
	   
	 
	 //calculate the Gamma (broadening) matrices for the leads
		for(i=0; i<cell1dim; i++)
		{
			for(j=0; j<cell1dim; j++)
			{
				Gamma[i][j] = I*(SigmaR[i][j] - SigmaA[i][j]);
				
				if(spindep==2)
				{
					Gamma[cell1dim+i][j] = I*(SigmaR[cell1dim+i][j] - SigmaA[cell1dim+i][j]);
				}
			}
		} 
		FreeMatrix(SigmaR); FreeMatrix(SigmaA); 
      
		
		
	//Calculate the transmissions between all the leads	
		//s1 and s2 keep track of cumulative dimensions (i.e starting indices)
		
		double tempt=0, tempt2=0, temptud=0, temptdu=0;
		//each lead pair has only a total electronic transmission
		if(spindep == 0)
		{
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
					
					tempt=0; tempt2=0; temptud=0; temptdu=0;
					for(k=0; k<d1; k++)
					{
						for(l=0; l<d2; l++)
						{
							tempt += creal( temp2[k][l] * gsa[l][k]);
						}
					}
					(tpara->transmissions)[i][j] = tempt;
					FreeMatrix(temp2);
					
					
					
// 					temp1=createSquareMatrix(d1);
// 					MatrixMultNS(temp2, gsa, temp1, d1, d2, d1);
// 					FreeMatrix(temp2);
// 					
// 					(tpara->transmissions)[i][j] = creal(MatrixTrace(temp1, d1));
// 					FreeMatrix(temp1);
					
					
					s2+=d2;
					FreeMatrix(gsr); FreeMatrix(Gamma2); FreeMatrix(gsa);
				}
				s1+=d1;
				FreeMatrix(Gamma1);
			}
		}
		
		
	
	
		//each lead pair has a total electronic transmission, and FOUR spin channel transmissions
		if(spindep == 1)
		{
			s1=0; 
			for(i=0; i<num_leads; i++)
			{
				d1=(cellinfo->lead_dims)[i];
				s2=0;
				
					
				//Gamma is now spin dependent!
				Gamma1 = createSquareMatrix(2*d1);
				MatrixCopyPart(Gamma, Gamma1, s1, s1, 0, 0, d1, d1);
				MatrixCopyPart(Gamma, Gamma1, tcell1dim+s1, tcell1dim+s1, d1, d1, d1, d1);	
				MatrixCopyPart(Gamma, Gamma1, s1, tcell1dim+s1, 0, d1, d1, d1);	
				MatrixCopyPart(Gamma, Gamma1, tcell1dim+s1, s1, d1, 0, d1, d1);
				
				//printEMatrix(Gamma1, 2*d1);
				
					
					for(j=0; j< num_leads; j++)
					{
						d2=(cellinfo->lead_dims)[j];
						gsr = createNonSquareMatrix(2*d1, 2*d2);
						MatrixCopyPart(g_sys_r, gsr, s1, s2, 0, 0, d1, d2);
						MatrixCopyPart(g_sys_r, gsr, tcell1dim+s1, tcell1dim+s2, d1, d2, d1, d2);
						MatrixCopyPart(g_sys_r, gsr, s1, tcell1dim+s2, 0, d2, d1, d2);
						MatrixCopyPart(g_sys_r, gsr, tcell1dim+s1, s2, d1, 0, d1, d2);
						
						//printEMatrix(gsr, 2*d1);
						
						Gamma2 = createSquareMatrix(2*d2);
						MatrixCopyPart(Gamma, Gamma2, s2, s2, 0, 0, d2, d2);
						MatrixCopyPart(Gamma, Gamma2, tcell1dim+s2, tcell1dim+s2, d2, d2, d2, d2);	
						MatrixCopyPart(Gamma, Gamma2, s2, tcell1dim+s2, 0, d2, d2, d2);	
						MatrixCopyPart(Gamma, Gamma2, tcell1dim+s2, s2, d2, 0, d2, d2);
						
						
						gsa = createNonSquareMatrix(2*d2, 2*d1);
						MatrixCopyPart(g_sys_a, gsa, s2, s1, 0, 0, d2, d1);
						MatrixCopyPart(g_sys_a, gsa, tcell1dim+s2 , tcell1dim+s1, d1, d2, d2, d1);
						MatrixCopyPart(g_sys_a, gsa, s2 , tcell1dim+s1, 0, d2, d2, d1);
						MatrixCopyPart(g_sys_a, gsa, tcell1dim+s2 , s1, d1, 0, d2, d1);
						
						
						temp1=createNonSquareMatrix(2*d1, 2*d2);
						MatrixMultNS(Gamma1, gsr, temp1, 2*d1, 2*d1, 2*d2);
						temp2=createNonSquareMatrix(2*d1, 2*d2);
						MatrixMultNS(temp1, Gamma2, temp2, 2*d1, 2*d2, 2*d2);
						FreeMatrix(temp1);
						
						
						tempt=0; tempt2 =0; temptud=0; temptdu=0;
						for(k=0; k<d1; k++)
						{
							for(l=0; l<2*d2; l++)
							{
								tempt += creal( temp2[k][l] * gsa[l][k]);
								tempt2 += creal( temp2[d1+k][l] * gsa[l][d1+k]);
								temptud += creal( temp2[k][l] * gsa[l][d1+k]);
								temptdu += creal( temp2[d1+k][l] *gsa[l][k]);
							}
						}
						(tpara->spin_trans)[i][j] = tempt;
						(tpara->spin_trans)[num_leads+i][num_leads+j] = tempt2;
						(tpara->spin_trans)[i][num_leads+j] = temptud;
						(tpara->spin_trans)[num_leads+i][j] = temptdu;
						(tpara->transmissions)[i][j] = tempt +tempt2 + temptud + temptdu;
						//printf("# %lf	%lf	%lf	%lf\n", tempt ,tempt2 , temptud , temptdu);
						FreeMatrix(temp2);
						
						
						s2+=d2;
						FreeMatrix(gsr); FreeMatrix(Gamma2); FreeMatrix(gsa);
						
					}
					s1+=d1;
					FreeMatrix(Gamma1);
					
				
			}
		}
		
		
		
		
		
		
	//Quantities needed for LDOS/current mapping
		double *bshifts = createDoubleArray(3);
		double current_mag;
		double _Complex hopmag;
		FILE *bigdump;
		char bigfile[200];
	
	

		
		//THIS IS NOT YET SPIN GENERALISED!
	//Calculate and output LDOS and current maps if required
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
			
			//loop over sites
			for(i=0; i<Ntot; i++)
			{
				ldoses[i] = -cimag(gii[i])/M_PI;
				s1=0; d1=0;
				
				//get contribution from each lead
				for(k=0; k<num_leads; k++)
				{
					d1=(cellinfo->lead_dims)[k];
					currents[k][i][0] = 0.0;
					currents[k][i][1] = 0.0;
			
					//sum over neighbours of the site
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
						// 		  fprintf(bigdump, "%lf	%lf	%lf	%lf	%.15e\n", (DeviceCell->pos)[i][0], 			(DeviceCell->pos)[i][1], (DeviceCell->pos)[j][0],(DeviceCell->pos)[j][1], current_mag);
						// 		fclose(bigdump);
				
			
						//only consider positive (outgoing) contributions
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
      
     
       
	//Free up memory
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

//spin channel: 0 (up), 1(down) or 2(both together)
void genDeviceGF(double _Complex En, RectRedux *DeviceCell, cnxProfile *cnxp, cellDivision *cellinfo, hoppingfunc *hoppingfn, gen_hop_params *hoppingparams, int mode, int mode2, int spinchannel, double _Complex **Gon, double _Complex *Gdiags, double _Complex **Goff, double _Complex **Sigma)
{

	
	//device parameters, modes, variables
		int geo = (DeviceCell->geo);
		int num_cells = (cellinfo->num_cells);
		int this_cell, dim, dim_old=0, dim_new;
		double _Complex **g_old,  **g00inv, **V12, **V21, **smallSigma, **temp1, **temp2;
		int i, j, k, l, index1, index2, this0, last0;
		int cell_start, cell_end, cell_iter, it_count;
		int dim1  = (cellinfo->cell_dims)[0];
		int Nrem = *(DeviceCell->Nrem);
		int Ntot = *(DeviceCell->Ntot);
		int spindep = (DeviceCell->spindep);
		int are_spin_pots = (DeviceCell->are_spin_pots);	//possibly not needed
		double **spin_pots = (DeviceCell->spin_pots);
		double *const_spin_pots= (DeviceCell->const_spin_pots);
		int spinA, spinB, spin1, spin2;
		int fac=1; //multiplication factor for matrix sizes due to spin

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
		
		
	//NORMAL, SPIN-UNPOLARISED RUN SETTINGS
		if(spindep==0)
		{
			
			
		}
		
	//TWO INDEPENDENT SPIN CHANNELS SETTINGS
		//spin channel chosen is given by spinchannel (0=spin_up, 1=spindown)
		if(spindep==2)
		{
		
			
		}
  
	//COMPLETELY SPIN DEPENDENT SETTINGS
		//spin channel chosen is given by spinchannel (0=spin_up, 1=spindown)
		if(spindep==1 || spinchannel == 2)
		{
			//ensure both params are correct for later loops!
			spindep =1;
			spinchannel = 2;
			
			fac=2;
		
			
		}
  
  
	
	//a large storage space for a lot of Green's functions used in double sweep recursive GF method
		if(mode>0)
		{
			//requests memory for an array of pointers to matrices
			//the individual matrices here can be spin-dependent, but does not effect this declaration
			allgs = (double _Complex ***)malloc(num_cells * sizeof(double _Complex **));
		}

		
	
	
	//SWEEP 1 - backwards recursive sweep from the final cell
		for(it_count=0; it_count<(cellinfo->num_cells); it_count ++)
		{
			this_cell = cell_start + it_count*cell_iter;
		
				//printf("########\n#cell %d\n", this_cell);
				//printf("### backwards sweep, cell %3d of %3d\n", this_cell, (cellinfo->num_cells));
			
			
			//generate cell GF
				dim = (cellinfo->cell_dims)[this_cell];
				g00inv = createSquareMatrix(fac*dim);
			
				this0=(cellinfo->starting_index)[this_cell];
				if(it_count>0)
				{
					last0=(cellinfo->starting_index)[this_cell-cell_iter];
					V21 = createNonSquareMatrix(fac*dim, fac*dim_old);
					V12 = createNonSquareMatrix(fac*dim_old, fac*dim);
				}
			
			
			//onsites 
				//contributions from site_pots, and onsite hopping correction terms
				for(i=0; i<dim; i++)
				{
					index1 = (cellinfo->cells_site_order)[this0 +i];
					
					//onsite terms for the single spin case
					if(spinchannel != 2)
					{
						(hoppingparams->spinA)=spinchannel;
						(hoppingparams->spinB)=spinchannel;
						
						//hopping onsites
						g00inv[i][i] = En - hoppingfn(DeviceCell, DeviceCell, index1, index1, bshifts, hoppingparams);
						
						//regular onsites
						g00inv[i][i] -= (DeviceCell->site_pots)[index1];
						
						//spin potential onsites
						if(spindep==1 && are_spin_pots==1)
						{
							if(spinchannel==0)
							{
								g00inv[i][i] -= (const_spin_pots[2] + spin_pots[index1][2]);
							}
							if(spinchannel==1)
							{
								g00inv[i][i] -= (-const_spin_pots[2] - spin_pots[index1][2]);
							}
						}
					}
					
					//onsite terms for both spins
					if(spinchannel == 2)
					{
						for(spin1 =0; spin1< 2; spin1++)
						{
							(hoppingparams->spinA)=spin1;
							
							//energy term
							g00inv[spin1*dim + i][spin1*dim + i] = En;	
							
							//spin-independent onsites
							g00inv[spin1*dim + i][spin1*dim + i] -= (DeviceCell->site_pots)[index1];
							
							for(spin2=0; spin2<2; spin2++)
							{
								(hoppingparams->spinB)=spin2;
								
								//hopping onsite terms, possibly spin dep
								g00inv[spin1*dim + i][spin2*dim + i] -= hoppingfn(DeviceCell, DeviceCell, index1, index1, bshifts, hoppingparams);
								
								
								//spin-dependent potentials
								
								if(are_spin_pots==1)
								{
									if(spin1==0 && spin2==0)
									{
										g00inv[spin1*dim + i][spin2*dim + i] -= (const_spin_pots[2] + spin_pots[index1][2]);
									}
									if(spin1==0 && spin2==1)
									{
										g00inv[spin1*dim + i][spin2*dim + i] -= (const_spin_pots[0] + spin_pots[index1][0] - I*(const_spin_pots[1] + spin_pots[index1][1]));
									}
									
									if(spin1==1 && spin2==0)
									{
										g00inv[spin1*dim + i][spin2*dim + i] -= (const_spin_pots[0] + spin_pots[index1][0] + I*(const_spin_pots[1] + spin_pots[index1][1]));
									}
									
									if(spin1==1 && spin2==1)
									{
										g00inv[spin1*dim + i][spin2*dim + i] -= (-const_spin_pots[2] - spin_pots[index1][2]);
									}
								}
							}
						}
								
					}
					
				}
// 				printEMatrix(g00inv, dim);

			
			
			//internal & external cell hoppings
				//loop over sites in cell
				for(i=0; i<dim; i++)
				{
					index1 = (cellinfo->cells_site_order)[this0 +i];
					
					//loop over neighbours of site
					for(j=0; j<(cnxp->site_cnxnum)[index1]; j++)
					{
						index2 = (cnxp->site_cnx)[index1][j];
					
						//if this neighbour in the same cell, calculate hopping param and place in unit cell Hamiltonian
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
							
							if(spinchannel != 2)
							{
								(hoppingparams->spinA)=spinchannel;
								(hoppingparams->spinB)=spinchannel;
								g00inv[i][k] -= hoppingfn(DeviceCell, DeviceCell, index1, index2, bshifts, hoppingparams);
							}
							
							if(spinchannel == 2)
							{
								for(spin1=0; spin1<2; spin1++)
								{
									(hoppingparams->spinA)=spin1;
									for(spin2=0; spin2<2; spin2++)
									{
										(hoppingparams->spinB)=spin2;
										g00inv[spin1*dim+i][spin2*dim+k] -= hoppingfn(DeviceCell, DeviceCell, index1, index2, bshifts, hoppingparams);
									}									
								}
							}
							
							
						}
					
					
						//if this neighbour in the previous cell, add hopping to the connection matrices 
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
							
							if(spinchannel != 2)
							{
								(hoppingparams->spinA)=spinA;
								(hoppingparams->spinB)=spinB;
								V21[i][k] = hoppingfn(DeviceCell, DeviceCell, index1, index2, bshifts, hoppingparams);
								V12[k][i] = hoppingfn(DeviceCell, DeviceCell, index2, index1, bshifts, hoppingparams);
							}
							
							
							if(spinchannel == 2)
							{
								for(spin1=0; spin1<2; spin1++)
								{
									(hoppingparams->spinA)=spin1;
									for(spin2=0; spin2<2; spin2++)
									{
										(hoppingparams->spinB)=spin2;
										
										V21[spin1*dim+i][spin2*dim_old+k] = hoppingfn(DeviceCell, DeviceCell, index1, index2, bshifts, hoppingparams);
										V12[spin1*dim_old+k][spin2*dim+i] = hoppingfn(DeviceCell, DeviceCell, index2, index1, bshifts, hoppingparams);
									}
								}
							}
							
						}
					}
				}

						//check non-zero elements at the cell stage
//  						printf("%d cell \n#g00inv\n", it_count);
//  						listNonZero(g00inv, dim, dim);
//  						if(it_count>0)
//  						{
//  							  printf("\n#V21\n");
//  							  listNonZero(V21, dim, dim_old);
// 							  printf("\n#V12\n");
//  							  listNonZero(V12, dim_old, dim);
//  						}
		
		
			

			
			//generate self energy term from previous cells
				smallSigma=createSquareMatrix(fac*dim);
				
				if(it_count>0)
				{
					//previous cells self energy: smallSigma = V21 g_old V12
					temp1 = createNonSquareMatrix(fac*dim, fac*dim_old);
					MatrixMultNS(V21, g_old, temp1, fac*dim, fac*dim_old, fac*dim_old);
					MatrixMultNS(temp1, V12, smallSigma, fac*dim, fac*dim_old, fac*dim);
					FreeMatrix(temp1);
				}
		
			
			
			
			//account for lead self energies if this is cell 0
				if(this_cell==0)
				{
					temp1 = createSquareMatrix(fac*dim);
					MatrixAdd(smallSigma, Sigma, temp1, fac*dim);
					MatrixCopy(temp1, smallSigma, fac*dim);
					FreeMatrix(temp1);
				}
				
				
				
			//add self energies to the GF and update the edge GF in g_old
			//free memory used in this cell iteration
				temp1 = createSquareMatrix(fac*dim);
				MatrixSubtract(g00inv, smallSigma, temp1, fac*dim);
				if(it_count>0)
				{
					FreeMatrix(g_old);
				}
				g_old = createSquareMatrix(fac*dim);
				InvertMatrixGSL(temp1, g_old, fac*dim);
				FreeMatrix(temp1); FreeMatrix(g00inv); FreeMatrix(smallSigma);
				if(it_count>0)
				{
					FreeMatrix(V12); FreeMatrix(V21);
				}
			
			
			
			//save the temporary edges if a dual sweep will be performed below
				if(mode>0)
				{
					allgs[this_cell] = createSquareMatrix(fac*dim);
					MatrixCopy(g_old, allgs[this_cell], fac*dim);
				}
			
			
			//update dimension of edge cell
				dim_old=dim;
			
		
		}//end of first sweep 
		
	//

	MatrixCopy(g_old, Gon, fac*dim1);
	FreeMatrix(g_old);

  
	//SWEEP2 - forward sweep for LDOS/current map details 
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
			

					if(mode2==1)
					{
						for(i=0; i<dim1; i++)
						{
							Gdiags[(cellinfo->cells_site_order)[i]] = Gon[i][i];
							for(j=0; j<dim1; j++)
							{
								Goff[(cellinfo->cells_site_order)[i]][j] = Gon[i][j];
							}
						}
					}
					
					//mode 2 is now pretty much redundant, but kept here just in case
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

			//LOOP OVER CELLS STARTS HERE
				for(it_count=1; it_count<(cellinfo->num_cells); it_count ++)
				{
					//dimensions and starting indices
						this_cell = cell_end - it_count*cell_iter;
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
						//it is filled at the end of the loop, and before the 1st iteration begins

					//off_old 
						//old version of off diagonal GF
						//is filled at the end of the iteration by copying off_new
						//has dimension (dim_old, dim1) (or vice versa) and is the connected version
					
					//off_new - to be calculated below
						if(mode2 == 1)
						{
							off_new = createNonSquareMatrix(dim, dim1);
						}
						if(mode2 == 2)
						{
							off_new = createNonSquareMatrix(dim1, dim);
						}
					
					
					//generate V matrices as in reverse loop
						V12 = createNonSquareMatrix(dim_old, dim);
						V21 = createNonSquareMatrix(dim, dim_old);
					
						//loop over sites in cell
						for(i=0; i<dim; i++)
						{
							index1 = (cellinfo->cells_site_order)[this0 +i];
							
							//loop over neighbours
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
									
									(hoppingparams->spinA)=spinA;
									(hoppingparams->spinB)=spinB;
									V21[i][k] = hoppingfn(DeviceCell, DeviceCell, index1, index2, bshifts, hoppingparams);
										
									//(hoppingparams->spinA)=spinB;
									//(hoppingparams->spinB)=spinA;
									V12[k][i] = hoppingfn(DeviceCell, DeviceCell, index2, index1, bshifts, hoppingparams);
								}
							}
						}
					
					
					//updates if mode 1
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
					
						
							//Freeing, and copying what needs to be saved
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
					
					//updates if mode 2
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
		}//SWEEP 2 COMPLETE

		
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
  int spinA = para->spinA;
  int spinB = para->spinB;
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

  
  //this hopping routine is zero for elements off-diagonal in spin-space
  if(spinA!=spinB)
  {
	  ans=0.0;
  }
  else
  {
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
	
  }
  
  return ans;
  
}


//returns hopping parameters for graphene with UNIFORM induced SOC terms.
//(initially just intrinsic SOC, rashba to come...) 
//(include possibility of Peierl's magnetic field phases also... (later!))
//hoppings are in the sequence:
//0=NNTB, 1=lambda_IA, 2=lambda_IB 

double _Complex grapheneSOCTB(RectRedux *aDeviceCell, RectRedux *bDeviceCell, int a, int b, double *bshifts, void *hoppingparams)
{
  gen_hop_params *para = (gen_hop_params *)hoppingparams; 
  int num_neigh = para->num_neigh;
  double _Complex *hops = para->hops;
  double  *NN_lowdis = para->NN_lowdis;
  double  *NN_highdis = para->NN_highdis;
  double  *NN_shifts = para->NN_shifts;
  double  *NN_zmin= para->NN_zmin;
  double  *NN_zmax = para->NN_zmax;
  int spinA = para->spinA;
  int spinB = para->spinB;
  double _Complex t0;
  double kpar = para->kpar;
  double _Complex ans=0.0;
  double x1 = (aDeviceCell->pos)[a][0], y1 = (aDeviceCell->pos)[a][1];
  double x2 = (bDeviceCell->pos)[b][0] + bshifts[0], y2 = (bDeviceCell->pos)[b][1] + bshifts[1];
  double y2p, distp;
  
  double dist = sqrt(pow(x2-x1, 2.0) + pow(y2-y1, 2.0));
  double zdiff = fabs((bDeviceCell->pos)[b][2] - (aDeviceCell->pos)[a][2]);
  double ycelldist;
  int i, j, k=0, l;
  int SOCsign, possign;
  int *possible = createIntArray(20);
  double midx=0.0, midy=0.0;
  double dista=0.0, distb=0.0;
  double v11, v12, v21, v22, crsp=0;
  
  cnxProfile *acnxp = (cnxProfile *)(aDeviceCell->cnxp);
  cnxProfile *bcnxp = (cnxProfile *)(bDeviceCell->cnxp);

  
  //this hopping routine is zero for elements off-diagonal in spin-space
  if(spinA!=spinB)
  {
	  ans=0.0;
  }
  else
  {
	//sign of SOC term coming from spin orientation
	SOCsign = 1;
	if(spinA == 1)
		SOCsign = -1;
	
	//sign of SOC term coming from clockwise/anticlockwise ((anti)clockwise from b to a) 
	
		
		//initially, multiplication factor zero
		//will be changed to 1 if atoms are 2NNs.
		//then sign changed according to clockwise or anticlockwise
			possign=0;
	
			
		//only run if we need SOC terms and if intrinsic params are non-zero
		if(num_neigh>0 && (hops[1] !=0.0 || hops[2] != 0.0) && (dist >= (para->NN_lowdis[1]) && dist < (para->NN_highdis[1])))
		{
			//case for both atoms in the same region (lead/device) AND bshifts=0.0
			//(indices will be the same here...)
			if ( ((aDeviceCell->islead) == (bDeviceCell->islead)) && bshifts[0]==0.0 &&  bshifts[1]==0.0)
			{
				k=0;
				//finds and counts common neighbours of two atoms in same device
				for(i=0; i< (acnxp->site_cnxnum)[a]; i++)
				{
					for(j=0; j<(bcnxp->site_cnxnum)[b]; j++)
					{
						if((acnxp->site_cnx)[a][i] == (bcnxp->site_cnx)[b][j])
						{
							possible[k] = (acnxp->site_cnx)[a][i];
							k++;
						}
					}
				}
				
				//
				i=0; j=0;
				while(j!=1 && i<k)
				{
					midx = (aDeviceCell->pos)[possible[i]][0];
					midy=  (aDeviceCell->pos)[possible[i]][1];
					dista = sqrt(pow(midx-x1, 2.0) + pow(midy-y1, 2.0));
					distb = sqrt(pow(midx-x2, 2.0) + pow(midy-y2, 2.0));
				
					if(dista >= (para->NN_lowdis[0]) && dista < (para->NN_highdis[0]) && distb >= (para->NN_lowdis[0]) && distb < (para->NN_highdis[0]))
					{
						j=1;
						possign == 1;
					}
					i++;
				}	
				
				if(j==1)
				{
					//v1=vector from b to midpoint, v2 from midpoint to a
					v11= midx-x2; v12=midy-y2;
					v21= x1-midx; v22=y1-midy;
					
					crsp = v11*v22 - v12*v21;
					
					if(crsp > 0)
						possign = 1;
					if(crsp < 0)
						possign = -1;
				}
			
			}
			
			
			//case for atoms in different regions (lead AND device) (with bshifts=0)
			//(check positions rather than indices to determine 3rd atom)
			if ( ((aDeviceCell->islead) != (bDeviceCell->islead)) || bshifts[0]!=0.0 || bshifts[1]!=0.0)
			{
				k=0; j=0; i=0;
							
				
				//check if an NN of a is an NN of b
				while(j!=1 && i< (acnxp->site_cnxnum)[a])
				{
					//is this neighbour an NN of a AND b?
					midx=(aDeviceCell->pos)[(acnxp->site_cnx)[a][i]][0];
					midy=(aDeviceCell->pos)[(acnxp->site_cnx)[a][i]][1];
					dista = sqrt(pow(midx-x1, 2.0) + pow(midy-y1, 2.0));
					distb = sqrt(pow(midx-x2, 2.0) + pow(midy-y2, 2.0));
					
					if(dista >= (para->NN_lowdis[0]) && dista < (para->NN_highdis[0]) && distb >= (para->NN_lowdis[0]) && distb < (para->NN_highdis[0]))
					{
						j=1;
					}
					else
						i++;
				}
				
				i=0;
				//if not, check if an NN of b is an NN of a
				while(j!=1 && i< (bcnxp->site_cnxnum)[b])
				{
					//is this neighbour an NN of a AND b?
					midx=(bDeviceCell->pos)[(bcnxp->site_cnx)[b][i]][0] + bshifts[0];
					midy=(bDeviceCell->pos)[(bcnxp->site_cnx)[b][i]][1] + bshifts[1];
					dista = sqrt(pow(midx-x1, 2.0) + pow(midy-y1, 2.0));
					distb = sqrt(pow(midx-x2, 2.0) + pow(midy-y2, 2.0));
					
					if(dista >= (para->NN_lowdis[0]) && dista < (para->NN_highdis[0]) && distb >= (para->NN_lowdis[0]) && distb < (para->NN_highdis[0]))
					{
						j=1;
					}
					else
					{
						i++;
					}
				}
				
				if(j==1)
				{
				//v1=vector from b to midpoint, v2 from midpoint to a
					v11= midx-x2; v12=midy-y2;
					v21= x1-midx; v22=y1-midy;
					
					crsp = v11*v22 - v12*v21;
					
					if(crsp > 0)
						possign = 1;
					if(crsp < 0)
						possign = -1;
				}

			}
			
		}
	
	
	//which coupling to use
	t0=0.0;
	
	
	//Nearest neighbour electronic hoppings terms	
		if(num_neigh>0)
		{
			i=0;
			if(dist >= (para->NN_lowdis[i]) && dist < (para->NN_highdis[i]))
			{
				if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
				{ 
					t0 = hops[i];
				}
			}
		
			if(dist == 0.0)
			{
				t0 += NN_shifts[i];
			}
			ans=t0;
		}
		
	//Intrinsic SOC - A sublattice
		if(num_neigh>=1)
		{
			i=1;
			
			if(dist >= (para->NN_lowdis[i]) && dist < (para->NN_highdis[i]))
			{
				if( (aDeviceCell->siteinfo)[a][1]==0 && (bDeviceCell->siteinfo)[b][1]==0)
				{
					if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
					{ 
						//printf("SOC! %d, %d %lf\n");
						t0 = possign*SOCsign*hops[i];
					}
				}
			}
		
			if(dist == 0.0)
			{
				t0 += NN_shifts[i];
			}
			ans=t0;
		}
		
	//Intrinsic SOC - B sublattice	
		if(num_neigh>=2)
		{
			i=2;
			
			if(dist >= (para->NN_lowdis[i]) && dist < (para->NN_highdis[i]))
			{
				if( (aDeviceCell->siteinfo)[a][1]==1 && (bDeviceCell->siteinfo)[b][1]==1)
				{
					if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
					{ 
						t0 = possign*SOCsign*hops[i];
					}
				}
			}
		
			if(dist == 0.0)
			{
				t0 += NN_shifts[i];
			}
			ans=t0;
		}
		
		
		
	
	

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
		
			//Nearest neighbour electronic hoppings terms	
			if(num_neigh>0)
			{
				i=0;
			
				if(distp > (para->NN_lowdis[i]) && distp < (para->NN_highdis[i]))
				{
					if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
					{
						t0 = hops[i];
					}
				}
			}
			
			
			possign=0;
			
			//check clockwise/anticlockwise ness for 2NNs
			if(distp > (para->NN_lowdis[1]) && distp < (para->NN_highdis[1]))
			{
				k=0; j=0; i=0;
				
				
				while(j!=1 && i< (acnxp->site_cnxnum)[a])
				{
					//is this neighbour an NN of a AND b?
					midx=(aDeviceCell->pos)[(acnxp->site_cnx)[a][i]][0];
					midy=(aDeviceCell->pos)[(acnxp->site_cnx)[a][i]][1];
					dista = sqrt(pow(midx-x1, 2.0) + pow(midy-y1, 2.0));
					distb = sqrt(pow(midx-x2, 2.0) + pow(midy-y2p, 2.0));
					
					if(dista >= (para->NN_lowdis[0]) && dista < (para->NN_highdis[0]) && distb >= (para->NN_lowdis[0]) && distb < (para->NN_highdis[0]))
					{
						j=1;
					}
					else
						i++;
				}
				
				i=0;
				
				//if not, check if an NN of b is an NN of a
				while(j!=1 && i< (bcnxp->site_cnxnum)[b])
				{
					//is this neighbour an NN of a AND b?
					midx=(bDeviceCell->pos)[(bcnxp->site_cnx)[b][i]][0] + bshifts[0];
					midy=(bDeviceCell->pos)[(bcnxp->site_cnx)[b][i]][1] + bshifts[1] + ycelldist;
					dista = sqrt(pow(midx-x1, 2.0) + pow(midy-y1, 2.0));
					distb = sqrt(pow(midx-x2, 2.0) + pow(midy-y2p, 2.0));
					
					if(dista >= (para->NN_lowdis[0]) && dista < (para->NN_highdis[0]) && distb >= (para->NN_lowdis[0]) && distb < (para->NN_highdis[0]))
					{
						j=1;
					}
					else
					{
						i++;
					}
				}
				
				if(j==1)
				{
				//v1=vector from b to midpoint, v2 from midpoint to a
					v11= midx-x2; v12=midy-y2p;
					v21= x1-midx; v22=y1-midy;
					
					crsp = v11*v22 - v12*v21;
					
					if(crsp > 0)
						possign = 1;
					if(crsp < 0)
						possign = -1;
				}
			
			//Intrinsic SOC - A sublattice
				if(num_neigh>=1)
				{
					i=1;
					if( (aDeviceCell->siteinfo)[a][1]==0 && (bDeviceCell->siteinfo)[b][1]==0)
					{
						if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
						{
							t0 = possign*SOCsign*hops[i];
						}
					}
					
				}
			//Intrinsic SOC - B sublattice
				if(num_neigh>=2)
				{
					i=2;
					
					if( (aDeviceCell->siteinfo)[a][1]==1 && (bDeviceCell->siteinfo)[b][1]==1)
					{
						if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
						{
							t0 = possign*SOCsign*hops[i];
						}
					}
					
				}
			
			
			}
			
			
		ans+=t0*cexp(-I*kpar);	
			
	
		
		
		//check if b is in cell below
		y2p = y2 - ycelldist;
		distp= sqrt(pow(x2-x1, 2.0) + pow(y2p-y1, 2.0));
		

		
		t0=0.0;
		
			//Nearest neighbour electronic hoppings terms	
			if(num_neigh>0)
			{
				i=0;
			
				if(distp > (para->NN_lowdis[i]) && distp < (para->NN_highdis[i]))
				{
					if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
					{
						t0 = hops[i];
					}
				}
			}
			
			
			possign=0;
			
			//check clockwise/anticlockwise ness for 2NNs
			if(distp > (para->NN_lowdis[1]) && distp < (para->NN_highdis[1]))
			{
				k=0; j=0; i=0;
				
				
				while(j!=1 && i< (acnxp->site_cnxnum)[a])
				{
					//is this neighbour an NN of a AND b?
					midx=(aDeviceCell->pos)[(acnxp->site_cnx)[a][i]][0];
					midy=(aDeviceCell->pos)[(acnxp->site_cnx)[a][i]][1];
					dista = sqrt(pow(midx-x1, 2.0) + pow(midy-y1, 2.0));
					distb = sqrt(pow(midx-x2, 2.0) + pow(midy-y2p, 2.0));
					
					if(dista >= (para->NN_lowdis[0]) && dista < (para->NN_highdis[0]) && distb >= (para->NN_lowdis[0]) && distb < (para->NN_highdis[0]))
					{
						j=1;
					}
					else
						i++;
				}
				
				i=0;
				
				//if not, check if an NN of b is an NN of a
				while(j!=1 && i< (bcnxp->site_cnxnum)[b])
				{
					//is this neighbour an NN of a AND b?
					midx=(bDeviceCell->pos)[(bcnxp->site_cnx)[b][i]][0] + bshifts[0];
					midy=(bDeviceCell->pos)[(bcnxp->site_cnx)[b][i]][1] + bshifts[1] - ycelldist;
					dista = sqrt(pow(midx-x1, 2.0) + pow(midy-y1, 2.0));
					distb = sqrt(pow(midx-x2, 2.0) + pow(midy-y2p, 2.0));
					
					if(dista >= (para->NN_lowdis[0]) && dista < (para->NN_highdis[0]) && distb >= (para->NN_lowdis[0]) && distb < (para->NN_highdis[0]))
					{
						j=1;
					}
					else
					{
						i++;
					}
				}
				
				if(j==1)
				{
				//v1=vector from b to midpoint, v2 from midpoint to a
					v11= midx-x2; v12=midy-y2p;
					v21= x1-midx; v22=y1-midy;
					
					crsp = v11*v22 - v12*v21;
					
					if(crsp > 0)
						possign = 1;
					if(crsp < 0)
						possign = -1;
				}
			
				//Intrinsic SOC - A sublattice
				if(num_neigh>=1)
				{
					i=1;
					if( (aDeviceCell->siteinfo)[a][1]==0 && (bDeviceCell->siteinfo)[b][1]==0)
					{
						if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
						{
							t0 = possign*SOCsign*hops[i];
						}
					}
					
				}
				//Intrinsic SOC - B sublattice
				if(num_neigh>=2)
				{
					i=2;
					
					if( (aDeviceCell->siteinfo)[a][1]==1 && (bDeviceCell->siteinfo)[b][1]==1)
					{
						if(zdiff>= NN_zmin[i] && zdiff< NN_zmax[i])
						{
							t0 = possign*SOCsign*hops[i];
						}
					}
					
				}
			
			
			}
			
		
		ans+=t0*cexp(I*kpar); 
	
	}
	
  }
  
  free(possible);
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
  int spinA = para->spinA;
  int spinB = para->spinB;
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
  
  
	//this hopping routine is zero for elements off-diagonal in spin-space
	if(spinA!=spinB)
	{
		ans=0.0;
	}
	
	else
	{
		
		
  
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
	  
	  
// 	    if(leadloop==0)
// 	  {
// 	    printf("DIM %d\n", dim1);
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
	}
      }
    free(bshifts);
}
  
 //general routine to prepare lead cells for Rubio-esque routine
//takes positions, cell sep. vector and hopping rule
//returns unit cell ginv, V12 and V21
//magnetic version
//returns spin channel of choice  (spin_channel = 0 or 1)
//or both in one bigger matrix (spin_channel = 2)
void lead_prep2_mag(double _Complex En, RectRedux *LeadCell, int leadindex, rib_lead_para *params, double _Complex **ginv, double _Complex **V12, double _Complex **V21, int spin_channel)
{
    int i, j, k, s1, s2;
    
    int dim = *(LeadCell->Nrem);
    
    hoppingfunc *hopfn = (hoppingfunc *)(params->hopfn);
    gen_hop_params *hoppingparams = (params->hoppara);
    double *bshifts = createDoubleArray(3);
    
	double **spin_pots = (LeadCell->spin_pots);
	double *const_spin_pots= (LeadCell->const_spin_pots);
	int are_spin_pots=(LeadCell->are_spin_pots);
    
	if(spin_channel == 0 || spin_channel ==1)
	{
		(hoppingparams->spinA)=spin_channel;
		(hoppingparams->spinB)=spin_channel;
		
		//g00
		for(i=0; i <dim; i++)
		{
			
			ginv[i][i] = En - (LeadCell->site_pots[i]) - (hopfn)(LeadCell, LeadCell, i, i, bshifts,  (params->hoppara) ) ;
			
			if(are_spin_pots == 1)
			{
				if(spin_channel ==0)
				{
					ginv[i][i] -= (const_spin_pots[2] + spin_pots[i][2]);
				}
				if(spin_channel ==1)
				{
					ginv[i][i] -= (-const_spin_pots[2] - spin_pots[i][2]);
				}
			}
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
			}
		}
	}
	
	if(spin_channel == 2)
	{
		
	   for(s1=0; s1<2; s1++)
	   {
		(hoppingparams->spinA)=s1;
		
		//common for both spins
		for(i=0; i <dim; i++)
		{
			ginv[s1*dim+i][s1*dim+i] = En- (LeadCell->site_pots[i]);
		}
		
		
		for(s2=0; s2<2; s2++)
		{
			(hoppingparams->spinB)=s2;
			
		
			
			for(i=0; i<3; i++)
				bshifts[i] = 0.0;
			
			//g00
			for(i=0; i <dim; i++)
			{
				
				ginv[s1*dim+i][s2*dim+i] -= (hopfn)(LeadCell, LeadCell, i, i, bshifts, (params->hoppara) ) ;
				
				
				//spin potentials
				if(are_spin_pots == 1)
				{
					if(s1 ==0 && s2==0)
					{
						ginv[s1*dim+i][s2*dim+i] -= (const_spin_pots[2] + spin_pots[i][2]);
					}
					
					if(s1==0 && s2==1)
					{
						ginv[s1*dim+i][s2*dim+i] -= (const_spin_pots[0] + spin_pots[i][0] - I*(const_spin_pots[1] + spin_pots[i][1]));
					}
					
					if(s1==1 && s2==0)
					{
						ginv[s1*dim+i][s2*dim+i] -= (const_spin_pots[0] + spin_pots[i][0] + I*(const_spin_pots[1] + spin_pots[i][1]));
					}
					
					
					if(s1 ==1 && s2==1)
					{
						ginv[s1*dim+i][s2*dim+i] -= (-const_spin_pots[2] - spin_pots[i][2]);
					}
				}
				
				
				for(j=0; j<dim; j++)
				{
					if(j!=i)
					{
						ginv[s1*dim+i][s2*dim+j] = - (hopfn)(LeadCell, LeadCell, i, j, bshifts, (params->hoppara) );
					}
				}
				
			}

			//Vs
			for(i=0; i<3; i++)
				bshifts[i] = (params->shift_vec)[i];
			for(i=0; i <dim; i++)
			{
				for(j=0; j<dim; j++)
				{
					V12[s1*dim+i][s2*dim+j] =  (hopfn)(LeadCell, LeadCell, i, j, bshifts, (params->hoppara) );
				}
				
			}
			
			for(i=0; i<3; i++)
				bshifts[i] = -(params->shift_vec)[i];
			
			for(i=0; i <dim; i++)
			{
				for(j=0; j<dim; j++)
				{
					V21[s1*dim+i][s2*dim+j] =  (hopfn)(LeadCell, LeadCell, i, j, bshifts, (params->hoppara) );
				}
			}
			
		}
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



void multipleCustomLeads (double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, lead_para *params, double _Complex **Sigma)
{
  
    int leadloop, dim1, dim1a,  dimcounta=0, dimcountb=0, lcount;
    double _Complex **smallSigma, **temp1;

    int num_leads = (cellinfo->num_leads);
    gen_hop_params *hopp = (params->hoppara);
    double sep;
    double _Complex *hops = (hopp->hops);
    int i, j, k;
    int iprime, jprime;
    int dim2, dim2a;
    int spindep = DeviceCell->spindep;
    int leadspindep;
    int cell1dim = (cellinfo->cell1dim);
    
    multiple_para *multiple;
	singleleadfunction *singlefn;
    void * singleparams;
    
    for(leadloop=0; leadloop < num_leads; leadloop++)
    {
		leadspindep = (LeadCells[leadloop]->spindep);
		multiple = (multiple_para *)(params-> multiple)[leadloop];
		singlefn = (singleleadfunction *)(multiple -> indiv_lead_fn);
		singleparams = (multiple -> indiv_lead_para);
		
            dim1a = (cellinfo->lead_dims)[leadloop];

            
	    if(spindep == 0)
	    {
		    dim2 =  dim1a;
		    dim2a = dim1a;
	    }
	    if(spindep == 1)
	    {
		dim2 = 2*dim1a;
		dim2a = 2*dim1a;
	    }
	    if(spindep == 2)
	    {
		dim2 = 2*dim1a;
		dim2a = dim1a;
	    }
	    
            smallSigma = createNonSquareMatrix(dim2, dim2a);
            
           //external function call of multiple.indiv_lead_fn to calculate sigma
	    (singlefn)(leadloop, En, DeviceCell, LeadCells, cellinfo, singleparams, smallSigma);
            
            MatrixCopyPart(smallSigma, Sigma, 0, 0, dimcounta, dimcounta, dim1a, dim1a);
	    
	    if(spindep==1)
	    {
		MatrixCopyPart(smallSigma, Sigma, 0, dim1a, dimcounta, cell1dim + dimcounta, dim1a, dim1a);
		MatrixCopyPart(smallSigma, Sigma, dim1a, dim1a, cell1dim + dimcounta, cell1dim + dimcounta, dim1a, dim1a);
		MatrixCopyPart(smallSigma, Sigma, dim1a, 0, cell1dim + dimcounta, dimcounta, dim1a, dim1a);
	    }
	    
	     
	    if(spindep==2)
	    {
		MatrixCopyPart(smallSigma, Sigma, dim1a, 0, cell1dim + dimcounta, dimcounta, dim1a, dim1a);
	    }
	    
            FreeMatrix(smallSigma); 
            
            dimcounta += dim1a;
	  
    }
    
 

}


//generate Sigma for a single Ribbon type lead
void singleRibbonLead (int leadnum, double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, void *params, double _Complex **Sigma)
{
    int leadloop, dim1, dim1a,  lcount, dimcounta=0, dim2, dim2a;
    double _Complex **ginv, **V12, **V21, **g00, **SL, **SR, **VLD, **VDL, **smallSigma, **temp1, **temp2, **temp3;
    double elemerr=1.0e-15;
	
    rib_lead_para *ribpara = (rib_lead_para *)params;
    
    int i, j, k, s1, s2;
  
    int spindep = (DeviceCell->spindep);
    int leadspindep = (LeadCells[leadnum] -> spindep);
    
    double *bshifts0 = createDoubleArray(3);
    hoppingfunc *hopfn = (hoppingfunc *)(ribpara->hopfn);
    gen_hop_params *hoppingparams = (gen_hop_params *)(ribpara->hoppara);
    
	for(leadloop=0; leadloop < leadnum; leadloop++)
	{
		dim1a = (cellinfo->lead_dims)[leadloop];
		dimcounta += dim1a;
	}
	
	dim1 = *(LeadCells[leadnum]->Nrem);
	dim1a = (cellinfo->lead_dims)[leadnum];
	
	if(spindep==0)
	{
		dim2 = dim1;
		dim2a = dim1a;
	}
	if(spindep==1 || spindep==2)
	{
		dim2 = 2*dim1;
		dim2a = 2*dim1a;
	}
		
	
      //generate leads SGFs
      if(spindep == 0 || spindep == 1)
      {
	  
	  ginv = createSquareMatrix(dim2);
	  V12 = createSquareMatrix(dim2);
	  V21 = createSquareMatrix(dim2);
	  g00 = createSquareMatrix(dim2);
	  SL = createSquareMatrix(dim2);

	  //generate the info required for Rubio method
	  
		if(spindep==0)
			lead_prep2(En, LeadCells[leadnum], leadnum, ribpara, ginv, V12, V21);
		if(spindep==1)
			lead_prep2_mag(En, LeadCells[leadnum], leadnum, ribpara, ginv, V12, V21, 2);
	  
	  
	 //printEMatrix(ginv, dim2);

	  InvertMatrixGSL(ginv, g00, dim2);
	  //printEMatrix(V21, dim2);
	  RubioSGF(SL, g00, V12, V21, dim2, &lcount, elemerr*dim2*dim2);
	  FreeMatrix(ginv); FreeMatrix(V12); FreeMatrix(V21); FreeMatrix (g00);
	//printEMatrix(SL, dim2);

    //connections of leads to device
	  //leads are connected to the device using the hopping rules of the lead region, 
	  //not the device region (if different)
	  
	  
	  VLD = createNonSquareMatrix(dim2, dim2a);
	  VDL = createNonSquareMatrix(dim2a, dim2);
	  
	  
	  for(s1=0; s1<spindep+1; s1++)
 	  {
		  (hoppingparams->spinA) = s1;
 		  for(s2=0; s2<spindep+1; s2++)
		  {
			(hoppingparams->spinB) = s2;
			for(i=0; i <dim1; i++)
			{
				for(j=0; j<dim1a; j++)
				{
					k=(cellinfo->lead_sites)[dimcounta + j];
					VLD[s1*dim1 + i][s2*dim1a + j] =  (hopfn)(LeadCells[leadnum], DeviceCell, i, k, bshifts0, hoppingparams );
					VDL[s1*dim1a +j][s2*dim1 +i] =  (hopfn)(DeviceCell, LeadCells[leadnum], k, i, bshifts0, hoppingparams );
				}
				
			}
		  }
	  }
	  temp1 = createNonSquareMatrix(dim2a, dim2);
	  MatrixMultNS(VDL, SL, temp1, dim2a, dim2, dim2);
	  MatrixMultNS(temp1, VLD, Sigma, dim2a, dim2, dim2a);
	  FreeMatrix(temp1); FreeMatrix(VLD); FreeMatrix(VDL);
      }
      
      if(spindep == 2)
      {
		for(s1=0; s1<2; s1++)
		{
			(hoppingparams->spinA) = s1;
			(hoppingparams->spinB) = s1;
	      
			ginv = createSquareMatrix(dim2);
			V12 = createSquareMatrix(dim2);
			V21 = createSquareMatrix(dim2);
			g00 = createSquareMatrix(dim2);
			SL = createSquareMatrix(dim2);

		//generate the info required for Rubio method
		
			lead_prep2_mag(En, LeadCells[leadnum], leadnum, ribpara, ginv, V12, V21, s1);
		
		

			InvertMatrixGSL(ginv, g00, dim2);
			RubioSGF(SL, g00, V12, V21, dim2, &lcount, elemerr*dim2*dim2);
			FreeMatrix(ginv); FreeMatrix(V12); FreeMatrix(V21); FreeMatrix (g00);
		

			//connections of leads to device
			//leads are connected to the device using the hopping rules of the lead region, 
			//not the device region (if different)
			
			
			VLD = createNonSquareMatrix(dim2, dim2a);
			VDL = createNonSquareMatrix(dim2a, dim2);
			
			
			
			for(i=0; i <dim1; i++)
			{
				for(j=0; j<dim1a; j++)
				{
					k=(cellinfo->lead_sites)[dimcounta + j];
					VLD[i][j] =  (hopfn)(LeadCells[leadnum], DeviceCell, i, k, bshifts0, hoppingparams );
					VDL[j][i] =  (hopfn)(DeviceCell, LeadCells[leadnum], k, i, bshifts0, hoppingparams );
				}
				
			}
			
			
			temp1 = createNonSquareMatrix(dim2a, dim2);
			temp2=createSquareMatrix(dim2a);
			MatrixMultNS(VDL, SL, temp1, dim2a, dim2, dim2);
			MatrixMultNS(temp1, VLD, temp2, dim2a, dim2, dim2a);
			
			
			MatrixCopyPart(temp2, Sigma, 0, 0, s1*dim1a, 0, dim1a, dim1a);
			
			FreeMatrix(temp1); FreeMatrix(VLD); FreeMatrix(VDL); FreeMatrix(temp2);
		
		}
	      
      }
	
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





