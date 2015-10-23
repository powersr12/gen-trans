
//General transport code
//Built from sections of ribbon-asymm code and others
//Accepts devices with RectRedux geometries, but calls specified routines for defining cells and setting hopping values
//This file is for developing routines - finished products should go in separate files to be included in compilation

//DO NOT TRY PERIODIC SYSTEMS WITH ODD RIBBON WIDTHS!

//command line version
#include "test.h"

#define eta 1.0e-4



//remove unneeded

main(int argc, char *argv[])
{
  
	//counters
	  int i, j,k,l;
	
	//length and index array related constants
	  int Ntot, Nrem, length1=2,  length2 = length1*3 , geo=0;
	  length2 = 30;
	
	//processor and energy loop related integers
	  int procs=1, this_proc=0, remainder, Epts_temp=11, en, numfin, Epoints=11, conf_num=0;
	  
	//energy loop related doubles 
	  double realE=0.0, Emax=-0.9, Emin=0.9, Estep=0.04, Emin_temp=0.0;
	  double _Complex En;
	  
	//k loop & periodicity related params
	  int isperiodic=1;
	  double kmin=0.0, kmax=2*M_PI;
	  int kpts=1;
	  double kstep;
	      if(kpts>1)
	      {
		kstep = (kmax-kmin)/(kpts-1);
	      }
	      if (kpts==1)
	      {
		kstep = 0.0;
	      }
	  cnxRulesFn *connectrules;
	    if(isperiodic==0)
	      connectrules = &zzacnn;
	    if(isperiodic==1)
	      connectrules = &zzacnnk;
	    
	  
	  
	  
	   int buffer_rows=0;
	//disorder / system config parameters - for sublattice
	  double suba_conc=1.0, suba_pot=0.062, subb_conc=1.0, subb_pot=-0.062;
	  
	//disorder / system config parameters - for antidots overwrites some of the above...
	  int AD_length=5, lat_width=1, lat_length=2;
	  double AD_rad=1.0;
	  
	//lead info
	  int num_leads=2;
	  
// 	  if(geo==0)
// 	  {
// 	    length1=2*lat_width*AD_length ; length2= 3*AD_length*lat_length + 2*buffer_rows;
// 	  }
// 	  
// 	  if(geo==1)
// 	  {
// 	    length1=6*lat_width*AD_length ; length2= AD_length*lat_length + 2*buffer_rows;
// 	  }

	    

	//misc  
	  char job_name[64] = "test";
	  int output_type=1;   //=0 no structure output, =1 just atoms by sublattice, =2 atoms and connectors
	  int mappings=0;		   
	  double cond2;

	  
	  
	//Read command line input - quit if incorrect input
	  //Maybe replace this with input file?
// 	    if(argc ==17)
// 	    {
// 		    sscanf(argv[1], "%s", job_name);
// 		    sscanf(argv[2], "%d", &geo);
// 		    sscanf(argv[3], "%d", &length1);
// 		    sscanf(argv[4], "%d", &length2);
// 		    sscanf(argv[5], "%d", &buffer_rows);
// 		    
// 		    sscanf(argv[6], "%lf", &suba_conc);
// 		    sscanf(argv[7], "%lf", &suba_pot);
// 		    
// 		    sscanf(argv[8], "%lf", &subb_conc);
// 		    sscanf(argv[9], "%lf", &subb_pot);
// 		    
// 		    sscanf(argv[10], "%lf", &Emin);
// 		    sscanf(argv[11], "%lf", &Emax);
// 		    sscanf(argv[12], "%d", &Epoints);
// 		    sscanf(argv[13], "%d", &mappings); 	//how frequently to make maps (0 = never, 1 = always, 
// 								//2 = every second energy etc, -n every nth energy but only for conf 0)
// 
// 		    sscanf(argv[14], "%d", &conf_num);
// 		    sscanf(argv[15], "%d", &procs);
// 		    sscanf(argv[16], "%d", &this_proc);
// 	    }
// 	    
// 	    else
// 	    {
// 		    exit(1);
// 	    }
//       

	//File I/O variables
	    FILE *output;
	    char filename[160], filename3[160], filename_temp[128], checkname[128], direcname[128], conffile[128], strucfile[128];
            char command[256];
	    
	//Create directory and filenaming convention
	    sprintf(direcname, "res/%s", job_name);
	    sprintf(command, "mkdir -p %s", direcname);
	    system(command);
	    
	    sprintf(filename_temp, "%s.conf%02d", job_name, conf_num); 
	    sprintf(strucfile, "%s/%s.struct", direcname, filename_temp);
	    sprintf(filename, "%s/.%s.part%02d", direcname, filename_temp, this_proc);
	
	//Create main output file
	    output =fopen(filename, "w");
	    fclose(output);

      
      
	//Status files for multiprocessor runs
      
	    //Check file - this keeps track of how many of the processes are complete
		FILE *check;
		sprintf(checkname, "%s/.%s.CHECK", direcname, filename_temp);

		if(this_proc ==0)
		{
			check = fopen(checkname, "w");
			fprintf(check, "%d", 0);
			fclose(check);
			output_type=1;
		}

	      Estep = (Emax - Emin) / (Epoints -1);
	      remainder = (Epoints % procs);
	      
	      if(remainder != 0 && this_proc < remainder)
	      {
		      Epts_temp = ((Epoints - remainder) / procs) + 1;
		      Emin_temp = Emin + this_proc * (Epts_temp) * Estep;
		      
	      }
	      else
	      {
		      Epts_temp = Epoints/procs;
		      Emin_temp = Emin + remainder * (Epts_temp + 1) * Estep + (this_proc - remainder) * (Epts_temp) * Estep ;

 	      }
      


    
      
//------------------------------ START OF THE ACTUAL INTERESTING COMPUTATIONAL NON-I/O BIT -----------------------------------//
      
           	//generate System Geometry

      
	    RectRedux System = {};
	    System.geo=geo;
	    System.length=length1;
	    System.length2=length2;
	    System.Nrem=&Nrem;
	    System.Ntot=&Ntot;
	    
	    
	    RectRedux *LeadCells[num_leads];

	    //timing
	    clock_t time;
	    time = clock();
	    

	    //generate, or read in, disorder configuration
		sprintf(conffile, "%s/.%s", direcname, filename_temp);
		
		if(this_proc==0)
		{
		  genSublatticeDevice (&System, buffer_rows, suba_conc, suba_pot, subb_conc, subb_pot, conf_num, 1, strucfile);
		//  genAntidotDevice (&System, buffer_rows, AD_length, AD_rad, lat_width, lat_length, conf_num, 1, strucfile);

		  
		  //export the disorder configuration for the other processes calculating the same configuration
		  if(procs>0)
		  {
		    exportRectConf(&System, conffile);
		  }
		}
		
		lead_para leadp={};
		
		
		time = clock() - time;
	  printf("#generated geometry in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();
		  
		  if(this_proc > 0)
		  {
		    //sleep(100);
		    importRectConf(&System, length1, length2, conffile);
		  }
		  
	    
	   // genSublatticeDevice (&System, 0, 0.00, -1.0, 0.0, 0.0, conf_num, 1, strucfile);
	    //exportRectConf(&System, "test");
	    double **pos = System.pos;
	    int **siteinfo = System.siteinfo;
	    double *site_pots = System.site_pots;
	    
	    

	  cnxProfile cnxp;
	  cnxp.max_neigh=3;
	  device_connectivity (&System, connectrules, NULL, &cnxp);
	  
		time = clock() - time;
		printf("#made connection profile in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();
   	  //printConnectivity (&System, &cnxp);
	  
	  cellDivision cellinfo;
	  genStartingCell(&System, &cellinfo, 2, NULL);
	  
		time = clock() - time;
		printf("#generated starting cell in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();
	  
	  cellSplitter(&System, &cnxp, &cellinfo);
	  
		time = clock() - time;
		printf("#split cells in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();

// 	  printf("# %d\n", cellinfo.group_cell);
	  
	  simpleTB_params hoppara ={-1.0, isperiodic, 0.0, 0.56, 0.59};
	  

	  peierlsTB_params maghoppara={};
	  maghoppara.t0=-1.0;
	  maghoppara.isperiodic=isperiodic;
	  maghoppara.kpar=0.0;
	  maghoppara.NN_lowdis=0.56;
	  maghoppara.NN_highdis=0.59;
	  maghoppara.gauge=0;
	  maghoppara.Btes=20;
	  
	  int *res = createIntArray(2);
	  res[0] = 1;
	  res[1] = 0; 
	  double **reslimits = createNonSquareDoubleMatrix(2, 2);
	  reslimits[0][0] = pos[0][0];
	  reslimits[0][1] = pos[Ntot-1][0];
	
	  maghoppara.restrics = res;
	  maghoppara.limits = reslimits;
	  
	  double _Complex **Sigma = createSquareMatrix(cellinfo.cell1dim);
	  double _Complex **G00 = createSquareMatrix(cellinfo.cell1dim);

 	 // genDeviceGF(0.2+eta*I, &System, &cnxp, &cellinfo, &simpleTB, &hoppara, 0, G00, NULL, Sigma);
	  
	  // simple2leads (0.2+eta*I, &System, LeadCells, &cellinfo, &leadp, Sigma);

// 	  leadp.hopfn = &simpleTB;
// 	  leadp.hoppara = &hoppara;
	  
	   leadp.hopfn = &peierlsTB;
	  leadp.hoppara = &maghoppara;
// 	  
	  leadp.leadsfn = &simple2leads;
	  
	  trans_params tpara = {};
	  tpara.num_leads=2;
	  tpara.TRsym=1;
	  double **transmissions = createNonSquareDoubleMatrix(num_leads, num_leads);

	  tpara.transmissions = transmissions;
	  		genLeads(&System, LeadCells, 2, 0, &leadp);

	  double kavg;
			
 	for(realE=0.002; realE<0.5; realE+=0.002)
	{
 	// realE=0.8;
	    //genTransmissions(realE+eta*I, &System, LeadCells, &cnxp, &cellinfo, &simpleTB, &hoppara, &leadp, &tpara);
	    kavg=0.0;

	   for(k=0; k<kpts; k++)
	   {
	     maghoppara.kpar = kmin + k*kstep;
	     	     hoppara.kpar = kmin + k*kstep;

	     genTransmissions(realE+eta*I, &System, LeadCells, &cnxp, &cellinfo, &peierlsTB, &maghoppara, &leadp, &tpara);
	     kavg += (transmissions[0][1]/kpts);
// 	     	    printf("%lf	%e	\n", maghoppara.kpar, transmissions[0][1]);

	   }
	   // printf("#%lf	%e	%e	%e	%e\n", realE, transmissions[0][0], transmissions[0][1], transmissions[1][0], transmissions[1][1]);
	    printf("%lf	%e	\n", realE, kavg);

	}
	  
// 	      double *ldoses = createDoubleArray(*(System.Nrem));
// 	      double **conds = createNonSquareDoubleMatrix(*(System.Nrem), 3);
// 
// 	      int *indices = createIntArray(2*System.length*System.length2);
// 	      int **neigh = createNonSquareIntMatrix(2*System.length*System.length2, 3);
	      
// 	      for(en=0; en < Epts_temp; en++)
// 	      {
// 
// 		      realE =  Emin_temp + en*Estep;
// 		      sprintf(filename3, "%s/%s.en_%.3lf", direcname, filename_temp, realE);
// 
// 		      if(mappings == 0)
// 		      {
// 			cond2 = genConduc5(realE + eta*I, &System, -1.0);
// 		      }
// 		      
// 		      if(mappings > 0)
// 		      {
// 			if((en % mappings) == 0)
// 			{
// 			  cond2 = genConduc4(realE + eta*I, &System, ldoses, conds, indices, neigh, -1.0);
// 			  CurrentMapOut(&System, conds, indices, neigh, filename3);
// 			  LDOSMapOut(&System, ldoses, filename3);
// 			}
// 			if((en % mappings) != 0)
// 			{
// 			  cond2 = genConduc5(realE + eta*I, &System, -1.0);
// 			}
// 		      }
// 		      
// 		      if(mappings < 0)
// 		      {
// 			if((en % mappings) == 0 && conf_num == 0)
// 			{
// 			  cond2 = genConduc4(realE + eta*I, &System, ldoses, conds, indices, neigh, -1.0);
// 			  CurrentMapOut(&System, conds, indices, neigh, filename3);
// 			  LDOSMapOut(&System, ldoses, filename3);
// 			}
// 			if((en % mappings) != 0 || conf_num != 0)
// 			{
// 			  cond2 = genConduc5(realE + eta*I, &System, -1.0);
// 			}
// 		      }
// 		      
// 		  //    printf("%lf	%e\n", realE,  genConduc5(realE + eta*I, &System, -1.0));
// 			
// 		      //cond2 =genConduc5(realE + eta*I, &System, -1.0);
// 		      
// 			//    printf("%lf	%e	%e\n", realE, genConduc5(realE + eta*I, &System, -1.0),genConduc4(realE + eta*I, &System, ldoses, conds, indices, neigh, -1.0) );
// 			    
// 		      output =fopen(filename, "a");
// 		      fprintf(output, "%lf	%e\n", realE, cond2);
// 		      fclose(output);
// 
// 		      //printf("%lf	%lf\n", realE, cond2);
// 
// 
//  	      }
	      
	      
// 	      genConduc4(realE + eta*I, &System, ldoses, conds, indices, neigh, -1.0);
// 	      CurrentMapOut(&System, conds, indices, neigh, filename);
// 	      LDOSMapOut(&System, ldoses, filename);



// 	      
// 	      
// 	      
// 	      
// //------------------------------ END OF THE ACTUAL INTERESTING COMPUTATIONAL NON-I/O BIT -----------------------------------//
// 
// 	      
// 	      
// 	      
// 	      
// 	      
	       //define this process as finished, and check how many others are
// 		check = fopen(checkname, "r");
// 		fscanf(check, "%d", &numfin);
// 		fclose(check);
// 
// 		numfin ++;
// 		
// 		
// 		check = fopen(checkname, "w");
// 		fprintf(check, "%d", numfin);
// 		fclose(check);
// 	      
// 		//if all processes finished then merge and sort outputs and purge temp files
// 		if(numfin == procs)
// 		{
// 		    sprintf(filename, "%s/.%s", direcname, filename_temp);
// 		    sprintf(command, "cat %s.part* | sort -n > %s/%s.dat", filename, direcname, filename_temp);
// 		    system(command);
// 		    
// 
// 		    
// 		    
// 		    if(conf_num != 0)
// 		    {
// 		      sprintf(command, "rm %s/.%s.*", direcname, filename_temp);
// 		      system(command);
// 		    }
// 	
// 		  
// 		}
		
}





