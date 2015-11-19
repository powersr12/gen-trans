
//General transport code
//Built from sections of ribbon-asymm code and others
//Accepts devices with RectRedux geometries, but calls specified routines for defining cells and setting hopping values
//This file is for developing routines - finished products should go in separate files to be included in compilation

//DO NOT TRY PERIODIC SYSTEMS WITH ODD RIBBON WIDTHS!

//command line version
#include "test.h"

// #define eta 1.0e-6


//remove unneeded

main(int argc, char *argv[])
{
    
	  //counters
	  int i, j,k,l;
	  
	  //general system inputs and default values
	      char systemtype[32];
	      char geotype[32], peritype[32], leadtype[32], sysinfo[80], loopinfo[64];
	      sprintf(systemtype, "SUBLATTICEPOT");
	      int length1=2, length2=3*length1, geo=0, isperiodic=1, ismagnetic=0;
	      int output_type=1;   //=0 no structure output, =1  atoms / holes
	      double eta = 1.0E-6;
	      
	  
	      //check for command line arguments which vary these
		  for(i=1; i<argc-1; i++)
		  {
		    if(strcmp("-systype", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%s", systemtype);
		    }
		    if(strcmp("-geo", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &geo);
		    }
		    if(strcmp("-length1", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &length1);
		    }
		    if(strcmp("-length2", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &length2);
		    }
		    if(strcmp("-periodic", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &isperiodic);
		    }
		     if(strcmp("-magnetic", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &ismagnetic);
		    }
		     if(strcmp("-output", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &output_type);
		    }
		     if(strcmp("-eta", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &eta);
		    }
		    
		  }
	  
	  
	 
		  if(geo==0)
		  {
		    sprintf(geotype, "ZZ");
		  }
		  if(geo==1)
		  {
		    sprintf(geotype, "AC");
		  }
		  
	  
  
	  generatefn *SysFunction;
	  void *SysPara;
	  void *hopfn;
	  
	  	  if(ismagnetic == 0)
		  {
		    hopfn = &simpleTB;   
		  }
		  else if(ismagnetic == 1)
		  {
		    hopfn = &peierlsTB;   
		  }
		  
	//hopping related
	double t0=-1.0, NNlowdis=0.56, NNhighdis=0.59;
	
	//check for command line arguments which vary these
	    for(i=1; i<argc-1; i++)
	    {
	      if(strcmp("-t0", argv[i]) == 0)
	      {
		  sscanf(argv[i+1], "%lf", &t0);
	      }
	      if(strcmp("-NNlowdis", argv[i]) == 0)
	      {
		  sscanf(argv[i+1], "%lf", &NNlowdis);
	      }
	      if(strcmp("-NNhighdis", argv[i]) == 0)
	      {
		  sscanf(argv[i+1], "%lf", &NNhighdis);
	      }
	    }	  
		  
		    
	
	
	//array related constants
	  int Ntot, Nrem;
	  
	//connectivity related constants
	  int max_neigh=3;
	
	//processor and energy/magnetic loop related integers
	  int procs=1, this_proc=0, remainder, en, numfin, conf_num=0;
	  int loop_pts=101, loop_pts_temp=0;
	  char loop_type[4];
	  sprintf(loop_type, "E"); //= 'E'; //E for energy loop, B for Bfield loop
	  
	  int mappings=0;
	  char job_name[64] = "";	//an addendum to folder names to 

	  
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-numofprocs", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &procs);
		      }
		      if(strcmp("-thisproc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &this_proc);
		      }
		      if(strcmp("-confnum", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &conf_num);
		      }
		      if(strcmp("-looptype", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s", loop_type);
		      }
		      if(strcmp("-looppts", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &loop_pts);
		      }
		      if(strcmp("-mappings", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &mappings);
		      }
		      if(strcmp("-jobname", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s", job_name);
		      }
		    }

	  
	//energy and magfield loop related doubles 
	  double realE=0.0, Emax=0.9, Emin=-0.9;
	  double Bfield=0.0, Bmax=20.0, Bmin=0.0;
	  double _Complex En;
	  
	  
	  //check for command line arguments which vary these
		  for(i=1; i<argc-1; i++)
		  {
		    if(strcmp("-Emin", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &Emin);
		    }
		    if(strcmp("-Emax", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &Emax);
		    }
		    if(strcmp("-Bmin", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &Bmin);
		    }
		    if(strcmp("-Bmax", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &Bmax);
		    }
		  }
	  

	  
	//k loop & periodicity related params
	  
	  double kmin=0.0, kmax=2*M_PI;
	  int kpts=1;
	  double kstep;
	  
		  //check for command line arguments which vary these
		  for(i=1; i<argc-1; i++)
		  {
		    if(strcmp("-kpts", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &kpts);
		    }
		    if(strcmp("-kmin", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &kmin);
		    }
		  }
		  
		  
	  
	    if(kpts>1)
	    {
	      kstep = (kmax-kmin)/(kpts-1);
	    }
	    if (kpts==1)
	    {
	      kstep = 0.0;
	    }
	    
	    cnxRulesFn *connectrules;
	    graph_conn_para cnxpara;
	    cnxpara.conn_sep_thresh_min = NNlowdis;
	    cnxpara.conn_sep_thresh_max = NNhighdis;
	    

	    
	    
	      if(isperiodic==0)
	      {
		connectrules = &graph_conn_sep;
		cnxpara.periodic = 0;
		sprintf(peritype, "RIBBON");
		kpts = 1;

	      }
	      if(isperiodic==1)
	      {
		connectrules = &zzacnnk;
		cnxpara.periodic = 1;
		sprintf(peritype, "PERIODIC");
	      }
		
	  
	  
	    
	  
	   int buffer_rows=0;
	//disorder / system config default parameters - for sublattice
	  double suba_conc=1.0, suba_pot=0.062, subb_conc=1.0, subb_pot=-0.062;
	  int xory=0;
	  double interface_width=0.0, interface_position;
	  
	  
	  
	    
	  
	  
	//disorder / system config parameters - for antidots overwrites some of the above...
	  int AD_length=5, AD_length2=5, lat_width=1, lat_length=2;
	  double AD_rad=1.0, xyfluc=0.0, radfluc=0.0;
	  char latgeo[24], dotgeo[24];
	  sprintf(latgeo, "trig");
	  sprintf(dotgeo, "circ");
	  
	//lead info
	  int num_leads=2;
	  
	  
	  //specific device settings, generation functions & parameters
		subl_para sublp = {};
		if(strcmp("SUBLATTICEPOT", systemtype) == 0)
		{	
		    //default values
		    
		    sublp.buffer_rows = buffer_rows;
		    sublp.a_conc = suba_conc;
		    sublp.a_pot = suba_pot;
		    sublp.b_conc = subb_conc;
		    sublp.b_pot = subb_pot;
		    sublp.seed = conf_num;

		    
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-subaconc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublp.a_conc));
		      }
		      if(strcmp("-subapot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublp.a_pot));
		      }
		      if(strcmp("-subbconc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublp.b_conc));
		      }
		      if(strcmp("-subbpot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublp.b_pot));
		      }
		      if(strcmp("-bufferrows", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(sublp.buffer_rows));
		      }
		      
		    }
		    
		    //set functions and params for use below
		    SysFunction = &genSublatticeDevice;
		    SysPara = &sublp;
		    
		    //set filename info - what to put in filename from these params
		    sprintf(sysinfo, "L2_%d_BUF_%d_SUBA_%.2lfx%.3lf_SUBB_%.2lfx%.3lf", length2, buffer_rows, (sublp.a_conc), (sublp.a_pot),(sublp.b_conc), (sublp.b_pot)); 
		}
		
		
		subint_para subintp = {};
		if(strcmp("SUBLATTICEINT", systemtype) == 0)
		{	
		    //default values
		    
		    subintp.buffer_rows = buffer_rows;
		    subintp.a_conc1 = suba_conc;
		    subintp.a_pot1 = suba_pot;
		    subintp.b_conc1 = subb_conc;
		    subintp.b_pot1 = subb_pot;
		    
		    subintp.xory = xory;
		    
		    
		    subintp.seed = conf_num;
		    
		

		    
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-subaconc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(subintp.a_conc1));
		      }
		      if(strcmp("-subapot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(subintp.a_pot1));
		      }
		      if(strcmp("-subbconc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(subintp.b_conc1));
		      }
		      if(strcmp("-subbpot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(subintp.b_pot1));
		      }
		       if(strcmp("-intxory", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(subintp.xory));
		      }
		       if(strcmp("-intwidth", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(subintp.int_width));
		      }
		      
		      if(strcmp("-bufferrows", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(subintp.buffer_rows));
		      }
		      
		    }
		    
			//default values of region 2 parameters are swaps of region 1
			subintp.a_conc2 = subintp.b_conc1;
			subintp.a_pot2 = subintp.b_pot1;
			subintp.b_pot2 = subintp.a_pot1;
			subintp.b_conc2 = subintp.a_conc1;
			
			//default interface positions
			if(geo==0)
			{
			  if(subintp.xory==0)
			    interface_position = (int)(length2/2) *1.0 - 0.5;
			  
			  if(subintp.xory==1)
			    interface_position = (length1*sqrt(3))/4;

			}
			
			if(geo==1)
			{
			  if(subintp.xory==0)
			    interface_position = (int)(length2/2) *sqrt(3) - 1/(2*sqrt(3));
			  
			  if(subintp.xory==1)
			    interface_position = (length1 - 1)/4;
			}
			
			subintp.int_pos = interface_position;		
			
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-subaconc2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(subintp.a_conc2));
		      }
		      if(strcmp("-subapot2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(subintp.a_pot2));
		      }
		      if(strcmp("-subbconc2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(subintp.b_conc2));
		      }
		      if(strcmp("-subbpot2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(subintp.b_pot2));
		      }
		      
		        if(strcmp("-intpos", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(subintp.int_pos));
		      }
		    }
		    
		    
		    
		    
		    
		    //set functions and params for use below
		    SysFunction = &genSublatticeInterface;
		    SysPara = &subintp;
		    
		    //set filename info - what to put in filename from these params
		    sprintf(sysinfo, "L2_%d_BUF_%d_SUBA_%.2lfx%.3lf_SUBB_%.2lfx%.3lf", length2, buffer_rows, (sublp.a_conc), (sublp.a_pot),(sublp.b_conc), (sublp.b_pot)); 
		}
		
		
		
		adot_para adotp = {};
		char lattice_info[35], dot_info[45];
		if(strcmp("ANTIDOTS", systemtype) == 0)
		{	
		    //default values
		    
		    adotp.buffer_rows = buffer_rows;
		    adotp.AD_length = AD_length;
		    adotp.AD_length2 = AD_length2;
		    adotp.latgeo = latgeo;
		    adotp.AD_rad = AD_rad;
		    adotp.AD_rad2 = AD_rad;
		    adotp.lat_width = lat_width;
		    adotp.lat_length = lat_length;
		    adotp.dotgeo = dotgeo;
		    adotp.seed = conf_num;
		    adotp.isperiodic = isperiodic;
		    adotp.xyfluc = xyfluc;
		    adotp.radfluc = radfluc;
		    
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-ADlength", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(adotp.AD_length));
		      }
		      if(strcmp("-ADlength2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(adotp.AD_length2));
		      }
		      if(strcmp("-ADrad", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(adotp.AD_rad));
			  adotp.AD_rad2 = adotp.AD_rad;
		      }
		       if(strcmp("-ADrad2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(adotp.AD_rad2));
		      }
		      if(strcmp("-latw", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(adotp.lat_width));
		      }
		      if(strcmp("-latl", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(adotp.lat_length));
		      }
		      if(strcmp("-bufferrows", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(adotp.buffer_rows));
		      }
		      if(strcmp("-latgeo", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s",  latgeo);
		      }
		      if(strcmp("-dotgeo", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s", dotgeo);
		      }
		      if(strcmp("-xyfluc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(adotp.xyfluc));
		      }
		      if(strcmp("-radfluc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(adotp.radfluc));
		      }
		      
		    }
		    
		    
		    //antidot system sizes calculated from lattice details
		    if(strcmp("trig", latgeo) == 0)
		    {
		      if(geo==0)
		      {
			length1=2*(adotp.lat_width)*(adotp.AD_length) ; 
			length2= 3*(adotp.AD_length)*(adotp.lat_length) + 2*(adotp.buffer_rows);
		      }
		      
		      if(geo==1)
		      {
			length1=6*(adotp.lat_width)*(adotp.AD_length) ; 
			length2= (adotp.AD_length)*(adotp.lat_length) + 2*(adotp.buffer_rows);
		      }
		      
		      sprintf(lattice_info, "%s_lat_L_%d", (adotp.latgeo),(adotp.AD_length)); 

		    }
		    
		    if(strcmp("rotrig", latgeo) == 0)
		    {
		      if(geo==0)
		      {
			length1= 2*((adotp.AD_length)+1)*(adotp.lat_width); 
			length2=  ((adotp.AD_length) + 1)*(adotp.lat_length) + 2*(adotp.buffer_rows);
		      }
		      
		      if(geo==1)
		      {
			length1=(2*(adotp.AD_length) + 2)*(adotp.lat_width) ; 
			length2= ((adotp.AD_length)+1)*(adotp.lat_length) + 2*(adotp.buffer_rows);
		      }
		      
		      sprintf(lattice_info, "%s_lat_L_%d", (adotp.latgeo),(adotp.AD_length)); 
		    }
		    
		    if(strcmp("rect", latgeo) == 0)
		    {
		      if(geo==0)
		      {
			length1=2*(adotp.lat_width)*(adotp.AD_length) ; 
			length2= (adotp.AD_length2)*(adotp.lat_length) + 2*(adotp.buffer_rows);
		      }
		      
		      if(geo==1)
		      {
			length1=2*(adotp.lat_width)*(adotp.AD_length2) ; 
			length2= (adotp.AD_length)*(adotp.lat_length) + 2*(adotp.buffer_rows);
		      }
		      sprintf(lattice_info, "%s_lat_L_%d_%d", (adotp.latgeo),(adotp.AD_length), (adotp.AD_length2)); 

		    }
		    
		    if(strcmp("circ", dotgeo) == 0 || strcmp("hexAC", dotgeo) == 0 || strcmp("hexZZ", dotgeo) == 0)
		    {
		      sprintf(dot_info, "%s_dot_R_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (adotp.dotgeo), (adotp.AD_rad), (adotp.lat_width),  (adotp.lat_length), (adotp.xyfluc), (adotp.radfluc)); 
		    }
		    
		    if(strcmp("rect", dotgeo) == 0 )
		    {
		      sprintf(dot_info, "%s_dot_R_%.1lf_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (adotp.dotgeo), (adotp.AD_rad), (adotp.AD_rad2), (adotp.lat_width),  (adotp.lat_length), (adotp.xyfluc), (adotp.radfluc)); 
		    }

		    
		    sprintf(sysinfo, "%s_%s", lattice_info, dot_info);
		    
		    //set functions and params for use below
		    SysFunction = &genAntidotDevice;
		    SysPara = &adotp;
		    
		    //set filename info? (circ and triglat to be generalised later!)
		}
		
		
		if(strcmp("CLEAN", systemtype) == 0)
		{	
		    
		    //set functions and params for use below
		    SysFunction = &simpleRibbonGeo;
		    SysPara = NULL;
		    
		    //set filename info?
		    	sprintf(sysinfo, "clean");

		}
		
	  


	  
	  
	//misc  
	  
	  double cond2;


	  //loop info
	      double loop_step, loop_min_temp, loopmax, loopmin;
	      
	      if(strcmp("E", loop_type) == 0)
	      {
		Bfield = Bmin;
		loopmin = Emin;
		loopmax = Emax;
		sprintf(loopinfo, "Eloop_%+.2lf_to_%+.2lf_Bfixed_%+.3lf", Emin, Emax, Bfield);
	      }

	      if(strcmp("B", loop_type) == 0)
	      {
		realE = Emin;
		loopmin = Bmin;
		loopmax = Bmax;
		sprintf(loopinfo, "Bloop_%+.2lf_to_%+.2lf_Efixed_%+.3lf", Bmin, Bmax, realE);

	      }

		loop_step = (loopmax-loopmin)/(loop_pts - 1);
		remainder = (loop_pts % procs);
		
	      
	      if(remainder != 0 && this_proc < remainder)
	      {
		      loop_pts_temp = ((loop_pts - remainder) / procs) + 1;
		      loop_min_temp = loopmin + this_proc * (loop_pts_temp) * loop_step;
		      
	      }
	      else
	      {
		      loop_pts_temp = loop_pts/procs;
		      loop_min_temp = loopmin + remainder * (loop_pts_temp + 1) * loop_step + (this_proc - remainder) * (loop_pts_temp) * loop_step ;

 	      }
	

	//File I/O variables
	    FILE *output;
	    char filename[280], filename3[280], filename_temp[100];
	    char checkname[280], direcname[160], conffile[280], strucfile[280];
            char command[400];
	    
	//Create directory and filenaming convention
	    sprintf(direcname, "../res/%s_%s_%.0e/%s%d_%s", systemtype, peritype, eta, geotype, length1, sysinfo);
	    sprintf(command, "mkdir -p %s", direcname);
	    system(command);
	    printf("# directory: %s\n", direcname);
	    
	    
	    sprintf(filename_temp, "%s_%s.conf%02d", loopinfo, job_name, conf_num); 
	    sprintf(strucfile, "%s/%s.struct", direcname, filename_temp);
	    sprintf(filename, "%s/.%s.part%02d", direcname, filename_temp, this_proc);
	
	//Create main output file
	    output =fopen(filename, "w");
	    fclose(output);

	    
 			  printf("#%s\n", systemtype);
	  	  	  printf("#loop type %s \n", loop_type);

	    
	    
      
	//Status files for multiprocessor runs
      
	    //Check file - this keeps track of how many of the processes are complete
		FILE *check;
		sprintf(checkname, "%s/.%s.CHECK", direcname, filename_temp);

		if(this_proc ==0)
		{
			check = fopen(checkname, "w");
			fprintf(check, "%d", 0);
			fclose(check);
			//output_type=1;
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
		  (SysFunction) ( &System, SysPara, output_type, strucfile);
		  //export the disorder configuration for the other processes calculating the same configuration
		  if(procs>0)
		  {
		    exportRectConf(&System, conffile);
		  }
		}
		

		
		time = clock() - time;
	  printf("#generated geometry in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();
		  
		  if(this_proc > 0)
		  {
		    sleep(8*length1);
		    importRectConf(&System, length1, length2, conffile);
		  }
		  

	    double **pos = System.pos;
	    int **siteinfo = System.siteinfo;
	    double *site_pots = System.site_pots;
	    
	    

	  cnxProfile cnxp;
	  cnxp.max_neigh=max_neigh;
	  device_connectivity (&System, connectrules, &cnxpara, &cnxp);
	  
		time = clock() - time;
		printf("#made connection profile in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();
   	 // printConnectivity (&System, &cnxp);
	  
	  //in theory cell Division could write to lead_para, so its defined here
	  lead_para leadp={};

		
	  cellDivision cellinfo;
	  genStartingCell(&System, &cellinfo, 2, NULL);
	  
		time = clock() - time;
		printf("#generated starting cell in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();
	  
	  cellSplitter(&System, &cnxp, &cellinfo);
	  
		time = clock() - time;
		printf("#split cells in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();


	  //hopping parameters and gauge info
	  gen_hop_params hoppara={};
	  hoppara.t0=t0;
	  hoppara.isperiodic=isperiodic;
	  hoppara.kpar=kmin;
	  hoppara.NN_lowdis=NNlowdis;
	  hoppara.NN_highdis=NNhighdis;
	  hoppara.gauge=0;
	  
	    //mag field cut off
	      int *res = createIntArray(2);
	      res[0] = 1;
	      res[1] = 0; 
	      double **reslimits = createNonSquareDoubleMatrix(2, 2);
	      reslimits[0][0] = pos[0][0];
	      reslimits[0][1] = pos[Ntot-1][0];
	
	  hoppara.restrics = res;
	  hoppara.limits = reslimits;
	  
	  double _Complex **Sigma = createSquareMatrix(cellinfo.cell1dim);
	  double _Complex **G00 = createSquareMatrix(cellinfo.cell1dim);


	  
	  leadp.hopfn = hopfn;
	  leadp.hoppara = &hoppara;
// 	  
	  leadp.leadsfn = &simple2leads;
	  
	  trans_params tpara = {};
	  tpara.num_leads=num_leads;
	  tpara.TRsym=ismagnetic;
	  double **transmissions = createNonSquareDoubleMatrix(num_leads, num_leads);

	  tpara.transmissions = transmissions;
	  genLeads(&System, LeadCells, 2, 0, &leadp);

	  double kavg;
	  
	  
	  for(en=0; en<loop_pts_temp; en++)
	  {
	      if(strcmp("E", loop_type) == 0)
	      {
		realE =  loop_min_temp + en*loop_step;
		sprintf(filename3, "%s/%s.en_%.3lf", direcname, filename_temp, realE);
	      }
	      
	      if(strcmp("B", loop_type) == 0)
	      {
		Bfield =  loop_min_temp + en*loop_step;
		sprintf(filename3, "%s/%s.bf_%.3lf", direcname, filename_temp, Bfield);

	      }
	    
	      hoppara.Btes=Bfield;
	      
	      //needs generalising for multiple transmission cases
	      kavg=0.0;
	      
	      for(k=0; k<kpts; k++)
	      {
		hoppara.kpar = kmin + k*kstep;
		
		genTransmissions(realE+eta*I, &System, LeadCells, &cnxp, &cellinfo, hopfn, &hoppara, &leadp, &tpara);
		kavg += (transmissions[0][1]/kpts);
	      }
	      
	      output =fopen(filename, "a");
	      if(strcmp("E", loop_type) == 0)
	      {
		fprintf(output, "%lf	%e\n", realE, kavg);
		printf( "%lf	%e\n", realE, kavg);

	      }
	      if(strcmp("B", loop_type) == 0)
	      {
		fprintf(output, "%lf	%e\n", Bfield, kavg);
		printf("%lf	%e\n", Bfield, kavg);
	      }
	      
	      fclose(output);

	    
	  }
			
//  	for(realE=0.002; realE<0.5; realE+=0.2)
// 	{
//  	realE=0.8;
// 	    genTransmissions(realE+eta*I, &System, LeadCells, &cnxp, &cellinfo, &simpleTB, &hoppara, &leadp, &tpara);
// 	    kavg=0.0;
// 
// 	   for(k=0; k<kpts; k++)
// 	   {
// 	     	     hoppara.kpar = kmin + k*kstep;
// 
// 	     genTransmissions(realE+eta*I, &System, LeadCells, &cnxp, &cellinfo, hopfn, &hoppara, &leadp, &tpara);
// 	     kavg += (transmissions[0][1]/kpts);
// 	     	    printf("%lf	%e	\n", maghoppara.kpar, transmissions[0][1]);
// 
// 	   }
// 	   printf("#%lf	%e	%e	%e	%e\n", realE, transmissions[0][0], transmissions[0][1], transmissions[1][0], transmissions[1][1]);
// 	    printf("%lf	%e	\n", realE, kavg);
// 
// 	}
	  
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
		check = fopen(checkname, "r");
		fscanf(check, "%d", &numfin);
		fclose(check);

		numfin ++;
		
		
		check = fopen(checkname, "w");
		fprintf(check, "%d", numfin);
		fclose(check);
	      
		//if all processes finished then merge and sort outputs and purge temp files
		if(numfin == procs)
		{
		    sprintf(filename, "%s/.%s", direcname, filename_temp);
 		    sprintf(command, "cat %s.part* | sort -n > %s/%s.dat", filename, direcname, filename_temp);
		    system(command);
		    

		    
		    
		    if(conf_num != 0)
		    {
		      sprintf(command, "rm %s/.%s.*", direcname, filename_temp);
		      system(command);
		    }
	
		  
		}
		
}





