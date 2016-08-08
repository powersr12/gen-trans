
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
	      char geotype[32], peritype[40], leadtype[32], sysinfo[120], loopinfo[80], disinfo[40];
	      sprintf(systemtype, "SUBLATTICEPOT");
	      int length1=2, length2=3*length1, geo=0, isperiodic=1, ismagnetic=0;
	      int makebands = 0, unfold=0, project=0, kxpoints=51, bandsonly=0, bandsminset=0, bandsmaxset=100000;
	      
	      int output_type=1;   //=0 no structure output, =1  atoms / holes
	      double eta = 1.0E-6;
	      int ishallbar=0;
	      int nngm = 1;	//default assumption is Nearest Neighbour graphene hopping parameters and lattice, 
				//changing this to 2 or 3 uses 2NN or 3NN models
				//future change - making this 0 should allow custom parameterisation somehow
	      
	      
	      int magsetup = 0;	//this indicates whether the magnetic field is present in the leads or not
				//(if the system is magnetic)
				//0 means no field in leads
				//1 means field everywhere
				//the correct gauge is then chosen depending on whether we have
				//a simple system or a hall bar
	      int gauge = 0; 	//default gauge choice, phase is along x
	      int potdis =0;	//is there an additional, random potential disorder in the system?
	      int splitgen =0;	//if =1, splits the system generation and calculation
	  
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
		    if(strcmp("-makebands", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &makebands);
		    }
		    if(strcmp("-unfold", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &unfold);
		    }
		    if(strcmp("-project", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &project);
		    }
		    if(strcmp("-kxpoints", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &kxpoints);
		    }
		    if(strcmp("-bandsonly", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &bandsonly);
		    }
		    if(strcmp("-bandsminset", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &bandsminset);
		    }
		      if(strcmp("-bandsmaxset", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &bandsmaxset);
		    }
		     if(strcmp("-output", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &output_type);
		    }
		     if(strcmp("-eta", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &eta);
		    }
		     if(strcmp("-hallbar", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &ishallbar);
		    }
		    if(strcmp("-NNmodel", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &nngm);
		    }
		    if(strcmp("-magsetup", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &magsetup);
		    }
		    if(strcmp("-potdis", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &potdis);
		    }
		  }
	  
	  
	 
		  if(geo==0)
		  {
		    sprintf(geotype, "ZZ");
		    if(nngm==2)
		    {
		      sprintf(geotype, "ZZ2N_");
		    }
		    if(nngm==3)
		    {
		      sprintf(geotype, "ZZ3N_");
		    }
		    
		  }
		  if(geo==1)
		  {
		    sprintf(geotype, "AC");
		    
		      if(nngm==2)
		      {
			sprintf(geotype, "AC2N_");
		      }
		      if(nngm==3)
		      {
			sprintf(geotype, "AC3N_");
		      }
		    
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
	int NNN=1;
	double t0=-1.0, NNlowdis=0.56, NNhighdis=0.59, onsitec=0.0;
	
	//another one for onsites should be added for 2NN
	
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
	    

	    double t1=0.0740741*t0, NNlowdis1=0.98, NNhighdis1=1.02, onsitec1=3*t1;
	    double t2=0.0666667*t0, NNlowdis2=1.13, NNhighdis2=1.17, onsitec2=0.0; 
	    //these can possibly be adjusted by command line later
		    
	
	
	//array related constants
	  int Ntot, Nrem;
	  
	//connectivity related constants
	  int max_neigh=3;  //for graphene NN_shifts
	  if(nngm == 2)
	    max_neigh = 9;
	  if(nngm == 3)
	    max_neigh = 12;
	  
	
	//processor and energy/magnetic loop related integers
	  int procs=1, this_proc=0, remainder, en, numfin, conf_num=0;
	  int loop_pts=101, loop_pts_temp=0;
	  char loop_type[4];
	  sprintf(loop_type, "E"); //= 'E'; //E for energy loop, B for Bfield loop
	  
	  //how often to map, map this particular energy/bfield/loop variable, mapping mode: ldos, currents, all
	  int mappings=0, mapnow=0, mapmode=0;
	  int map_all_leads=1, map_lead=0;
	  
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
		      if(strcmp("-mapmode", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &mapmode);
		      }
		      if(strcmp("-mapall", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &map_all_leads);
		      }
		      if(strcmp("-maplead", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &map_lead);
		      }
		      if(strcmp("-splitgen", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &splitgen);
		      }
		      	      
		      if(strcmp("-jobname", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s", job_name);
		      }
		    }

		    int buffer_rows=0;
		    for(i=1; i<argc-1; i++)
		    {
			if(strcmp("-bufferrows", argv[i]) == 0)
			{
			    sscanf(argv[i+1], "%d", &buffer_rows);
			}
		    }
		    
		    
	//energy and magfield loop related doubles 
	  double realE=0.0, Emax=0.9, Emin=-0.9;
	  double Bfield=0.0, Bmax=20.0, Bmin=0.0;
	  double VG=0.0, VGmin = 0.0, VGmax=10.0;
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
		    if(strcmp("-VGmin", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &VGmin);
		    }
		    if(strcmp("-VGmax", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &VGmax);
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
	    
	    if(nngm==2)
	    {
	      cnxpara.conn_sep_thresh_max = NNhighdis1;
	    }
	     if(nngm==3)
	    {
	      cnxpara.conn_sep_thresh_max = NNhighdis2;
	    }
	    
	    

	    //lead info
	  int num_leads=2;
	    
	      if(isperiodic==0)
	      {
		connectrules = &graph_conn_sep;
		cnxpara.periodic = 0;
		sprintf(peritype, "RIBBON");
		kpts = 1;

	      }
	      if(isperiodic==1)
	      {
		connectrules = &graph_conn_sep;
		cnxpara.periodic = 1;
		sprintf(peritype, "PERIODIC");
	      }
		
	      //being a hall bar overwrites other settings which may be incorrect!
	      //these settings are used later for standard hall probe placements
		double geo_renorm=2;
		int hall_denom=6, hall_rel_width=1, hall_start=1;
		int hall_second_start = hall_denom - hall_start - hall_rel_width; 
		int ntop=2, nbot=2;
		int hall_num_y_cells=10;
	      
	      	hallbpara hallp={};
	
	      
	     
		
	  
	  
	    
	  
	
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
		    
// 		    printf("#BUFFERS %d %d\n", buffer_rows, adotp.buffer_rows);
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
		
		
		
		int edge_buffer_rows=10, sruns=20;
		double smax=2.0, minper=40.0;
		int vacruns =0;
		double vacprob=0.0;
		
		edgedis_para edgedisp = {};
		if(strcmp("EDGEDIS", systemtype) == 0)
		{	
		    //default values
		    
		    edgedisp.buffer_rows = edge_buffer_rows;
		    edgedisp.sruns = sruns;
		    edgedisp.smax = smax;
		    edgedisp.minper = minper;
		    edgedisp.vacruns = vacruns;
		    edgedisp.vacprob = vacprob;

		    edgedisp.seed = conf_num;

		    
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-ebufrows", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(edgedisp.buffer_rows));
		      }
		      if(strcmp("-sruns", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(edgedisp.sruns));
		      }
		      if(strcmp("-smax", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(edgedisp.smax));
		      }
		      if(strcmp("-minper", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(edgedisp.minper));
		      }
		      if(strcmp("-vacruns", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(edgedisp.vacruns));
		      }
		      if(strcmp("-vacprob", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(edgedisp.vacprob));
		      }
		      
		      
		    }
		    
		    //set functions and params for use below
		    SysFunction = &genEdgeDisorderedDevice;
		    SysPara = &edgedisp;
		    
		    //set filename info - what to put in filename from these params
		    sprintf(sysinfo, "L2_%d_DIS%dx%.2lf_per_%.0lf_vac%dx%.2lf", length2, (edgedisp.sruns), (edgedisp.smax),(edgedisp.minper), (edgedisp.vacruns), (edgedisp.vacprob)); 
		}
		
		
		
		if(strcmp("CLEAN", systemtype) == 0)
		{	
		    
		    //set functions and params for use below
		    SysFunction = &simpleRibbonGeo;
		    SysPara = NULL;
		    
		    //set filename info?
		    	sprintf(sysinfo, "clean_l2_%d", length2);

		}
		
	  


	  
	  
	//misc  
	  
	  double cond2;
	  
	  
	  //potential disorder settings
	      double pdconc=0.0, pddelta=0.0, pdxi=1.0;
	      int pdmap=0;
	      potDis_para dispara = {};
	      
	      if(potdis==1)
	      {
		for(i=1; i<argc-1; i++)
		{
		  if(strcmp("-pdconc", argv[i]) == 0)
		  {
		      sscanf(argv[i+1], "%lf", &pdconc);
		  }
		  if(strcmp("-pddelta", argv[i]) == 0)
		  {
		      sscanf(argv[i+1], "%lf", &pddelta);
		  }
		  if(strcmp("-pdxi", argv[i]) == 0)
		  {
		      sscanf(argv[i+1], "%lf", &pdxi);
		  }
		  if(strcmp("-pdmap", argv[i]) == 0)
		  {
		      sscanf(argv[i+1], "%d", &pdmap);
		  }
		 }
		 
		 dispara.conc = pdconc;
		 dispara.delta = pddelta;
		 dispara.xi = pdxi;
		 dispara.seed = conf_num;
		
		 sprintf(disinfo, "POTDIS_c_%.2lf_d_%.3lf_xi_%.1lf", pdconc, pddelta, pdxi);
	      }
	  
	  
	  
	  
	  //standard hall ssettings
	   if(ishallbar==1)
	      {
		connectrules = &graph_conn_sep;
		isperiodic=0;
		cnxpara.periodic = 0;
		kpts = 1;
		
		 for(i=1; i<argc-1; i++)
		 {
		    if(strcmp("-ntop", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &ntop);
		    }
		    if(strcmp("-nbot", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &nbot);
		    }
		    if(strcmp("-halldenom", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &hall_denom);
		    }
		    if(strcmp("-hallrelwidth", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &hall_rel_width);
		    }
		    if(strcmp("-hallstart", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &hall_start);
		    }
		    hall_second_start = hall_denom - hall_start - hall_rel_width;
		    //possibilities for fancier more probe geometries here
		    if(strcmp("-hallsecstart", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &hall_second_start);
		    }
		    if(strcmp("-hallycells", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &hall_num_y_cells);
		    }
		
		 }
		  num_leads =2 + ntop + nbot;
		  
		  sprintf(peritype, "HALLBAR_%d_%d_rw_%d", ntop, nbot, hall_denom);

		  
		  hallp.num_top_probes=ntop, hallp.num_bot_probes=nbot;
		  hallp.toppx = createIntArray(ntop);
		  hallp.toppw = createIntArray(ntop);
		  hallp.toppc = createIntArray(ntop);
		  hallp.botpx = createIntArray(nbot);
		  hallp.botpw = createIntArray(nbot);
		  hallp.botpc = createIntArray(nbot);
		  
		  //THIS ASSUMES ntop and nbot = 2
		  //NEEDS TO BE GENERALISED IF USED FOR OTHER CASES!
		  hallp.toppx[0] = (int) ((hall_start)*(length2-2*buffer_rows)/hall_denom) + buffer_rows;
		  hallp.toppw[0] = (int) (hall_rel_width*geo_renorm*(length2-2*buffer_rows)/hall_denom) ;
		  hallp.toppx[1] = (int) ((hall_second_start)*(length2-2*buffer_rows)/hall_denom) + buffer_rows;
		  hallp.toppw[1] = (int) (hall_rel_width*geo_renorm*(length2-2*buffer_rows)/hall_denom) ;
		  hallp.botpx[0] = (int) ((hall_start)*(length2-2*buffer_rows)/hall_denom) + buffer_rows;
		  hallp.botpw[0] = (int) (hall_rel_width*geo_renorm*(length2-2*buffer_rows)/hall_denom) ;
		  hallp.botpx[1] = (int) ((hall_second_start)*(length2-2*buffer_rows)/hall_denom) + buffer_rows;
		  hallp.botpw[1] = (int) (hall_rel_width*geo_renorm*(length2-2*buffer_rows)/hall_denom) ;
		  
		  hallp.toppc[0] = hall_num_y_cells; hallp.toppc[1] = hall_num_y_cells;
		  hallp.botpc[0] = hall_num_y_cells; hallp.botpc[1] = hall_num_y_cells;
		
	      }
	  


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
	      
	      
	      double edge_cut=5.0, subs_thick=0.1E-6, subs_epsr=3.9;
	      int vgtype=0; 
	      char vgtname[20];
	      if(strcmp("VG", loop_type) == 0)
	      {
		Bfield = Bmin;
		realE = Emin;
		loopmin = VGmin;
		loopmax = VGmax;
		
		 //additional parameters for gate voltage type loop.
		 for(i=1; i<argc-1; i++)
		 {
		    if(strcmp("-vgtype", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &vgtype);
		    }
		    if(strcmp("-vgedgecut", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &edge_cut);
		    }
		    if(strcmp("-vgsubsthick", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &subs_thick);
		    }
		    if(strcmp("-vgsubsepsr", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &subs_epsr);
		    }
		 }
		if(vgtype==0)
		{
		  sprintf(vgtname, "constpot");
		}
		else if (vgtype ==1)
		{
		  sprintf(vgtname, "efetov");
		}
		else if (vgtype ==2)
		{
		  sprintf(vgtname, "efetov2");
		}
		
	
		
		sprintf(loopinfo, "VG_%s_loop_%+.2lf_to_%+.2lf_Bfixed_%+.3lf_Efixed_%+.3lf", vgtname, VGmin, VGmax, Bfield, realE);
	      }

	      
	      
	    //loop details
	    int calc_procs = procs;
	    int this_calc_proc = this_proc;
	    if(splitgen == 1)
	    {
	      calc_procs = procs-1;
	      this_calc_proc = this_proc-1;
	    }
	      
	      
		if(loop_pts>1)
		{
		  loop_step = (loopmax-loopmin)/(loop_pts - 1);
		  remainder = (loop_pts % calc_procs);
		}
		if(loop_pts==1)
		{
		   loop_step=0.0;
		   remainder=0;
		}
		
		
		if(remainder != 0 && this_calc_proc < remainder)
		{
		  loop_pts_temp = ((loop_pts - remainder) / calc_procs) + 1;
		  loop_min_temp = loopmin + this_calc_proc * (loop_pts_temp) * loop_step;
		  
		}
		else
		{
		  loop_pts_temp = loop_pts/calc_procs;
		  loop_min_temp = loopmin + remainder * (loop_pts_temp + 1) * loop_step + (this_calc_proc - remainder) * (loop_pts_temp) * loop_step ;
		}
		
		if(splitgen == 1 && this_proc == 0)
		{
		  loop_pts_temp=0;
		}
	

	//File I/O variables
	    FILE *output;
	    char filename[300], filename3[300], filename_temp[300];
	    char checkname[300], direcname[300], conffile[350], strucfile[300], disorderfile[300];
	    char bandname1[300], bandname2[300], bandname3[300], mapname[400];
            char command[600];
	    FILE *bandfile;
	    FILE *mapfile;
	    
	//Create directory and filenaming convention
	    sprintf(direcname, "../res/%s_%s_%.0e/%s%d_%s", systemtype, peritype, eta, geotype, length1, sysinfo);
	    
	    if(potdis==1)
	    {
	      sprintf(direcname, "%s/%s", direcname, disinfo);
	    }
	    
	    sprintf(command, "mkdir -p %s", direcname);
	    system(command);
	    printf("# directory: %s\n", direcname);
	    
	    
	    sprintf(filename_temp, "%s_%s.conf%02d", loopinfo, job_name, conf_num); 

	    sprintf(strucfile, "%s/%s.struct", direcname, filename_temp);
	    sprintf(disorderfile, "%s/%s.disprof", direcname, filename_temp);

	    sprintf(filename, "%s/.%s.part%02d", direcname, filename_temp, this_proc);
	    
	    sprintf(bandname1, "%s/%s.bands", direcname, filename_temp);
	    sprintf(bandname3, "%s/%s.wbands", direcname, filename_temp);
	    if(makebands ==1)
	    {
	      bandfile = fopen(bandname1, "w"); fclose(bandfile);
	      if(unfold==1)
	      {
		bandfile = fopen(bandname3, "w"); fclose(bandfile);
	      }
	    }
	    

	
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
			fprintf(check, "%d", -1);
			fclose(check);
			//output_type=1;
		}

		
	
		srand(this_proc); // ? why is this here??


    
      
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
	    
	    int can_start_yet=-1;
	    

	    //generate, or read in, disorder configuration
		sprintf(conffile, "%s/.%s", direcname, filename_temp);
		
		if(this_proc==0)
		{
		  (SysFunction) ( &System, SysPara, output_type, strucfile);
		  
		  //additional, general disorder routine(s) here.
		  if(potdis == 1)
		  {
		    potentialDisorder (&System, &dispara, pdmap, disorderfile );
		  }
		  
		  //export the disorder configuration for the other processes calculating the same configuration
		  if(procs>0)
		  {
		    exportRectConf(&System, conffile);
		  }
		  
		  check = fopen(checkname, "w");
		  fprintf(check, "%d", 0);
		  fclose(check);
		  
		}
		

// 		exit(0);
		

		  if(this_proc > 0)
		  {
		    while(can_start_yet<0)
		    {
			if(splitgen != 1)
			{
			  sleep( (int) myRandNum(1.0, 45.0));
			}
			check = fopen(checkname, "r");
			fscanf(check, "%d", &can_start_yet);
		    }
		    importRectConf(&System, length1, length2, conffile);
		  }
		  

		  time = clock() - time;
		  printf("#generated geometry in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		  time = clock();
		  
	    double **pos = System.pos;
	    int **siteinfo = System.siteinfo;
	    double *site_pots = System.site_pots;
	    
	     //in theory cell Division could write to lead_para, so its defined here
	  lead_para leadp={};
	  
	  

	  	  int halloutput;
		  
	    
		  
	  if(magsetup==0)
	  {
	    gauge = 0;
	  }
	  if(magsetup==1)
	  {
	    gauge = 1;
	  }
		  
		  
	  if(ishallbar == 1)
	  {
		if(strcmp("VG", loop_type) == 0)
		{
		  printf("hall bar and edge potential models incompatible!");
		  exit(1);
		}
	    
	    
		if(output_type == 1 && this_proc == 0)
		  halloutput = 1;
		
		else
		  halloutput = 0;
		
		
		HallBarify (&System, LeadCells, &hallp, &leadp, halloutput, strucfile);
		
			time = clock() - time;
		printf("#converted to Hall Bar in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();
		
		if(magsetup==0)
		{
		  gauge = 3;
		}
		if(magsetup==1)
		{
		  gauge = 2;
		}
		
		
		pos = System.pos;
		
		
		

	  }

		
	  
	  
	  
	  if(ishallbar != 1)
	  {
	    genLeads(&System, LeadCells, num_leads, 0, &leadp);
	  }
	  

  
	  
	  
	  

	  cnxProfile cnxp;
	  cnxp.max_neigh=max_neigh;
	  device_connectivity (&System, connectrules, &cnxpara, &cnxp);
	  
		time = clock() - time;
		printf("#made connection profile in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();
   	 // printConnectivity (&System, &cnxp);
	  
		
		
	  gen_start_params start_p ={};
	  start_p.rule = &graph_conn_sep2;
	  start_p.rule_params = &cnxpara;
	  start_p.num_leads = num_leads;
	  start_p.Leads = LeadCells;
	  
		
	  cellDivision cellinfo;
	 // genStartingCell(&System, &cellinfo, 2, NULL);
	  
	  	  genStartingCell(&System, &cellinfo, 3, &start_p);

	  
		time = clock() - time;
		printf("#generated starting cell in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();
	  
	  cellSplitter(&System, &cnxp, &cellinfo);
	  
		time = clock() - time;
		printf("#split cells in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();

// exit(0);
		
		
	  //hopping parameters and gauge info

//GAUGE NEEDS TO BE GENERALISED FOR HALL BAR CASE
//TIDY THIS UP IN GENERAL TO ALLOW FIELD IN LEADS
	  gen_hop_params hoppara={};
	  hoppara.num_neigh = nngm;
	  hoppara.hops=createCompArray(hoppara.num_neigh);
	  hoppara.NN_lowdis=createDoubleArray(hoppara.num_neigh);
	  hoppara.NN_highdis=createDoubleArray(hoppara.num_neigh);
	  hoppara.NN_shifts=createDoubleArray(hoppara.num_neigh);
	  
	  hoppara.hops[0]=t0;
	  hoppara.NN_lowdis[0] = NNlowdis;
	  hoppara.NN_highdis[0] = NNhighdis;
	  
	  if(nngm>1)
	  {
	    hoppara.hops[1]=t1;
	    hoppara.NN_lowdis[1]=NNlowdis1;
	    hoppara.NN_highdis[1]=NNhighdis1;
	    hoppara.NN_shifts[1]= onsitec1;
	  }
	    
	  if(nngm>2)
	  {
	    hoppara.hops[2]=t2;
	    hoppara.NN_lowdis[2]=NNlowdis2;
	    hoppara.NN_highdis[2]=NNhighdis2;
	    hoppara.NN_shifts[2]= onsitec2;
	    
	  }
// 	  hoppara.t0=t0;
	  hoppara.isperiodic=isperiodic;
	  hoppara.kpar=kmin;

	  hoppara.gauge=gauge;
	  int *res = createIntArray(2);
	  double **reslimits = createNonSquareDoubleMatrix(2, 6);
	  
	  
	  //mag field cut offs
	  if(gauge==0)
	  {
	      res[0] = 1;
	      reslimits[0][0] = pos[0][0];
	      reslimits[0][1] = pos[Ntot-1][0];
	  }
	      
	  
	  //this assumes bottom and top probes are collinear, I guess is can be generalised later
	  if(gauge == 2 || gauge == 3)
	  {
	    reslimits[0][0] = pos[2*length1*2][0];
	    reslimits[0][1] = pos[2*length1*(hallp.toppx[0] - (int)(buffer_rows/2) )][0];
	    reslimits[0][2] = pos[2*length1*(hallp.toppx[1] + (int)((hallp.toppw[1]+1)/2)  + (int)(buffer_rows/2))][0];
	    reslimits[0][3] = pos[2*length1*(length2 - 2)][0];
	    
	    reslimits[0][4] = pos[2*length1*(hallp.toppx[0] + (int)((hallp.toppw[0]+1)/2)  +1)][0];
	    reslimits[0][5] = pos[2*length1*(hallp.toppx[1] -1)][0];

	    
	    
	    
	    printf("# gauge changes points: %lf, %lf, %lf, %lf, %lf, %lf\n", reslimits[0][0], reslimits[0][1], reslimits[0][2], reslimits[0][3], reslimits[0][4], reslimits[0][5]);
	    
	  }
	  
	  if(gauge == 3)
	  {

	    
	    reslimits[1][0] = pos[2*length1*length2 + 2*(hallp.toppc[0]*hallp.toppw[0] +  hallp.toppc[1]*hallp.toppw[1] + hallp.botpc[0]*hallp.botpw[0])-1 ][1];
	    reslimits[1][3] = pos[2*length1*length2 + 2*hallp.toppc[0]*hallp.toppw[0]  -1][1];

	    

	    
	  }

	
	  hoppara.restrics = res;
	  hoppara.limits = reslimits;
	  
	  double _Complex **Sigma = createSquareMatrix(cellinfo.cell1dim);
	  double _Complex **G00 = createSquareMatrix(cellinfo.cell1dim);


	  
	  leadp.hopfn = hopfn;
	  leadp.hoppara = &hoppara;
// 	  
	 // leadp.leadsfn = &simple2leads;
	  	  leadp.leadsfn = &multipleLeads;

	  
	  trans_params tpara = {};
	  tpara.num_leads=num_leads;
	  tpara.TRsym=ismagnetic;
	  double **transmissions = createNonSquareDoubleMatrix(num_leads, num_leads);
	  double **ttildas=createNonSquareDoubleMatrix(num_leads, num_leads);
	  double **mtrans = createNonSquareDoubleMatrix(num_leads-1, num_leads-1);
	  double *vecx = createDoubleArray(num_leads-1);
	  double *vecb = createDoubleArray(num_leads-1);
	  

	  tpara.transmissions = transmissions;

	  double kavg;
	  
	  double kxl;
	  double *bands;
	  double **weights;
	  double **projections;
	  int bandmode;
	  double kxstep;
	  
	  double *ldoses ;
	  double ***currents;
	 
	  
	  
	  if(mappings != 0)
	  {
	    ldoses = createDoubleArray(Ntot);
	    currents = (double ***)malloc(num_leads * sizeof(double **));
	    for(i=0; i<num_leads; i++)
	    {
	      currents[i] = createNonSquareDoubleMatrix(Ntot, 2);
	    }
	  }
	  
	 double *engdeppots,  *origpots , **lead_origpots, **lead_engdeppots;
	 
	 //back up energy independent on-site potentials
	 if(strcmp("VG", loop_type) == 0)
	 {
	   engdeppots = createDoubleArray(Ntot);
	   origpots = createDoubleArray(Ntot);
	   for(i=0; i<Ntot; i++)
	   {
	     origpots[i] = System.site_pots[i];
	   }
	   
	   lead_engdeppots = (double **)malloc(num_leads*sizeof(double *));
	   lead_origpots = (double **)malloc(num_leads*sizeof(double *));
	   for(j=0; j<num_leads; j++)
	   {
	     lead_engdeppots[j] = createDoubleArray( *(LeadCells[j]->Ntot));
	     lead_origpots[j] = createDoubleArray( *(LeadCells[j]->Ntot));
	     
	     for(i=0; i<*(LeadCells[j]->Ntot); i++)
	     {
	      lead_origpots[j][i] = (LeadCells[j]->site_pots)[i];
	     }
	   }
	 }

	  
	 if(makebands == 1)
	 {
	   bands = createDoubleArray(Nrem);
	   weights = createNonSquareDoubleMatrix(Nrem, length2);
	   projections = createNonSquareDoubleMatrix(Nrem, Nrem);
	   bandmode=0;
	   kxstep = (2*M_PI/length2)/(kxpoints -1);
	   
	   
	    if(ismagnetic == 1)
	    {
	      hoppara.Btes=Bmin;
	    }
	    
	    //generate engdeppots if required
	    gate_induced_pot (vgtype, &System, engdeppots, VGmin, edge_cut, subs_thick, subs_epsr);

	    
	    
	    //set total onsite potentials
	    for(i=0; i<Ntot; i++)
	    {
	      System.site_pots[i] = origpots[i] + engdeppots[i];
	    }
	    
	    
	      for(kxl=0; kxl< 2*M_PI/length2; kxl += kxstep)
	      {
      //  	   kxl=0.2;
		 
		  
		  
		  if(project == 1)
		  {
		      
		      bandmode=1;
		  }
		  
		  if(unfold == 1)
		  {
		      bandmode =2;
		  }
			  
		  genKXbandproj(&System, hopfn, &hoppara, bandmode, kxl, bands, projections, weights);

		  
		  bandfile = fopen(bandname1, "a"); 
		  fprintf(bandfile, "%lf\t", kxl);
		  for (i=0; i<Nrem; i++)
		  {
		    fprintf(bandfile, "%e\t", bands[i]);
		  }
		  fprintf(bandfile, "\n");
		  fclose(bandfile);
		  
		  if(unfold ==1)
		  {
		    bandfile = fopen(bandname3, "a"); 
		    for(j=0; j<length2; j++)
		    {
		      for (i=0; i<Nrem; i++)
		      {
			fprintf(bandfile, "%lf \t %e	%e\n",  kxl + j*2*M_PI/length2, bands[i], weights[i][j]);
		      }
		    }
		    fclose(bandfile);
		  }
		  
		  if(project==1)
		  {
		      for (i=bandsminset; i<mymin(bandsmaxset, Nrem); i++)
		      {
			sprintf(bandname2, "%s/%s.proj_kx_%.2lf_%d", direcname, filename_temp, kxl, i);
			bandfile = fopen(bandname2, "w"); 
			
		      		      
			for (j=0; j<Nrem; j++)
			{
			  fprintf(bandfile, "%lf	%lf	%e\n", pos[j][0], pos[j][1], projections[i][j]);
			}
			
			fclose(bandfile);
		      }
		    
		  }
		  
		  
	      }
	      
	      
		if(bandsonly == 1)
		{
		  
		  //delete temporary files
		  if(conf_num != 0)
		    {
		      sprintf(command, "rm %s/.%s.*", direcname, filename_temp);
		      system(command);
		    }
		  exit(0);
		}
		
	    //reset onsite potentials if they have been gate-voltage model altered
	    
	      if(strcmp("VG", loop_type) == 0)
	      {
		for(i=0; i<Ntot; i++)
		{
		  System.site_pots[i] = origpots[i];
		}
	      }
	    
	    
	 }
	  
	  
	  
	  
	    
	  
	  
	  for(en=0; en<loop_pts_temp; en++)
	  {
	      if(strcmp("E", loop_type) == 0)
	      {
		realE =  loop_min_temp + en*loop_step;
	      }
	      
	      if(strcmp("B", loop_type) == 0)
	      {
		Bfield =  loop_min_temp + en*loop_step;

	      }
	      
	      if(strcmp("VG", loop_type) == 0)
	      {
		VG =  loop_min_temp + en*loop_step;
		
		//generate gate-depedent onsites for this VG
		  //in the device:
		    gate_induced_pot (vgtype, &System, engdeppots, VG, edge_cut, subs_thick, subs_epsr);

		  //in the leads
		    for(j=0; j<num_leads; j++)
		    {
		      	gate_induced_pot (vgtype, LeadCells[j], lead_engdeppots[j], VG, edge_cut, subs_thick, subs_epsr);
		    }

		
		//set total onsite pots
		    for(i=0; i<Ntot; i++)
		    {
		      System.site_pots[i] = origpots[i] + engdeppots[i];
		    }
		    
		    for(j=0; j<num_leads; j++)
		    {
		      for(i=0; i<*(LeadCells[j]->Ntot); i++)
		      {
			(LeadCells[j]->site_pots)[i]= lead_origpots[j][i] + lead_engdeppots[j][i];
		      }
		    }
				

	      }
	      
	      
	      if(mappings != 0)
	      {
		mapnow=0;
		
		if( (en % mappings) == 0)
		{
		  mapnow = 1;
		}
		  
	      }
	      mapmode = mapnow;
	      
// 	      printf("#map mode %d\n", mapmode);
	      
	      if(mapmode == 1)
	      {
		for (i=0; i<Ntot; i++)
		    { 
		      ldoses[i] = 0;
		      
		      for(j=0; j< num_leads; j++)
		      {
			for(k=0; k<2; k++)
			{
			  currents[j][i][k] = 0;
			}
		      }
		    }
	      }
	      
	     
	    
	      hoppara.Btes=Bfield;
	      
	      //needs generalising for multiple transmission cases
	      kavg=0.0;
	      
	      for(k=0; k<kpts; k++)
	      {
		hoppara.kpar = kmin + k*kstep;
		
		genTransmissions(realE+eta*I, &System, LeadCells, &cnxp, &cellinfo, hopfn, &hoppara, &leadp, &tpara, mapmode, ldoses, currents);
		kavg += (transmissions[0][1]/kpts);
	      }

	      output =fopen(filename, "a");
	      
	      if(ishallbar == 0)
	      {
		if(strcmp("E", loop_type) == 0)
		{
		  fprintf(output, "%lf	%e\n", realE, kavg);		
  // 		printf("%lf	%e	%e	%e	%e	%e\n", realE, transmissions[0][1], transmissions[0][2], transmissions[0][3], transmissions[0][4], transmissions[0][5]);

		  //printf( "%lf	%e\n", realE, kavg);

		}
		if(strcmp("B", loop_type) == 0)
		{
		  fprintf(output, "%lf	%e\n", Bfield, kavg);
  // 		printf("%lf	%e	%e	%e	%e	%e\n", Bfield, transmissions[0][1], transmissions[0][2], transmissions[0][3], transmissions[0][4], transmissions[0][5]);
		}
		
		if(strcmp("VG", loop_type) == 0)
		{
		  fprintf(output, "%lf	%e\n", VG, kavg);
		}
	      }
	      
	      
	      //print maps
	      double xc, yc;
	      double probe_pots[num_leads];
	      
	      if(mapmode==1)
	      {
		     sprintf(mapname, "%s/E_%+.2lf_B_%+.3lf_%s.conf%02d.ldos", direcname, realE, Bfield, job_name, conf_num);
		     

		      mapfile = fopen(mapname, "w");
		      
		      for(i=0; i<Ntot; i++)
		      {
			fprintf(mapfile, "%lf	%lf	%e\n", (pos)[i][0], (pos)[i][1], ldoses[i]);
		      }
		      fclose(mapfile);
		      

	      
		      for(j=0; j<num_leads; j++)
		      {
			
			sprintf(mapname, "%s/E_%+.2lf_B_%+.3lf_%s.conf%02d.cmaps_l%d", direcname, realE, Bfield, job_name, conf_num, j);


			  mapfile = fopen(mapname, "w");
			  for(i=0; i<Ntot; i++)
			  {
			    fprintf(mapfile, "%lf	%lf	%e %e\n", (pos)[i][0], (pos)[i][1], currents[j][i][0], currents[j][i][1] );
			  }
			  fclose(mapfile);
			  
		      }
		      
		      
		  
			   

		      
		    
	      }
	      
	      
	      
	      
	      //Standard 6 probe hall bar, outputs Rxy, Rxx and T_LR
		//consider alternative hall bar arrangements later?
		//generalise ishallbar=2, =3 etc for SHE setups?
	      if (ishallbar==1)
	      {
		EmptyDoubleMatrix(ttildas, num_leads, num_leads);
		EmptyDoubleMatrix(mtrans, num_leads-1, num_leads-1);
		
		for(i=0; i<num_leads; i++)
		{
		  for(j=0; j<num_leads; j++)
		  {
		    if(i!=j)
		    {
			    ttildas[i][i] += transmissions[j][i] ; 
			    ttildas[i][j] = -transmissions[i][j] ;
		    }
		  }
		}
		
		mtrans[0][0] = 1 / ttildas[0][0];
		
		for(i=1; i<num_leads-1; i++)
		{
		    mtrans[0][i] = - ttildas[0][i+1] / ttildas[0][0];
		    mtrans[i][0] =  ttildas[i+1][0] / ttildas[0][0];
		    
		    for(j=1; j< num_leads-1; j++)
		    {
		      mtrans[i][j] = ttildas[i+1][j+1] - (ttildas[i+1][0] * ttildas[0][j+1])/ttildas[0][0] ;
		      
		    }
		}
		vecb[0] = 1.0; vecb[1] = 0.0;  vecb[2] = 0.0; vecb[3] = 0.0; vecb[4] = 0.0; 
		vecx[0] = 0.0; vecx[1] = 0.0;  vecx[2] = 0.0; vecx[3] = 0.0; vecx[4] = 0.0; 
		
		LinEqnDouble (mtrans, vecb, vecx, num_leads-1);
		
		
		
		if(strcmp("E", loop_type) == 0)
		{
		  fprintf(output, "%lf	%.12e	%.12e	%.12e	%.12e	%.12e	%.12e	%.12e\n", realE, (vecx[1]-vecx[3])/vecx[0], (vecx[1]-vecx[2])/vecx[0], kavg, vecx[1], vecx[2], vecx[3], vecx[4]);		
  // 		printf("%lf	%e	%e	%e	%e	%e\n", realE, transmissions[0][1], transmissions[0][2], transmissions[0][3], transmissions[0][4], transmissions[0][5]);

		  //printf( "%lf	%e\n", realE, kavg);

		}
		if(strcmp("B", loop_type) == 0)
		{
		  fprintf(output, "%lf	%.12e	%.12e	%.12e	%.12e	%.12e	%.12e	%.12e\n", Bfield, (vecx[1]-vecx[3])/vecx[0], (vecx[1]-vecx[2])/vecx[0], kavg, vecx[1], vecx[2], vecx[3], vecx[4]);		
  // 		printf("%lf	%e	%e	%e	%e	%e\n", Bfield, transmissions[0][1], transmissions[0][2], transmissions[0][3], transmissions[0][4], transmissions[0][5]);
		}
		
		
		
		    //make a composite map for multiprobe systems
		      if(mapmode == 1)
		      {
			probe_pots[0] = 1.0;
			probe_pots[1] = 0.0;
			probe_pots[2] = vecx[1];
			probe_pots[3] = vecx[2];
			probe_pots[4] = vecx[3];
			probe_pots[5] = vecx[4];
			
			printf("#curr: %.10e\n", vecx[0]);

			sprintf(mapname, "%s/E_%+.2lf_B_%+.3lf_%s.conf%02d.cmaps_multi", direcname, realE, Bfield, job_name, conf_num);
			mapfile = fopen(mapname, "w");
			for(i=0; i<Ntot; i++)
			{ 
			  xc=0.0; yc=0.0;
 			  k=1;
//  			  for(k=0; k<num_leads; k++)
//  			  {
			    for(j=0; j<num_leads; j++)
			    {
				if(j!=k)
				{
				  if(probe_pots[k] > probe_pots[j])
				  {
				    xc += (probe_pots[k] - probe_pots[j])*currents[k][i][0];
				    yc += (probe_pots[k] - probe_pots[j])*currents[k][i][1];
				  }
				  
				  if(probe_pots[j] > probe_pots[k])
				  {
				    xc += (probe_pots[j] - probe_pots[k])*currents[j][i][0];
				    yc += (probe_pots[j] - probe_pots[k])*currents[j][i][1];
				  }
				}
			      
//  			    }
			  }
			  
			  
		
		
			  
			    fprintf(mapfile, "%lf	%lf	%.12e %.12e\n", (pos)[i][0], (pos)[i][1], xc, yc);
			}
			fclose(mapfile);
		      }
			  
		

	      }
	      
	      fclose(output);

	      
	      
	  }
			


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
// 	       printf("%
		if(conf_num>0)
		{
		  sleep((int) myRandNum(0.0, 10.0));
		}
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





