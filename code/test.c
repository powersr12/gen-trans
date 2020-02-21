
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
	      char geotype[32], peritype[40], leadtype[32], sysinfo[120], loopinfo[80], loopmininfo[80], disinfo[40], cedgeinfo[40], dedgeinfo[40], eedgeinfo[40];
	      sprintf(systemtype, "SUBLATTICEPOT");
	      int length1=2, length2=3*length1, geo=0, isperiodic=0, ismagnetic=0, ispatched=0;
	      int makebands = 0, unfold=0, project=0, kxpoints=51, bandsonly=0, bandsminset=0, bandsmaxset=100000;
              int numpatches = 0;
              char temp_in_string[40];

	      
	      int output_type=1;   //=0 no structure output, =1  atoms / holes
	      double eta = 1.0E-6;
	      int ishallbar=0;
                                //common complex lead geometries can be chosen here
                                //each have special settings throughout the code
                                //be sure that, e.g., gauges, are carefully configured for each....
                                //=0 (default) simple 2 lead left/right - works fine with periodicity
                                    //output is TLR
                                //=1 6 lead HALL BAR  - Rxy and Rxx are calculated and output based on left to right injected current
                                    //bond currents in 'multi' show net current flow
                                //=2 6 lead SPIN HALL BAR - same as Hall bar, but current is injected top to bottom, and the non local resistance (& spin hall angle if system is SPIN POLARISED) are output
                                    //bond currents in 'multi' show net current flow
                                    //shows total, and each spin version (will!)
				//= 3 / four metal leads in a Hanle type configuration.
				//=4 CUSTOM MODE
				//allows each lead to be defined separately.
              
              
	      int nngm = 1;	//default assumption is Nearest Neighbour graphene hopping parameters and lattice, 
				//changing this to 2 or 3 uses 2NN or 3NN models
				//future change - making this 0 should allow custom parameterisation somehow
	      
	      int num_neigh; 	//this is connected to nngm, but is the index used in hopping arrays.
				//used to allow the nngm=1 mode to have more than one hopping paramaters
				//e.g. in bilayers
	      
	      int magsetup = 0;	//this indicates whether the magnetic field is present in the leads or not
				//(if the system is magnetic)
				//0 means no field in leads
				//1 means field everywhere
				//the correct gauge is then chosen depending on whether we have
				//a simple system or a hall bar
				//2 means no field in right/left leads in hall bar geometry
				//this forces the choice gauge=0
	      int gauge = 0; 	//default gauge choice, phase is along x
	      int potdis =0;	//is there an additional, random potential disorder in the system?
	      int cedge = 0;	//is/are there (a) custom potential(s) near the ribbon edges
	      int dedge = 0;    //secondary edge potential
      	      int eedge = 0;    //tertiary edge potential (is this taking the piss a bit?)

	      
	      
	      int splitgen =0;	//if =1, splits the system generation and calculation
	      int metal_leads = 0;      //if this is set to 1, metallic leads are used instead of ribbons -
                                        //these need to be accounted for properly with positions etc...
                                        
	      
              int take_first_chain_pots = 0;          
              
              int abs_pots = 0;         // = 1 use absorbing potentials on top and bottom of nanoribbons
                                        // = 2 use absorbing potentials on all four sides of device
                                        // BE VERY CAREFUL WITH LEAD PLACEMENT!! (suggested to use METALXY, far from edges only)
                                        
              double abs_pots_width = 2.0; //width of absorbing region
                                        
                                        
	     //lead_para moved here to allow additional_params to be specified in system generation
	     //e.g. particularly for bilayers
	     
		lead_para leadp={};
		int gen_leads_mode=0;
	  
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
		    if(strcmp("-cedge", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &cedge);
		    }
		    if(strcmp("-dedge", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &dedge);
		    }
		    if(strcmp("-eedge", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &eedge);
		    }
		    if(strcmp("-takefirst", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &take_first_chain_pots);
		    }
		    if(strcmp("-abspots", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &abs_pots);
		    }
		    if(strcmp("-abspotsw", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%lf", &abs_pots_width);
		    }
		    if(strcmp("-patched", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &ispatched);
		    }
		    if(strcmp("-numpatches", argv[i]) == 0)
		    {
			sscanf(argv[i+1], "%d", &numpatches);
		    }
		    
		  }
			//default value of num_neigh for single layer graphene.
			//each mode increases the num of neighbour terms by 1
			num_neigh=nngm;
			 
                
                        
		
		  
	  
  
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
	gen_hop_params hoppara={};
	gen_hop_params *hop_to_load = &graphene_NNTB_hop_params;
	int max_neigh=(hop_to_load->max_neigh)[num_neigh-1];
	
		    
	
	
	//array related constants
	  int Ntot, Nrem;
	  
	
	//processor and energy/magnetic loop related integers
	  int procs=1, this_proc=0, remainder, en, numfin, conf_num=0;
	  int loop_pts=101, loop_pts_temp=0;
	  char loop_type[4];
	  sprintf(loop_type, "E"); //= 'E'; //E for energy loop, B for Bfield loop
	  
	  //how often to map, map this particular energy/bfield/loop variable, mapping mode: ldos, currents, all
	  int mappings=0, mapnow=0, mapmode=1;
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
		    
			  if(magsetup==2)
			  {
			    sprintf(job_name, "%s.ms2", job_name);
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
	    
	    
	    int starting_cell_mode=3;
	    cnxRulesFn *connectrules;
	    graph_conn_para cnxpara;
	    
	    cnxRulesFn *default_connect_rule =  &graph_conn_sep;
	    cnxRulesFn2 *defdouble_connect_rule =  &graph_conn_sep2;
	    void *default_connection_params = &cnxpara;
	   
	    //default graphene settings - can be changed for BLG, etc
	    
		cnxpara.conn_sep_thresh_min = (hop_to_load->NN_lowdis)[0];
		cnxpara.conn_sep_thresh_max = (hop_to_load->NN_highdis)[num_neigh-1];
		
	    

	    //lead info
	  int num_leads=2;
	  
		for(i=1; i<argc-1; i++)
		  {
			if(strcmp("-numleads", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%d", &num_leads);
			}
		  }
	  
	  
	      if(isperiodic==0)
	      {
// 		connectrules = default_connect_rule;
		cnxpara.periodic = 0;
		sprintf(peritype, "RIBBON");
		kpts = 1;

	      }
	      if(isperiodic==1)
	      {
// 		connectrules = default_connect_rule;
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
	  double interface_width=0.0, interface_position, interface_position2;
	  double intreltoedge=2.0;
	  int xleadsavgpot=0;
	  
	  
	  
	//disorder / system config parameters - for antidots overwrites some of the above...
	  int AD_length=5, AD_length2=5, lat_width=1, lat_length=2;
	  double AD_rad=1.0, xyfluc=0.0, radfluc=0.0;
	  char latgeo[24], dotgeo[24], bubgeo[24];
	  sprintf(latgeo, "trig");
	  sprintf(dotgeo, "circ");
          sprintf(bubgeo, "membrane");
	  
	
	  
	  
	  
	  
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
			  suba_conc = sublp.a_conc;
		      }
		      if(strcmp("-subapot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublp.a_pot));
			  suba_pot = sublp.a_pot;
		      }
		      if(strcmp("-subbconc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublp.b_conc));
			  subb_conc = sublp.b_conc;
		      }
		      if(strcmp("-subbpot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublp.b_pot));
			  subb_pot = sublp.b_pot;
		      }
		      if(strcmp("-bufferrows", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(sublp.buffer_rows));
		      }
		      if(strcmp("-xleadsavgpot", argv[i]) == 0)   //puts an average potential on the left/right leads
		      {
			  sscanf(argv[i+1], "%d", &xleadsavgpot);
		      }
		      
		    }
		    
			if(xleadsavgpot==1)
			{
				sprintf(job_name, "%s.xlavgp", job_name);
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
		    sprintf(sysinfo, "L2_%d_BUF_%d_SUBA_%.2lfx%.3lf_SUBB_%.2lfx%.3lf", length2, buffer_rows, (subintp.a_conc1), (subintp.a_pot1),(subintp.b_conc1), (subintp.b_pot1)); 
		}
		
		sub2int_para sub2intp = {};
		if(strcmp("SUBTWOINT", systemtype) == 0)
		{	
		    //default values
		    
		    sub2intp.buffer_rows = buffer_rows;
		    sub2intp.a_conc1 = suba_conc;
		    sub2intp.a_pot1 = suba_pot;
		    sub2intp.b_conc1 = subb_conc;
		    sub2intp.b_pot1 = subb_pot;
		    
		    sub2intp.xory = xory;
		    
		    
		    sub2intp.seed = conf_num;
		    
		

		    
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-subaconc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.a_conc1));
		      }
		      if(strcmp("-subapot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.a_pot1));
		      }
		      if(strcmp("-subbconc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.b_conc1));
		      }
		      if(strcmp("-subbpot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.b_pot1));
		      }
		       if(strcmp("-intxory", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(sub2intp.xory));
		      }
                        if(strcmp("-intwidth", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.int_width1));
		      }
		      if(strcmp("-intreltoedge", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &intreltoedge);
		      }
		      if(strcmp("-bufferrows", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(sub2intp.buffer_rows));
		      }
		      
		    }
		    
		      sub2intp.int_width2 =  sub2intp.int_width1;
			//default values of region 2 parameters are swaps of region 1
			sub2intp.a_conc2 = sub2intp.b_conc1;
			sub2intp.a_pot2 = sub2intp.b_pot1;
			sub2intp.b_pot2 = sub2intp.a_pot1;
			sub2intp.b_conc2 = sub2intp.a_conc1;
			
			//default values of region 3 parameters are same as region 1
			sub2intp.a_conc3 = sub2intp.a_conc1;
			sub2intp.a_pot3 = sub2intp.a_pot1;
			sub2intp.b_pot3 = sub2intp.b_pot1;
			sub2intp.b_conc3 = sub2intp.b_conc1;
			
			//default interface positions
			if(geo==0)
			{
			  if(sub2intp.xory==0)
			  {
			    interface_position = (int)(length2/4) *1.0 - 0.5;
			    interface_position2 = (int)(length2/4) *3.0 - 0.5;
			  }
			  
			  if(sub2intp.xory==1)
			  {
			    interface_position = intreltoedge;
			    interface_position2 = (length1*sqrt(3))/2 - intreltoedge;
			  }

			}
			
			if(geo==1)
			{
			  if(sub2intp.xory==0)
			  {
			    interface_position = (int)(length2/4) *sqrt(3) - 1/(2*sqrt(3));
			    interface_position2 = (int)(3*length2/4) *sqrt(3) - 1/(2*sqrt(3));
			  }
			  
			  if(sub2intp.xory==1)
			  {
			    interface_position = intreltoedge;
			    interface_position2 = length1/2 -0.5 - intreltoedge;
			  }
			}
			
			sub2intp.int_pos1 = interface_position;		
			sub2intp.int_pos2 = interface_position2;	
			
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-subaconc2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.a_conc2));
		      }
		      if(strcmp("-subapot2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.a_pot2));
		      }
		      if(strcmp("-subbconc2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.b_conc2));
		      }
		      if(strcmp("-subbpot2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.b_pot2));
		      }
		      
		        if(strcmp("-intpos", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.int_pos1));
		      }
		      
		          if(strcmp("-intpos2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.int_pos2));
		      }
		      
		          if(strcmp("-intwidth2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sub2intp.int_width2));
		      }
		    }
		    
		    
		    
		    printf("#interface 1: %lf, 2: %lf\n", (sub2intp.int_pos1), (sub2intp.int_pos2));
			

		    
		    //set functions and params for use below
		    SysFunction = &genSublatticeTwoInts;
		    SysPara = &sub2intp;
		    
		    //set filename info - what to put in filename from these params
		    sprintf(sysinfo, "L2_%d_%dINT_W%.1lf_SUBA_%.2lfx%.3lf_SUBB_%.2lfx%.3lf", length2, sub2intp.xory, sub2intp.int_width1, (sub2intp.a_conc2), (sub2intp.a_pot2),(sub2intp.b_conc2), (sub2intp.b_pot2)); 
		}
		
		
		
		submoire_para sublmp = {};
		double mev = 0.001/2.7;
		//mev=1.0;
		char moirename[40];
		if(strcmp("SUBLMOIRE", systemtype) == 0)
		{	
		    //default values
		    
			//mass and potential vales at high symmetry points
				sublmp.AA_mass = -40*mev;
				sublmp.AA_pot = 0*mev;
				sublmp.AB_mass = 30*mev;
				sublmp.AB_pot = 40*mev;
				sublmp.BA_mass = 15*mev;
				sublmp.BA_pot = -40*mev;
		    
			//Moire lattice vector zigzag direction length	
				sublmp.lM = 59;
				
			//Moire pattern origin offset
				sublmp.x_offset=0.0;
				sublmp.y_offset=0.0;
		    
				sprintf(moirename, "default");

		    
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-moireaamass", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublmp.AA_mass));
		      }
		      if(strcmp("-moireabmass", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublmp.AB_mass));
		      }
		      if(strcmp("-moirebamass", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublmp.BA_mass));
		      }
		       if(strcmp("-moireaapot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublmp.AA_pot));
		      }
		      if(strcmp("-moireabpot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublmp.AB_pot));
		      }
		      if(strcmp("-moirebapot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sublmp.BA_pot));
		      }
		      if(strcmp("-moirename", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s", moirename);
		      }
		      
		    }
		    
			
			
		    //set functions and params for use below
		    SysFunction = &genSublatticeMoire;
		    SysPara = &sublmp;
		    
		    //set filename info - what to put in filename from these params
		    sprintf(sysinfo, "L2_%d_%s_off_%.1lf_%.1lf", length2, moirename, sublmp.x_offset, sublmp.y_offset); 
		}
		
		
		adot_para adotp = {};
		char lattice_info[35], dot_info[45];
		if( strcmp("ANTIDOTS", systemtype) == 0 ||  strcmp("DOTS", systemtype) == 0)
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
                    adotp.orientation = 0;
		    
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
                          adotp.latgeo = latgeo;
		      }
		      if(strcmp("-dotgeo", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s", dotgeo);
                          adotp.dotgeo = dotgeo;
		      }
		      if(strcmp("-xyfluc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(adotp.xyfluc));
		      }
		      if(strcmp("-radfluc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(adotp.radfluc));
		      }
		      if(strcmp("-ADflip", argv[i]) == 0) //=0 "up", =1 "down", =2 random mix
		      {
			  sscanf(argv[i+1], "%d", &(adotp.orientation));
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
		    
		    if( strcmp("triAC", dotgeo) == 0 || strcmp("triZZ", dotgeo) == 0 )
		    {
                        if(adotp.orientation == 0)
                            sprintf(dot_info, "%s_dot_up_R_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (adotp.dotgeo), (adotp.AD_rad), (adotp.lat_width),  (adotp.lat_length), (adotp.xyfluc), (adotp.radfluc)); 
                        if(adotp.orientation == 1)
                            sprintf(dot_info, "%s_dot_down_R_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (adotp.dotgeo), (adotp.AD_rad), (adotp.lat_width),  (adotp.lat_length), (adotp.xyfluc), (adotp.radfluc)); 
                        if(adotp.orientation == 2)
                            sprintf(dot_info, "%s_dot_mix_R_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (adotp.dotgeo), (adotp.AD_rad), (adotp.lat_width),  (adotp.lat_length), (adotp.xyfluc), (adotp.radfluc)); 
		    }

		    
		    sprintf(sysinfo, "%s_%s", lattice_info, dot_info);
		    
		    //set functions and params for use below
                    
                    if( strcmp("ANTIDOTS", systemtype) == 0)
                        SysFunction = &genAntidotDevice;
                    
                    if( strcmp("DOTS", systemtype) == 0)
                        SysFunction = &genDotDevice;
                    
		    SysPara = &adotp;
		    
		}
		
		
		sldot_para sldotp = {};
		if(strcmp("SUBLDOTS", systemtype) == 0)
		{	
		    //default values -- use antidot defaults
		    
		    sldotp.buffer_rows = buffer_rows;
		    sldotp.SD_length = AD_length;
		    sldotp.SD_length2 = AD_length2;
		    sldotp.latgeo = latgeo;
		    sldotp.SD_rad = AD_rad;
		    sldotp.SD_rad2 = AD_rad;
		    sldotp.lat_width = lat_width;
		    sldotp.lat_length = lat_length;
		    sldotp.dotgeo = dotgeo;
		    sldotp.seed = conf_num;
		    sldotp.isperiodic = isperiodic;
		    sldotp.xyfluc = xyfluc;
		    sldotp.radfluc = radfluc;
                    adotp.orientation = 0;
		    
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-SDlength", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(sldotp.SD_length));
		      }
		      if(strcmp("-SDlength2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(sldotp.SD_length2));
		      }
		      if(strcmp("-SDrad", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sldotp.SD_rad));
			  sldotp.SD_rad2 = sldotp.SD_rad;
		      }
		       if(strcmp("-SDrad2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sldotp.SD_rad2));
		      }
		      if(strcmp("-latw", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(sldotp.lat_width));
		      }
		      if(strcmp("-latl", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(sldotp.lat_length));
		      }
		      if(strcmp("-bufferrows", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(sldotp.buffer_rows));
		      }
		      if(strcmp("-latgeo", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s",  latgeo);
                          sldotp.latgeo = latgeo;
		      }
		      if(strcmp("-dotgeo", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s", dotgeo);
                          sldotp.dotgeo = dotgeo;
		      }
		      if(strcmp("-xyfluc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sldotp.xyfluc));
		      }
		      if(strcmp("-radfluc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sldotp.radfluc));
		      }
		      
                      if(strcmp("-subaconc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sldotp.a_conc));
		      }
		      if(strcmp("-subapot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sldotp.a_pot));
		      }
		      if(strcmp("-subbconc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sldotp.b_conc));
		      }
		      if(strcmp("-subbpot", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(sldotp.b_pot));
		      }
		      if(strcmp("-SDflip", argv[i]) == 0) //=0 "up", =1 "down", =2 random mix
		      {
			  sscanf(argv[i+1], "%d", &(sldotp.orientation));
		      }
		      
		      
		      
		    }
		    
// 		    printf("#BUFFERS %d %d\n", buffer_rows, adotp.buffer_rows);
		    //antidot system sizes calculated from lattice details
		    if(strcmp("trig", latgeo) == 0)
		    {
		      if(geo==0)
		      {
			length1=2*(sldotp.lat_width)*(sldotp.SD_length) ; 
			length2= 3*(sldotp.SD_length)*(sldotp.lat_length) + 2*(sldotp.buffer_rows);
		      }
		      
		      if(geo==1)
		      {
			length1=6*(sldotp.lat_width)*(sldotp.SD_length) ; 
			length2= (sldotp.SD_length)*(sldotp.lat_length) + 2*(sldotp.buffer_rows);
		      }
		      
		      sprintf(lattice_info, "%s_lat_L_%d", (sldotp.latgeo),(sldotp.SD_length)); 

		    }
		    
		    if(strcmp("rotrig", latgeo) == 0)
		    {
		      if(geo==0)
		      {
			length1= 2*((sldotp.SD_length)+1)*(sldotp.lat_width); 
			length2=  ((sldotp.SD_length) + 1)*(sldotp.lat_length) + 2*(sldotp.buffer_rows);
		      }
		      
		      if(geo==1)
		      {
			length1=(2*(sldotp.SD_length) + 2)*(sldotp.lat_width) ; 
			length2= ((sldotp.SD_length)+1)*(sldotp.lat_length) + 2*(sldotp.buffer_rows);
		      }
		      
		      sprintf(lattice_info, "%s_lat_L_%d", (sldotp.latgeo),(sldotp.SD_length)); 
		    }
		    
		    if(strcmp("rect", latgeo) == 0)
		    {
		      if(geo==0)
		      {
			length1=2*(sldotp.lat_width)*(sldotp.SD_length) ; 
			length2= (sldotp.SD_length2)*(sldotp.lat_length) + 2*(sldotp.buffer_rows);
		      }
		      
		      if(geo==1)
		      {
			length1=2*(sldotp.lat_width)*(sldotp.SD_length2) ; 
			length2= (sldotp.SD_length)*(sldotp.lat_length) + 2*(sldotp.buffer_rows);
		      }
		      sprintf(lattice_info, "%s_lat_L_%d_%d", (sldotp.latgeo),(sldotp.SD_length), (sldotp.SD_length2)); 

		    }
		    
		    if(strcmp("circ", dotgeo) == 0 || strcmp("hexAC", dotgeo) == 0 || strcmp("hexZZ", dotgeo) == 0)
		    {
		      sprintf(dot_info, "%s_dot_R_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (sldotp.dotgeo), (sldotp.SD_rad), (sldotp.lat_width),  (sldotp.lat_length), (sldotp.xyfluc), (sldotp.radfluc)); 
		    }
		    
		    if(strcmp("rect", dotgeo) == 0 )
		    {
		      sprintf(dot_info, "%s_dot_R_%.1lf_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (sldotp.dotgeo), (sldotp.SD_rad), (sldotp.SD_rad2), (sldotp.lat_width),  (sldotp.lat_length), (sldotp.xyfluc), (sldotp.radfluc)); 
		    }
		    
		    if( strcmp("triAC", dotgeo) == 0 || strcmp("triZZ", dotgeo) == 0 )
		    {
                        if(sldotp.orientation == 0)
                            sprintf(dot_info, "%s_dot_up_R_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (sldotp.dotgeo), (sldotp.SD_rad), (sldotp.lat_width),  (sldotp.lat_length), (sldotp.xyfluc), (sldotp.radfluc)); 
                        if(sldotp.orientation == 1)
                            sprintf(dot_info, "%s_dot_down_R_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (sldotp.dotgeo), (sldotp.SD_rad), (sldotp.lat_width),  (sldotp.lat_length), (sldotp.xyfluc), (sldotp.radfluc)); 
                        if(sldotp.orientation == 2)
                            sprintf(dot_info, "%s_dot_mix_R_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (sldotp.dotgeo), (sldotp.SD_rad), (sldotp.lat_width),  (sldotp.lat_length), (sldotp.xyfluc), (sldotp.radfluc)); 
		    }


		    
		    sprintf(sysinfo, "%s_SUBA_%.2lfx%.3lf_SUBB_%.2lfx%.3lf_%s", lattice_info, (sldotp.a_conc), (sldotp.a_pot),(sldotp.b_conc), (sldotp.b_pot), dot_info);
		    
		    //set functions and params for use below
		    SysFunction = &genSublatticeDots;
		    SysPara = &sldotp;
		    
		}		
		
		
		
		bubble_para bubp = {};
                double bub_height = 1.0, zfluc=0.0;
                double ifluc=0.0;   //internal fluctuation in bubble perturbations 

 		char bub_info[66];   //(use lattice_info from antidots)
		if(strcmp("BUBBLES", systemtype) == 0)
		{	
		    //default values
		    
		    bubp.buffer_rows = buffer_rows;
		    bubp.cell_length = AD_length;
		    bubp.cell_length2 = AD_length2;
		    bubp.latgeo = latgeo;
                    bubp.bubgeo = bubgeo;
		    bubp.bub_rad = AD_rad;
		    bubp.bub_height = bub_height;
		    bubp.lat_width = lat_width;
		    bubp.lat_length = lat_length;
		    
		    bubp.seed = conf_num;
		    bubp.isperiodic = isperiodic;
		    bubp.xyfluc = xyfluc;
		    bubp.radfluc = radfluc;
                    bubp.zfluc = zfluc;
                    bubp.ifluc = ifluc;
		    
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-celllength", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(bubp.cell_length));
		      }
		      if(strcmp("-celllength2", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(bubp.cell_length2));
		      }
		      if(strcmp("-bubrad", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(bubp.bub_rad));
		      }
		       if(strcmp("-bubheight", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(bubp.bub_height));
		      }
		      if(strcmp("-latw", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(bubp.lat_width));
		      }
		      if(strcmp("-latl", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(bubp.lat_length));
		      }
		      if(strcmp("-bufferrows", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(bubp.buffer_rows));
		      }
		      if(strcmp("-latgeo", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s",  latgeo);
                          bubp.latgeo = latgeo;
		      }
		      if(strcmp("-bubgeo", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s", bubgeo);
                          bubp.bubgeo = bubgeo;
		      }
		      if(strcmp("-xyfluc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(bubp.xyfluc));
		      }
		      if(strcmp("-radfluc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(bubp.radfluc));
		      }
		      if(strcmp("-zfluc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(bubp.zfluc));
		      }
		      if(strcmp("-ifluc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(bubp.ifluc));
		      }
		      
		    }
		    
// 		    printf("#BUFFERS %d %d\n", buffer_rows, adotp.buffer_rows);
		    //antidot system sizes calculated from lattice details
		    if(strcmp("trig", latgeo) == 0)
		    {
		      if(geo==0)
		      {
			length1=2*(bubp.lat_width)*(bubp.cell_length) ; 
			length2= 3*(bubp.cell_length)*(bubp.lat_length) + 2*(bubp.buffer_rows);
		      }
		      
		      if(geo==1)
		      {
			length1=6*(bubp.lat_width)*(bubp.cell_length) ; 
			length2= (bubp.cell_length)*(bubp.lat_length) + 2*(bubp.buffer_rows);
		      }
		      
		      sprintf(lattice_info, "%s_lat_L_%d", (bubp.latgeo),(bubp.cell_length)); 

		    }
		    
		    if(strcmp("rotrig", latgeo) == 0)
		    {
		      if(geo==0)
		      {
			length1= 2*((bubp.cell_length)+1)*(bubp.lat_width); 
			length2=  ((bubp.cell_length) + 1)*(bubp.lat_length) + 2*(bubp.buffer_rows);
		      }
		      
		      if(geo==1)
		      {
			length1=(2*(bubp.cell_length) + 2)*(bubp.lat_width) ; 
			length2= ((bubp.cell_length)+1)*(bubp.lat_length) + 2*(bubp.buffer_rows);
		      }
		      
		      sprintf(lattice_info, "%s_lat_L_%d", (bubp.latgeo),(bubp.cell_length)); 
		    }
		    
		    if(strcmp("rect", latgeo) == 0)
		    {
		      if(geo==0)
		      {
			length1=2*(bubp.lat_width)*(bubp.cell_length) ; 
			length2= (bubp.cell_length2)*(bubp.lat_length) + 2*(bubp.buffer_rows);
		      }
		      
		      if(geo==1)
		      {
			length1=2*(bubp.lat_width)*(bubp.cell_length2) ; 
			length2= (bubp.cell_length)*(bubp.lat_length) + 2*(bubp.buffer_rows);
		      }
		      sprintf(lattice_info, "%s_lat_L_%d_%d", (bubp.latgeo),(bubp.cell_length), (bubp.cell_length2)); 

		    }
		    
		    //if(strcmp("circ", dotgeo) == 0 || strcmp("hexAC", dotgeo) == 0 || strcmp("hexZZ", dotgeo) == 0)
		    //{
		      sprintf(bub_info, "%s_bub_R_%.1lf_H_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf_zf_%.1lf_if_%.02lf", (bubp.bubgeo), (bubp.bub_rad), (bubp.bub_height), (bubp.lat_width),  (bubp.lat_length), (bubp.xyfluc), (bubp.radfluc), (bubp.zfluc), (bubp.ifluc)); 
		    //}
		    
// 		    if(strcmp("rect", dotgeo) == 0 )
// 		    {
// 		      sprintf(dot_info, "%s_dot_R_%.1lf_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (adotp.dotgeo), (adotp.AD_rad), (adotp.AD_rad2), (adotp.lat_width),  (adotp.lat_length), (adotp.xyfluc), (adotp.radfluc)); 
// 		    }

		    
		    sprintf(sysinfo, "%s_%s", lattice_info, bub_info);
		    
		    //set functions and params for use below
		    SysFunction = &genBubbleDevice;
		    SysPara = &bubp;
                    hopfn = &strainedTB;
		    
		}
		
		
		
		fold_para foldp = {};
                

 		char foldtype[66], farrange[66];   
        sprintf(foldtype, "gaussfold");
        sprintf(farrange, "equaly");
        double foldh=4.0, fw1=1000.0, fw2=4.0, fangle=0.0, fangf=0.0;
        double foldedge= 30.0, foldsteep= 0.2;
        int numfolds = 10;
                
		if(strcmp("FOLDS", systemtype) == 0)
		{	
		    //default values
		    
		    foldp.foldgeo = foldtype;
		    foldp.numfolds = numfolds;
			foldp.arrange = farrange;
			foldp.width1=fw1;
			foldp.width2=fw2;
			foldp.height=foldh;
			foldp.angle=fangle*M_PI;
			foldp.w1fluc=0.0; foldp.w2fluc=0.0;  foldp.hfluc=0.0; foldp.angfluc=0.0; foldp.posfluc=0.0;
			foldp.seed = conf_num;
			foldp.edgepos = foldedge;
			foldp.edgesteep = foldsteep;
                    
		    
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-foldgeo", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%s", (foldp.foldgeo));
		      }
		      if(strcmp("-farrange", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%s", (foldp.arrange));
		      }
		      if(strcmp("-numfolds", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%d", &(foldp.numfolds));
		      }
		      if(strcmp("-fw1", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%lf", &(foldp.width1));
		      }
		      if(strcmp("-fw2", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%lf", &(foldp.width2));
		      }
		      if(strcmp("-fheight", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%lf", &(foldp.height));
		      }
		      if(strcmp("-fangle", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%lf", &fangle);
				foldp.angle = M_PI * fangle;
		      }
		      if(strcmp("-fedge", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%lf", &(foldp.edgepos));
		      }
		      if(strcmp("-fsteep", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%lf", &(foldp.edgesteep));
		      }
		      
		      if(strcmp("-fw1fluc", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%lf", &(foldp.w1fluc));
		      }
		      if(strcmp("-fw2fluc", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%lf", &(foldp.w2fluc));
		      }
		      if(strcmp("-fhfluc", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%lf", &(foldp.hfluc));
		      }
		      if(strcmp("-fposfluc", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%lf", &(foldp.posfluc));
		      }
		      if(strcmp("-fangfluc", argv[i]) == 0)
		      {
				sscanf(argv[i+1], "%lf", &fangf);
				foldp.angfluc = M_PI * fangf;
		      }
	
		    }
		    
 		    
			if( (foldp.w1fluc) == 0.0 && (foldp.w2fluc) == 0.0 && (foldp.hfluc) == 0.0 && (foldp.angfluc) == 0.0 && (foldp.posfluc)==0)
			{
				sprintf(sysinfo, "L2_%d_N%d_%s_%s_%.1lf_x_%.1lf_x_%.1lf_%.1lfPI", length2, (foldp.numfolds), (foldp.foldgeo),(foldp.arrange), (foldp.width1), (foldp.width2), (foldp.height), (foldp.angle) );
			}
			else 
			{
				sprintf(sysinfo, "L2_%d_N%d_%s_%s_%.1lfpm%.1lf_x_%.1lfpm%.1lf_x%.1lfpm%.1lf_%.1lfpm%.1lfPI_pf%.1lf", 
				length2, (foldp.numfolds), (foldp.foldgeo),(foldp.arrange), (foldp.width1), (foldp.w1fluc), (foldp.width2), (foldp.w2fluc),
				 (foldp.height), (foldp.hfluc), (foldp.angle), (foldp.angfluc),  (foldp.posfluc) );
			}
       
		    
		    //set functions and params for use below
		    SysFunction = &genStrainFolds;
		    SysPara = &foldp;
            hopfn = &strainedTB;
		    
		}		
		
		
		
		
		
		
		symstrain_para symsp = {};
                double symsmag = 1.0, symswid=1.0;
                int leadstrain=0;

 		char symstype[66];   
                sprintf(symstype, "gaussfold");
                
		if(strcmp("SYMSTRAIN", systemtype) == 0)
		{	
		    //default values
		    
		    symsp.straingeo = symstype;
                    symsp.strain_mag = symsmag;
                    symsp.strain_width = symswid;
                    symsp.location = 0;
                    
		    
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-symsgeo", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s", (symsp.straingeo));
		      }
		      if(strcmp("-symsmag", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(symsp.strain_mag));
		      }
		      if(strcmp("-symswidth", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(symsp.strain_width));
		      }
		       if(strcmp("-symsloc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(symsp.location));
		      }
		       if(strcmp("-leadstrain", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &leadstrain);
		      }
		    }
		    
// 		    
                    if( (symsp.location) == 0)
                    {
                        sprintf(sysinfo, "L2_%d_%s_centre_%.1lf_%.1lf", length2, (symsp.straingeo),(symsp.strain_mag), (symsp.strain_width));
                    }
                    if( (symsp.location) == 1)
                    {
                        sprintf(sysinfo, "L2_%d_%s_edges_%.1lf_%.1lf", length2, (symsp.straingeo),(symsp.strain_mag), (symsp.strain_width));
                    }
                    
                    if(leadstrain==0)
                        sprintf(sysinfo, "%s_nl", sysinfo);
		    
		    //set functions and params for use below
		    SysFunction = &genSymStrain;
		    SysPara = &symsp;
                    hopfn = &strainedTB;
		    
		}
		
		
		
		randstrain_para randsp = {};
                double randsmag = 1.0, randswid=1.0;

		if(strcmp("RANDOMSTRAIN", systemtype) == 0)
		{	
		    //default values
		    
                    randsp.strain_mag = 0.05;
                    randsp.strain_width = 2;
                    randsp.location = 0;
                    randsp.buffer_rows = 1;
                    randsp.seed = conf_num;

		    
		    //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      
		      if(strcmp("-randsmag", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(randsp.strain_mag));
		      }
		      if(strcmp("-randswidth", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%lf", &(randsp.strain_width));
		      }
		       if(strcmp("-randsloc", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(randsp.location));
		      }
		      if(strcmp("-bufferrows", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%d", &(randsp.buffer_rows));
		      }
		    }
		    
// 		    
                    if( (randsp.location) == 0)
                    {
                        sprintf(sysinfo, "L2_%d_BUF_%d_S_%.3lf", length2, randsp.buffer_rows,(randsp.strain_mag));
                    }
                    if( (randsp.location) == 1)
                    {
                        sprintf(sysinfo, "L2_%d_BUF_%d_S_edges_%.3lf_%.1lf", length2, randsp.buffer_rows, (randsp.strain_mag), (randsp.strain_width));
                    }
                    
                  
		    
		    //set functions and params for use below
		    SysFunction = &genRandStrain;
		    SysPara = &randsp;
                    hopfn = &strainedTB;
		    
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
		
		simple558_para p558 = {};
		p558.cellswide=4;
		p558.GBpos=2;
		p558.anddis=0;
		p558.andD=0.0;
		p558.andW=0.0;
		p558.vacdis=0;
		p558.vacW=0.0;
		p558.vacP=0.0;
		p558.ycut=0.0;
		p558.seed = conf_num;
		if(strcmp("GB558", systemtype) == 0)
		{	
		    	
			
		for(i=1; i<argc-1; i++)
		{
			if(strcmp("-GBcells", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%d", &(p558.cellswide));
			}
			if(strcmp("-GBpos", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%d", &(p558.GBpos));
			}
			if(strcmp("-GBand", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%d", &(p558.anddis));
			}
			if(strcmp("-GBandD", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%lf", &(p558.andD));
			}
			if(strcmp("-GBandW", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%lf", &(p558.andW));
			}
			if(strcmp("-GBvacs", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%d", &(p558.vacdis));
			}
			if(strcmp("-GBvacW", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%lf", &(p558.vacW));
			}
			if(strcmp("-GBvacP", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%lf", &(p558.vacP));
			}
			if(strcmp("-GBycut", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%lf", &(p558.ycut));
			}
		}
			
			
			
		    //set functions and params for use below
		    SysFunction = &simple558GB;
		    SysPara = &p558;
		    
		    
		    
		    length1=4*(p558.cellswide); 
		    geo=1;
		    nngm=1;
		    
		    //set filename info?
			sprintf(sysinfo, "GB558_l2_%d", length2);
			
			if(p558.anddis==1)
			{
				sprintf(sysinfo, "GB558_l2_%d_ANDL_%.1lf_ANDW_%.2lf_yc_%.0lf", length2, p558.andD, p558.andW, p558.ycut);
			}
			if(p558.vacdis==1)
			{
				sprintf(sysinfo, "GB558_l2_%d_VACW_%.1lf_VACP_%.4lf_yc_%.0lf", length2, p558.vacW, p558.vacP, p558.ycut);
			}
			if(p558.anddis==1 && p558.vacdis==1)
			{
				sprintf(sysinfo, "GB558_l2_%d_ANDL_%.1lf_ANDW_%.2lf_VACW_%.1lf_VACP_%.4lf_yc_%.0lf", length2, p558.andD, p558.andW, p558.vacW, p558.vacP, p558.ycut);
			}
			

		}
		
		
		
		if(strcmp("CLEAN", systemtype) == 0)
		{	
		    
		    //set functions and params for use below
		    SysFunction = &simpleRibbonGeo;
		    SysPara = NULL;
		    
		    //set filename info?
		    	sprintf(sysinfo, "clean_l2_%d", length2);

		}
		
	  
		//BILAYER SYSTEMS
		bilayer_para bilp = {};
		blg_conn_para blg_cnxpara;
		double *subpots = createDoubleArray(4);
		double blgbias=0.0;
                int blg_ind_pots = 0;  //whether individual lattice potentials are set
		
		if(strcmp("BLG", systemtype) == 0)
		{	
			//this system uses the bilayer graphene hopping parameters
			//these are defined in 'useful_hops.c'
			hop_to_load = &BLG_NNTB_hop_params;
			
			
			
			//default values
			nngm=1; //mode 1 includes NN intra layer and direct inter layer.
				//mode 2 adds skew hopping terms to and from nondimer sites
				//mode 0 (never called), would decouple the two layers completely
				//(would it??)
				
			
			bilp.type_shift=1;	//default is AB
			bilp.shift_angle=0;
			char blgtype[6];
			double blg_shift[2];
			bilp.shift_vec = blg_shift;
			bilp.zsep = 1.362;
			
			//default is to not include skew hopping terms between layers
			bilp.skews = 0;
			bilp.subpots = subpots;
			
			
			
			//check for command line arguments which vary these
			for(i=1; i<argc-1; i++)
			{
				if(strcmp("-blgtype", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &(bilp.type_shift));
					//0=AA, 1=AB, 2=CUSTOM
				}
				if(strcmp("-blgshiftx", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &(bilp.shift_vec)[0]);
				}
				if(strcmp("-blgshifty", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &(bilp.shift_vec)[1]);
				}
				//note shift angle not yet implemented
				if(strcmp("-blgshiftangle", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &(bilp.shift_angle));
				}
				if(strcmp("-blgskews", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &(bilp.skews));
				}
				if(strcmp("-blgbias", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &blgbias);
					subpots[0] = -blgbias;
					subpots[1] = -blgbias;
					subpots[2] = blgbias;
					subpots[3] = blgbias;
				}
				if(strcmp("-blgA1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subpots[0]);
					blg_ind_pots = 1;
				}
                                if(strcmp("-blgB1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subpots[1]);
					blg_ind_pots = 1;
				}
                                if(strcmp("-blgA2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subpots[2]);
					blg_ind_pots = 1;
				}
                                if(strcmp("-blgB2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subpots[3]);
					blg_ind_pots = 1;
				}
			}
			
			
			if(bilp.type_shift ==1)
			{
				sprintf(blgtype, "AB");
			}
			else if(bilp.type_shift ==12)
			{
				sprintf(blgtype, "ABal");
			}
			else if (bilp.type_shift ==0)
			{
				sprintf(blgtype, "AA");
			}
			else if (bilp.type_shift ==2)
			{
				sprintf(blgtype, "CUST");
			}

			//"blgskews" from command line turns on skew hopping terms if 1NN is being used
			//else nngm in command line takes precedence.
			if(bilp.skews == 1)
			{
				nngm = 2;
			}
			num_neigh=nngm+1;
			max_neigh=(hop_to_load->max_neigh)[num_neigh-1];
			
			
			//basic BLG connection stuff
				blg_cnxpara.intra_thresh_min = BLG_NNTB_hop_params.NN_lowdis[0];
				blg_cnxpara.intra_thresh_max = BLG_NNTB_hop_params.NN_highdis[0];
				blg_cnxpara.inter_thresh_min = 0.000;
				blg_cnxpara.inter_thresh_max = BLG_NNTB_hop_params.NN_highdis[num_neigh-1];
				blg_cnxpara.zthresh1 = 0.1;
				blg_cnxpara.zthresh2 = bilp.zsep + 0.1;
		
			
			
			blg_cnxpara.periodic = cnxpara.periodic;
			

			
			
			//set functions and params for use below
			SysFunction = &simpleBilayerGeo;
			SysPara = &bilp;
			leadp.additional_params = &bilp;
			gen_leads_mode=1;
			
			default_connect_rule = &blg_conn_sep;
			defdouble_connect_rule = &blg_conn_sep2;
			default_connection_params = &blg_cnxpara;
			
			
			
			//set filename info - what to put in filename from these params
			sprintf(sysinfo, "%s_L2_%d", blgtype, length2);
			if(blgbias !=0.0)
			{
				sprintf(sysinfo, "%s_b%.2lf_L2_%d", blgtype, blgbias, length2);
			}
			
			if(blg_ind_pots != 0)
			{
				sprintf(sysinfo, "%s_A1_%.3lf_B1_%.3lf_A1_%.3lf_B1_%.3lf_L2_%d", blgtype, subpots[0], subpots[1], subpots[2], subpots[3], length2);
			}
			
			if(bilp.skews == 1)
			{
				sprintf(sysinfo, "%s.s_L2_%d", blgtype, length2); 
				if(blgbias !=0.0)
				{
					sprintf(sysinfo, "%s.s_b%.2lf_L2_%d", blgtype, blgbias, length2);
				}
			}
			
			
			
		}

	  
	  
                //Create a bilayer device with random distributions of (sublattice dependent) potentials
                blgpots_para blgp = {};
                int use_av_pots=0;  //use average potentials in the leads / buffer regions, otherwise zero
                double *subpotsc = createDoubleArray(4);
                double *subconcs = createDoubleArray(4);
                for(i=0;i<4;i++)
                {
                    subconcs[i]=1.0;
                    subpotsc[i]=0.0;
                }

                if(strcmp("BLGPOTS", systemtype) == 0)
		{	
			//this system uses the bilayer graphene hopping parameters
			//these are defined in 'useful_hops.c'
			hop_to_load = &BLG_NNTB_hop_params;
			
			//default values
			nngm=1; //mode 1 includes NN intra layer and direct inter layer.
				//mode 2 adds skew hopping terms to and from nondimer sites
				//mode 0 (never called), would decouple the two layers completely
				//(would it??)
				
			bilp.type_shift=1;	//default is AB
			bilp.shift_angle=0;
			char blgtype[6];
			double blg_shift[2];
			bilp.shift_vec = blg_shift;
			bilp.zsep = 1.362;
			
			//default is to not include skew hopping terms between layers
			bilp.skews = 0;
			bilp.subpots = subpotsc;
                        blgp.subpots = subpots;
                        blgp.subconcs = subconcs;
                        blgp.seed = conf_num;

                        
			
			//check for command line arguments which vary these
			for(i=1; i<argc-1; i++)
			{
				if(strcmp("-blgtype", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &(bilp.type_shift));
					//0=AA, 1=AB, 2=CUSTOM, 12 for alpha-aligned AGNRs
				}
				if(strcmp("-blgshiftx", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &(bilp.shift_vec)[0]);
				}
				if(strcmp("-blgshifty", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &(bilp.shift_vec)[1]);
				}
				//note shift angle not yet implemented
				if(strcmp("-blgshiftangle", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &(bilp.shift_angle));
				}
				if(strcmp("-blgskews", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &(bilp.skews));
				}
				if(strcmp("-blgA1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subpots[0]);
				}
                                if(strcmp("-blgB1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subpots[1]);
				}
                                if(strcmp("-blgA2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subpots[2]);
				}
                                if(strcmp("-blgB2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subpots[3]);
				}
				if(strcmp("-blgA1c", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subconcs[0]);
				}
                                if(strcmp("-blgB1c", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subconcs[1]);
				}
                                if(strcmp("-blgA2c", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subconcs[2]);
				}
                                if(strcmp("-blgB2c", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &subconcs[3]);
				}
				if(strcmp("-useavpots", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &use_av_pots);
				}
			}
			
			if(use_av_pots == 1)
                        {
                            for(i=0; i<4;i++)
                            {
                                subpotsc[i] = subpots[i] *subconcs[i];
                            }
                        }
			
			
			if(bilp.type_shift ==1)
			{
				sprintf(blgtype, "AB");
			}
			else if(bilp.type_shift ==12)
			{
				sprintf(blgtype, "ABal");
			}
			else if (bilp.type_shift ==0)
			{
				sprintf(blgtype, "AA");
			}
			else if (bilp.type_shift ==2)
			{
				sprintf(blgtype, "CUST");
			}

			//"blgskews" from command line turns on skew hopping terms if 1NN is being used
			//else nngm in command line takes precedence.
			if(bilp.skews == 1)
			{
				nngm = 2;
			}
			num_neigh=nngm+1;
			max_neigh=(hop_to_load->max_neigh)[num_neigh-1];
			
			
			//basic BLG connection stuff
				blg_cnxpara.intra_thresh_min = BLG_NNTB_hop_params.NN_lowdis[0];
				blg_cnxpara.intra_thresh_max = BLG_NNTB_hop_params.NN_highdis[0];
				blg_cnxpara.inter_thresh_min = 0.000;
				blg_cnxpara.inter_thresh_max = BLG_NNTB_hop_params.NN_highdis[num_neigh-1];
				blg_cnxpara.zthresh1 = 0.1;
				blg_cnxpara.zthresh2 = bilp.zsep + 0.1;
		
			blg_cnxpara.periodic = cnxpara.periodic;
			

                        blgp.BLGpara = &bilp;
						
			//set functions and params for use below
			SysFunction = &BLGPotentials;
			SysPara = &blgp;
			leadp.additional_params = &bilp;
			gen_leads_mode=1;
			
			default_connect_rule = &blg_conn_sep;
			defdouble_connect_rule = &blg_conn_sep2;
			default_connection_params = &blg_cnxpara;
			
			
			
// 			set filename info - what to put in filename from these params
                         sprintf(sysinfo, "%s_A1_%.3lf_%.3lf_B1_%.3lf_%.3lf_A2_%.3lf_%.3lf_B2_%.3lf_%.3lf_L2_%d", blgtype, subpots[0], subconcs[0], subpots[1], subconcs[1], subpots[2], subconcs[2], subpots[3], subconcs[3], length2);
                        
			if(use_av_pots ==1)
                            sprintf(sysinfo, "%s_avlp_A1_%.3lf_%.3lf_B1_%.3lf_%.3lf_A2_%.3lf_%.3lf_B2_%.3lf_%.3lf_L2_%d", blgtype, subpots[0], subconcs[0], subpots[1], subconcs[1], subpots[2], subconcs[2], subpots[3], subconcs[3], length2);
			
			
		}
		
		blgints_para blgintp = {};
                //int use_av_pots=0;  //use average potentials in the leads / buffer regions, otherwise zero
               // double **bisubpots, **bisubconcs ;

                
                double tempintwidth=0;
		char blgintconfig[100];
		//bilayer graphene sublattice (doping) interfaces
		if(strcmp("BLGINTS", systemtype) == 0)
		{	
			//this system uses the bilayer graphene hopping parameters
			//these are defined in 'useful_hops.c'
			hop_to_load = &BLG_NNTB_hop_params;
			
			//default values
			nngm=1; //mode 1 includes NN intra layer and direct inter layer.
				//mode 2 adds skew hopping terms to and from nondimer sites
				//mode 0 (never called), would decouple the two layers completely
				//(would it??)
				
			bilp.type_shift=1;	//default is AB
			bilp.shift_angle=0;
			char blgtype[6];
			double blg_shift[2];
			bilp.shift_vec = blg_shift;
			bilp.zsep = 1.362;
			
			//default is to not include skew hopping terms between layers
			bilp.skews = 0;
			bilp.subpots = subpotsc;
                        
                        blgintp.num = 1;
                        blgintp.seed = conf_num;
                        blgintp.xory = xory;
                        
                        
			//check for command line arguments which vary these
			for(i=1; i<argc-1; i++)
			{
				if(strcmp("-blgtype", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &(bilp.type_shift));
					//0=AA, 1=AB, 2=CUSTOM, 12 for alpha-aligned AGNRs
				}
				if(strcmp("-blgshiftx", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &(bilp.shift_vec)[0]);
				}
				if(strcmp("-blgshifty", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &(bilp.shift_vec)[1]);
				}
				//note shift angle not yet implemented
				if(strcmp("-blgshiftangle", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &(bilp.shift_angle));
				}
				if(strcmp("-blgskews", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &(bilp.skews));
				}
				if(strcmp("-useavpots", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &use_av_pots);
				}
				if(strcmp("-intnum", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &(blgintp.num));
				}
                                if(strcmp("-intxory", argv[i]) == 0)
                                {
                                    sscanf(argv[i+1], "%d", &(blgintp.xory));
                                }
                                if(strcmp("-intwidth", argv[i]) == 0)
                                {
                                    sscanf(argv[i+1], "%lf", &tempintwidth);
                                }
                                //default -- swap 
                                if(strcmp("-blgA1", argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &subpots[0]);
                                }
                                if(strcmp("-blgB1", argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &subpots[1]);
                                }
                                if(strcmp("-blgA2", argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &subpots[2]);
                                }
                                if(strcmp("-blgB2", argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &subpots[3]);
                                }
                                if(strcmp("-blgA1c", argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &subconcs[0]);
                                }
                                if(strcmp("-blgB1c", argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &subconcs[1]);
                                }
                                if(strcmp("-blgA2c", argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &subconcs[2]);
                                }
                                if(strcmp("-blgB2c", argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &subconcs[3]);
                                }
				 if(strcmp("-blintname", argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%s", blgintconfig);
                                }
                        }
                        
                        blgintp.subpots = createNonSquareDoubleMatrix(blgintp.num + 1, 4);
                        blgintp.subconcs = createNonSquareDoubleMatrix(blgintp.num + 1, 4);
                        blgintp.int_pos=createDoubleArray(blgintp.num+1);   //last not used
                        blgintp.int_width=createDoubleArray(blgintp.num+1); //last not used
                        
                        //default values and positions of interfaces and regions
                        for(j=0; j<blgintp.num +1; j++)
                        {
                            for(k=0;k<4;k++)
                            {
                                blgintp.subconcs[j][k] = subconcs[k];
                                
                                if(j%2 ==0)
                                    blgintp.subpots[j][k] = subpots[k];
                                if(j%2 ==1)
                                    blgintp.subpots[j][k] = -subpots[k];
                            }
                            blgintp.int_width[j] = tempintwidth;
                            
                            //default interface positions -- check these, not necessarily perfect 
                            if(geo==0)
                            {
                                if(blgintp.xory==0)
                                {
                                    blgintp.int_pos[j] = (int)(length2/(blgintp.num +1)) *(j+1)*1.0 - 0.5;
                                }
                                
                                if(blgintp.xory==1)
                                {
                                    blgintp.int_pos[j] = (int)(length1/(blgintp.num +1)) *(j+1)*(sqrt(3)/2);
                                }

                            }
                            
                            if(geo==1)
                            {
                                if(blgintp.xory==0)
                                {
                                    blgintp.int_pos[j] = (int)(length2/(blgintp.num +1)) *(j+1)*sqrt(3) - 1/(2*sqrt(3));
                                }
                                
                                if(blgintp.xory==1)
                                {
                                    blgintp.int_pos[j] = ((length1*1.0-0.5)/(blgintp.num +1)) *(j+1)*0.5 ;
                                }
                            }
                            
                            
                            
                        }
                        
                        
                        for(i=1; i<argc-1; i++)
			{
                            for(j=0; j<blgintp.num+1; j++)
                            {
                                
                                sprintf(temp_in_string, "-blgA1_%d", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &blgintp.subpots[j][0]   );
                                }
                                sprintf(temp_in_string, "-blgB1_%d", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &blgintp.subpots[j][1]   );
                                }
                                sprintf(temp_in_string, "-blgA2_%d", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &blgintp.subpots[j][2]   );
                                }
                                sprintf(temp_in_string, "-blgB2_%d", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &blgintp.subpots[j][3]   );
                                }
                                
                                sprintf(temp_in_string, "-blgA1c_%d", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &blgintp.subconcs[j][0]   );
                                }
                                sprintf(temp_in_string, "-blgB1c_%d", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &blgintp.subconcs[j][1]   );
                                }
                                sprintf(temp_in_string, "-blgA2c_%d", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &blgintp.subconcs[j][2]   );
                                }
                                sprintf(temp_in_string, "-blgB2c_%d", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &blgintp.subconcs[j][3]   );
                                }   
                                
                                sprintf(temp_in_string, "-intpos%d", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &blgintp.int_pos[j]   );
                                }   
                                
                                sprintf(temp_in_string, "-intwidth%d", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &blgintp.int_width[j]   );
                                }  
                                
                            }
				
			}
			
			if(use_av_pots == 1)
                        {
                            for(i=0; i<4;i++)
                            {
                                subpotsc[i] = subpots[i] *subconcs[i];
                            }
                        }
			
			
			if(bilp.type_shift ==1)
			{
				sprintf(blgtype, "AB");
			}
			else if(bilp.type_shift ==12)
			{
				sprintf(blgtype, "ABal");
			}
			else if (bilp.type_shift ==0)
			{
				sprintf(blgtype, "AA");
			}
			else if (bilp.type_shift ==2)
			{
				sprintf(blgtype, "CUST");
			}

			//"blgskews" from command line turns on skew hopping terms if 1NN is being used
			//else nngm in command line takes precedence.
			if(bilp.skews == 1)
			{
				nngm = 2;
			}
			num_neigh=nngm+1;
			max_neigh=(hop_to_load->max_neigh)[num_neigh-1];
			
			
			//basic BLG connection stuff
				blg_cnxpara.intra_thresh_min = BLG_NNTB_hop_params.NN_lowdis[0];
				blg_cnxpara.intra_thresh_max = BLG_NNTB_hop_params.NN_highdis[0];
				blg_cnxpara.inter_thresh_min = 0.000;
				blg_cnxpara.inter_thresh_max = BLG_NNTB_hop_params.NN_highdis[num_neigh-1];
				blg_cnxpara.zthresh1 = 0.1;
				blg_cnxpara.zthresh2 = bilp.zsep + 0.1;
		
			blg_cnxpara.periodic = cnxpara.periodic;
			

                        blgintp.BLGpara = &bilp;
						
			//set functions and params for use below
			SysFunction = &BLGInterfaces;
			SysPara = &blgintp;
			leadp.additional_params = &bilp;
			gen_leads_mode=1;
			
			default_connect_rule = &blg_conn_sep;
			defdouble_connect_rule = &blg_conn_sep2;
			default_connection_params = &blg_cnxpara;
			
 			sprintf(blgintconfig, "A1_%.3lf_%.3lf_B1_%.3lf_%.3lf_A2_%.3lf_%.3lf_B2_%.3lf_%.3lf", subpots[0], subconcs[0], subpots[1], subconcs[1], subpots[2], subconcs[2], subpots[3], subconcs[3]);
                        
                        
                        
                        
                        for(i=1; i<argc-1; i++)
                        {
                            if(strcmp("-blintname", argv[i]) == 0)
                            {
                                sscanf(argv[i+1], "%s", blgintconfig);
                            }
                        }
                        
                        
			
			
			//set filename info - what to put in filename from these params
                        sprintf(sysinfo, "%s_%s_L2_%d", blgtype, blgintconfig, length2);
			
			if(use_av_pots ==1)
                            sprintf(sysinfo, "%s_avlp_%s_L2_%d", blgtype, blgintconfig, length2);
			
			
		}
		
		
		
		//test device for generic multilayers
		multilayer_para mlp = {};
 		char ml_info[66];   //(filename gunk)
                char **layertype;
 		
		if(strcmp("MLTEST", systemtype) == 0)
		{	
                        
                        default_connect_rule = &mlg_conn_sep;
			defdouble_connect_rule = &mlg_conn_sep2;
			default_connection_params = &blg_cnxpara;
                        nngm=1;
                        num_neigh=nngm+1;
                    
                        //default values
                        
                        mlp.num_layers = 2;
                        mlp.hoptype = 0;	//by default use simplest (BLG-like) hoppings
											//hoptype = 1 uses a more complex, Moon-Koshino-esque approach for twisted systems
                        mlp.length = createIntArray(mlp.num_layers);
                        mlp.length2 = createIntArray(mlp.num_layers);
                        mlp.geo = createIntArray(mlp.num_layers);
                        mlp.theta = createDoubleArray(mlp.num_layers);
                        mlp.delta = createNonSquareDoubleMatrix(mlp.num_layers, 3);
                        mlp.epsilon = createDoubleArray(mlp.num_layers);
                        mlp.origin = 0; //0 means (0,0), 1 means centre of l1 x l2 device
                        mlp.layerfn = (generatefn **)malloc(mlp.num_layers * sizeof(generatefn *));
                        mlp.layerpara = (void **)malloc(mlp.num_layers * sizeof(void *));
                        layertype=malloc(mlp.num_layers * sizeof(char *));
                        for(i=0; i< mlp.num_layers ; i++)
                            layertype[i] = (char *)malloc(40  * sizeof(char)); 
					/*allocates memory for whole matrix, 
					returned pointe r= pointer to 1st array*/
					
					
							
                        

                        for(i=0; i< mlp.num_layers; i++)
                        {
                            (mlp.length)[i] = length1;
                            (mlp.length2)[i] = length2;
                            (mlp.geo)[i] = geo;
                            (mlp.theta)[i] = 0;
                            (mlp.epsilon)[i] = 1.0;

                            (mlp.delta)[i][0] = 0.0;
                            (mlp.delta)[i][1] = 0.0 + (1/sqrt(3)) * (i %2);
                            (mlp.delta)[i][2] = i * 1.362;
                            (mlp.layerfn)[i] = &simpleRibbonGeo;
                            (mlp.layerpara)[i]  = NULL;
                            sprintf(layertype[i], "CLEAN");
                            
                        }
                        
                 
                        
                        //check for command line arguments which vary these
                        for(i=1; i<argc-1; i++)
                        {
                            if(strcmp("-numlayers", argv[i]) == 0)
                            {
                                    sscanf(argv[i+1], "%d", &(mlp.num_layers));
                            }
                        }
                            
                        for(i=1; i<argc-1; i++)
                        {   
                            for(j=0; j<mlp.num_layers; j++)
                            {
                                sprintf(temp_in_string, "-lay%dl1", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%d", &(mlp.length)[j]);
                                }
                                sprintf(temp_in_string, "-lay%dl2", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%d", &(mlp.length2)[j]);
                                }
                                sprintf(temp_in_string, "-lay%dgeo", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%d", &(mlp.geo)[j]);
                                }
                                sprintf(temp_in_string, "-lay%ddeltax", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &(mlp.delta)[j][0]);
                                }
                                sprintf(temp_in_string, "-lay%ddeltay", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &(mlp.delta)[j][1]);
                                }
                                sprintf(temp_in_string, "-lay%ddeltaz", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &(mlp.delta)[j][2]);
                                }
                                sprintf(temp_in_string, "-lay%dtheta", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &(mlp.theta)[j]);
                                }
                                sprintf(temp_in_string, "-lay%deps", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &(mlp.epsilon)[j]);
                                }
                                
                                
                            }
                            //center of rotation and common origin of layers
                            if(strcmp("-mlorigin", argv[i]) == 0)
                            {
                                    sscanf(argv[i+1], "%d", &(mlp.origin));
                            }
                            // simple or moon-koshino hopping
                            if(strcmp("-mlhoptype", argv[i]) == 0)
                            {
                                    sscanf(argv[i+1], "%d", &(mlp.hoptype));
                            }
                        
                        }
                        
                        if (mlp.hoptype == 1)
							hop_to_load = &MLG_NNTB_hop_params;
								
						if (mlp.hoptype == 0)
							hop_to_load = &MLG2_NNTB_hop_params;
							
						max_neigh=(hop_to_load->max_neigh)[num_neigh-1];

    
                        //set functions and params for use below
                        SysFunction = &customMultilayer;
                        SysPara = &mlp;
                        hopfn = &multilayerTB;
                        
                        //(amended from bilayer to make connection profiles... (somewhat temporary.....)
                        blg_cnxpara.intra_thresh_min = (hop_to_load->NN_lowdis)[0];
                        blg_cnxpara.intra_thresh_max = (hop_to_load->NN_highdis)[0];
                        blg_cnxpara.inter_thresh_min = (hop_to_load->NN_lowdis)[1];
                        blg_cnxpara.inter_thresh_max = (hop_to_load->NN_highdis)[1];
                        blg_cnxpara.zthresh1 = (hop_to_load->NN_zmin)[1];
                        blg_cnxpara.zthresh2 = (hop_to_load->NN_zmax)[1];
                        
                        
                        sprintf(sysinfo, "L2_%d", length2); 
                  
		}
		
		
		//Bilayer device (using customMultilayer) combining a clean layer and a "Dots" type layer
		
                if(strcmp("BLDOTS", systemtype) == 0)
		{	
                        //multilayer connectivity
                        //hop_to_load = &MLG_NNTB_hop_params;
                        default_connect_rule = &mlg_conn_sep;
			defdouble_connect_rule = &mlg_conn_sep2;
			default_connection_params = &blg_cnxpara;
                        nngm=1;
                        num_neigh=nngm+1;
			
                    
                        //multilayer setup
                        mlp.num_layers = 2;
                        mlp.hoptype = 0;	//by default use simplest (BLG-like) hoppings
											//hoptype = 1 uses a more complex, Moon-Koshino-esque approach for twisted systems
                        mlp.length = createIntArray(mlp.num_layers);
                        mlp.length2 = createIntArray(mlp.num_layers);
                        mlp.geo = createIntArray(mlp.num_layers);
                        mlp.theta = createDoubleArray(mlp.num_layers);
                        mlp.delta = createNonSquareDoubleMatrix(mlp.num_layers, 3);
                        mlp.epsilon = createDoubleArray(mlp.num_layers);
                        mlp.origin = 1; //0 means (0,0), 1 means centre of l1 x l2 device
                        mlp.layerfn = (generatefn **)malloc(mlp.num_layers * sizeof(generatefn *));
                        mlp.layerpara = (void **)malloc(mlp.num_layers * sizeof(void *));
                        layertype=malloc(mlp.num_layers * sizeof(char *));
                        for(i=0; i< mlp.num_layers ; i++)
                            layertype[i] = (char *)malloc(40  * sizeof(char)); 
					/*allocates memory for whole matrix, 
					returned pointe r= pointer to 1st array*/
					
				
                        
                                        
                        //dotted layer setup
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
                        adotp.orientation = 0;
                        
                        
                        
                        
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
                                adotp.latgeo = latgeo;
                            }
                            if(strcmp("-dotgeo", argv[i]) == 0)
                            {
                                sscanf(argv[i+1], "%s", dotgeo);
                                adotp.dotgeo = dotgeo;
                            }
                            if(strcmp("-xyfluc", argv[i]) == 0)
                            {
                                sscanf(argv[i+1], "%lf", &(adotp.xyfluc));
                            }
                            if(strcmp("-radfluc", argv[i]) == 0)
                            {
                                sscanf(argv[i+1], "%lf", &(adotp.radfluc));
                            }
                            if(strcmp("-ADflip", argv[i]) == 0) //=0 "up", =1 "down", =2 random mix
                            {
                                sscanf(argv[i+1], "%d", &(adotp.orientation));
                            }
                            
                            // simple or moon-koshino hopping
                            if(strcmp("-mlhoptype", argv[i]) == 0)
                            {
                                    sscanf(argv[i+1], "%d", &(mlp.hoptype));
                            }
                        
                        }
		    
                        //system sizes calculated from lattice details
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
                        
                        if( strcmp("triAC", dotgeo) == 0 || strcmp("triZZ", dotgeo) == 0 )
                        {
                            if(adotp.orientation == 0)
                                sprintf(dot_info, "%s_dot_up_R_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (adotp.dotgeo), (adotp.AD_rad), (adotp.lat_width),  (adotp.lat_length), (adotp.xyfluc), (adotp.radfluc)); 
                            if(adotp.orientation == 1)
                                sprintf(dot_info, "%s_dot_down_R_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (adotp.dotgeo), (adotp.AD_rad), (adotp.lat_width),  (adotp.lat_length), (adotp.xyfluc), (adotp.radfluc)); 
                            if(adotp.orientation == 2)
                                sprintf(dot_info, "%s_dot_mix_R_%.1lf_%dx%d_xyf_%.1lf_rf_%.1lf", (adotp.dotgeo), (adotp.AD_rad), (adotp.lat_width),  (adotp.lat_length), (adotp.xyfluc), (adotp.radfluc)); 
                        }
		    
                        sprintf(sysinfo, "%s_%s", lattice_info, dot_info);
                                        

                        //Layer setup (defaults to AB, unrotated)
                        //AB stacking assumed -- can be changed
                        for(i=0; i< mlp.num_layers; i++)
                        {
                            (mlp.length)[i] = length1;
                            (mlp.length2)[i] = length2;
                            (mlp.geo)[i] = geo;
                            (mlp.theta)[i] = 0;
                            (mlp.epsilon)[i] = 1.0;

                            (mlp.delta)[i][0] = 0.0;
                            (mlp.delta)[i][1] = 0.0 + (1/sqrt(3)) * (i %2);
                            (mlp.delta)[i][2] = i * 1.362;
                        }
                        
                        (mlp.layerfn)[0] = &simpleRibbonGeo;
                        (mlp.layerpara)[0]  = NULL;
                        sprintf(layertype[0], "CLEAN");
                        
                        (mlp.layerfn)[1] = &genDotDevice;
                        (mlp.layerpara)[1]  = &adotp;
                        sprintf(layertype[1], "DOTS");
                            
                        
                        
                 
                        
                        //check for command line arguments which vary these
                      
                        for(i=1; i<argc-1; i++)
                        {   
                            for(j=0; j<mlp.num_layers; j++)
                            {
                                sprintf(temp_in_string, "-lay%ddeltax", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &(mlp.delta)[j][0]);
                                }
                                sprintf(temp_in_string, "-lay%ddeltay", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &(mlp.delta)[j][1]);
                                }
                                sprintf(temp_in_string, "-lay%ddeltaz", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &(mlp.delta)[j][2]);
                                }
                                sprintf(temp_in_string, "-lay%dtheta", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &(mlp.theta)[j]);
                                }
                                sprintf(temp_in_string, "-lay%deps", j);
                                if(strcmp(temp_in_string, argv[i]) == 0)
                                {
                                        sscanf(argv[i+1], "%lf", &(mlp.epsilon)[j]);
                                }
                                
                                
                            }
                            //center of rotation and common origin of layers
                            if(strcmp("-mlorigin", argv[i]) == 0)
                            {
                                    sscanf(argv[i+1], "%d", &(mlp.origin));
                            }
                        
                        }
                        
                        
                        if (mlp.hoptype == 1)
							hop_to_load = &MLG_NNTB_hop_params;
								
						if (mlp.hoptype == 0)
							hop_to_load = &MLG2_NNTB_hop_params;
							
						max_neigh=(hop_to_load->max_neigh)[num_neigh-1];
    
                        //set functions and params for use below
                        SysFunction = &customMultilayer;
                        SysPara = &mlp;
                        hopfn = &multilayerTB;
                        
                        //(amended from bilayer to make connection profiles... (somewhat temporary.....)
                        blg_cnxpara.intra_thresh_min = (hop_to_load->NN_lowdis)[0];
                        blg_cnxpara.intra_thresh_max = (hop_to_load->NN_highdis)[0];
                        blg_cnxpara.inter_thresh_min = (hop_to_load->NN_lowdis)[1];
                        blg_cnxpara.inter_thresh_max = (hop_to_load->NN_highdis)[1];
                        blg_cnxpara.zthresh1 = (hop_to_load->NN_zmin)[1];
                        blg_cnxpara.zthresh2 = (hop_to_load->NN_zmax)[1];
                        
                        
                        //sprintf(sysinfo, "L2_%d", length2); 
                  
		}
	  
	  
	  
	  
	  //END OF DEVICE DEFINITIONS!!
	  
	  
	  
	 
	  
	  
	  
	  
	  
	  
	  double kxmin=0.0, kxmax =(2*M_PI/length2);
		for(i=1; i<argc-1; i++)
		{
			if(strcmp("-kxmin", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%lf", &kxmin);
			}
			if(strcmp("-kxmax", argv[i]) == 0)
			{
				sscanf(argv[i+1], "%lf", &kxmax);
			}

		}
	  
	  
	//misc  
	  
	  double cond2;
	  
	  
	  double startloc, endloc, starty, endy;
	  
	  //approximate size values for lead generation...
		if(geo==0)
		{
			endloc = length2*1.0 -0.5;
			startloc = 0.0;
			starty = 0.0;
			endy = length1*sqrt(3)/2;
		}
		if(geo==1)
		{
			endloc = length2*sqrt(3) - (1/sqrt(3));
			startloc= 0.0;
			starty = -0.5;
			endy = length1*1.0/2;
		}
		
	  
	   //ADDITIONAL EDGE POTENTIALS, DISORDER ETC
	  
	  
	  
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
		
		 sprintf(disinfo, "POTDIS_c_%.4lf_d_%.3lf_xi_%.1lf", pdconc, pddelta, pdxi);
	      }
	      
	      
	//potential profiles near edges
		int cedge_leads=0, cedge_posonly=0;
		int cedge_type=0;
		char cedgeconf[40], cedgetypename[40];
		//strength (1) and coefficients (2,3,4) for each sublattice at top T and bottom B edges
		//(2) generally shows the extent of the potential (decay constant, etc)
		//(3) is generally a distance modifying the edge position from the edgemost atom to elsewhere
		double cedge_AT1=0.0, cedge_AT2=0, cedge_AT3=0.0, cedge_AT4=0;
		double cedge_BT1=0.0, cedge_BT2=0, cedge_BT3=0.0, cedge_BT4=0;
		double cedge_AB1=0.0, cedge_AB2=0, cedge_AB3=0.0, cedge_AB4=0;
		double cedge_BB1=0.0, cedge_BB2=0, cedge_BB3=0.0, cedge_BB4=0;
		
			sprintf(cedgeconf, "misc");
			sprintf(cedgetypename, "");
		
		cedgepot_para cedgepara = {};
		
		if(cedge==1)
		{
			for(i=1; i<argc-1; i++)
			{
				if(strcmp("-celeads", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &cedge_leads);
				}
				if(strcmp("-ceposonly", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &cedge_posonly);
				}
				if(strcmp("-cetype", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &cedge_type);
				}
				//shorthand arguments for multiple params
				//sets that there is no dependence on particlar sublattice or edge
				if(strcmp("-ce1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AT1);
					cedge_BT1 = cedge_AT1;
					cedge_AB1 = cedge_AT1;
					cedge_BB1 = cedge_AT1;
				}
				if(strcmp("-ce2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AT2);
					cedge_BT2 = cedge_AT2;
					cedge_AB2 = cedge_AT2;
					cedge_BB2 = cedge_AT2;
				}
				if(strcmp("-ce3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AT3);
					cedge_BT3 = cedge_AT3;
					cedge_AB3 = cedge_AT3;
					cedge_BB3 = cedge_AT3;
				}
				if(strcmp("-ce4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AT4);
					cedge_BT4 = cedge_AT4;
					cedge_AB4 = cedge_AT4;
					cedge_BB4 = cedge_AT4;
				}
				
				//individual sublattice/edge params
				if(strcmp("-ceAT1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AT1);
				}
				if(strcmp("-ceAT2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AT2);
				}
				if(strcmp("-ceAT3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AT3);
				}
				if(strcmp("-ceAT4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AT4);
				}
				
				if(strcmp("-ceBT1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_BT1);
				}
				if(strcmp("-ceBT2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_BT2);
				}
				if(strcmp("-ceBT3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_BT3);
				}
				if(strcmp("-ceBT4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_BT4);
				}
				
				if(strcmp("-ceAB1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AB1);
				}
				if(strcmp("-ceAB2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AB2);
				}
				if(strcmp("-ceAB3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AB3);
				}
				if(strcmp("-ceAB4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_AB4);
				}
				
				if(strcmp("-ceBB1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_BB1);
				}
				if(strcmp("-ceBB2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_BB2);
				}
				if(strcmp("-ceBB3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_BB3);
				}
				if(strcmp("-ceBB4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &cedge_BB4);
				}
				
				if(strcmp("-cename", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%s", cedgeconf);
				}
				
				if(cedge_type==0)
					sprintf(cedgetypename, "cnst");  //const pot within distance of edge
					
				if(cedge_type==1)
					sprintf(cedgetypename, "expdcy");  //power law from edge
				
				if(cedge_type==2)
					sprintf(cedgetypename, "pwrlaw");  //exponential decay from edge
					
				if(cedge_type==3)
					sprintf(cedgetypename, "lindrop");  //power law from edge
				
				if(cedge_type==4)
					sprintf(cedgetypename, "efetov");  //exponential decay from edge	
				
			
			}
		
			cedgepara.type=cedge_type;
			cedgepara.AT1=cedge_AT1;
			cedgepara.AT2=cedge_AT2;
			cedgepara.AT3=cedge_AT3;
			cedgepara.AT4=cedge_AT4; 
			cedgepara.BT1=cedge_BT1;
			cedgepara.BT2=cedge_BT2;
			cedgepara.BT3=cedge_BT3;
			cedgepara.BT4=cedge_BT4;
			cedgepara.AB1=cedge_AB1;
			cedgepara.AB2=cedge_AB2;
			cedgepara.AB3=cedge_AB3;
			cedgepara.AB4=cedge_AB4;
			cedgepara.BB1=cedge_BB1;
			cedgepara.BB2=cedge_BB2;
			cedgepara.BB3=cedge_BB3;
			cedgepara.BB4=cedge_BB4;
			
			sprintf(cedgeinfo, "CEDGE_%s_%s", cedgetypename, cedgeconf);
			
			if(cedge_leads==0)
				sprintf(cedgeinfo, "CEDGEnl_%s_%s", cedgetypename, cedgeconf);
		}
	      
	      
	      
	//secondary potential profiles near edges
		int dedge_leads=0, dedge_posonly=0;
		int dedge_type=0;
		char dedgeconf[40], dedgetypename[40];
		//strength (1) and coefficients (2,3,4) for each sublattice at top T and bottom B edges
		//(2) generally shows the extent of the potential (decay constant, etc)
		//(3) is generally a distance modifying the edge position from the edgemost atom to elsewhere
		double dedge_AT1=0.0, dedge_AT2=0, dedge_AT3=0.0, dedge_AT4=0;
		double dedge_BT1=0.0, dedge_BT2=0, dedge_BT3=0.0, dedge_BT4=0;
		double dedge_AB1=0.0, dedge_AB2=0, dedge_AB3=0.0, dedge_AB4=0;
		double dedge_BB1=0.0, dedge_BB2=0, dedge_BB3=0.0, dedge_BB4=0;
		
			sprintf(dedgeconf, "misc");
			sprintf(dedgetypename, "");
		
		cedgepot_para dedgepara = {};
		
		if(dedge==1)
		{
			for(i=1; i<argc-1; i++)
			{
				if(strcmp("-deleads", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &dedge_leads);
				}
				if(strcmp("-deposonly", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &dedge_posonly);
				}
				if(strcmp("-detype", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &dedge_type);
				}
				//shorthand arguments for multiple params
				//sets that there is no dependence on particlar sublattice or edge
				if(strcmp("-de1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AT1);
					dedge_BT1 = dedge_AT1;
					dedge_AB1 = dedge_AT1;
					dedge_BB1 = dedge_AT1;
				}
				if(strcmp("-de2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AT2);
					dedge_BT2 = dedge_AT2;
					dedge_AB2 = dedge_AT2;
					dedge_BB2 = dedge_AT2;
				}
				if(strcmp("-de3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AT3);
					dedge_BT3 = dedge_AT3;
					dedge_AB3 = dedge_AT3;
					dedge_BB3 = dedge_AT3;
				}
				if(strcmp("-de4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AT4);
					dedge_BT4 = dedge_AT4;
					dedge_AB4 = dedge_AT4;
					dedge_BB4 = dedge_AT4;
				}
				
				//individual sublattice/edge params
				if(strcmp("-deAT1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AT1);
				}
				if(strcmp("-deAT2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AT2);
				}
				if(strcmp("-deAT3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AT3);
				}
				if(strcmp("-deAT4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AT4);
				}
				
				if(strcmp("-deBT1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_BT1);
				}
				if(strcmp("-deBT2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_BT2);
				}
				if(strcmp("-deBT3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_BT3);
				}
				if(strcmp("-deBT4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_BT4);
				}
				
				if(strcmp("-deAB1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AB1);
				}
				if(strcmp("-deAB2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AB2);
				}
				if(strcmp("-deAB3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AB3);
				}
				if(strcmp("-deAB4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_AB4);
				}
				
				if(strcmp("-deBB1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_BB1);
				}
				if(strcmp("-deBB2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_BB2);
				}
				if(strcmp("-deBB3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_BB3);
				}
				if(strcmp("-deBB4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &dedge_BB4);
				}
				
				if(strcmp("-dename", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%s", dedgeconf);
				}
				
				if(dedge_type==0)
					sprintf(dedgetypename, "cnst");  //const pot within distance of edge
					
				if(dedge_type==1)
					sprintf(dedgetypename, "expdcy");  //power law from edge
				
				if(dedge_type==2)
					sprintf(dedgetypename, "pwrlaw");  //exponential decay from edge
					
				if(dedge_type==3)
					sprintf(dedgetypename, "lindrop");  //power law from edge
				
				if(dedge_type==4)
					sprintf(dedgetypename, "efetov");  //exponential decay from edge	
				
			
			}
		
			dedgepara.type=dedge_type;
			dedgepara.AT1=dedge_AT1;
			dedgepara.AT2=dedge_AT2;
			dedgepara.AT3=dedge_AT3;
			dedgepara.AT4=dedge_AT4; 
			dedgepara.BT1=dedge_BT1;
			dedgepara.BT2=dedge_BT2;
			dedgepara.BT3=dedge_BT3;
			dedgepara.BT4=dedge_BT4;
			dedgepara.AB1=dedge_AB1;
			dedgepara.AB2=dedge_AB2;
			dedgepara.AB3=dedge_AB3;
			dedgepara.AB4=dedge_AB4;
			dedgepara.BB1=dedge_BB1;
			dedgepara.BB2=dedge_BB2;
			dedgepara.BB3=dedge_BB3;
			dedgepara.BB4=dedge_BB4;
			
			sprintf(dedgeinfo, "CEDGE_%s_%s", dedgetypename, dedgeconf);
			
			if(dedge_leads==0)
				sprintf(dedgeinfo, "CEDGEnl_%s_%s", dedgetypename, dedgeconf);
		}
	  
	  
	//tertiary potential profiles near edges
		int eedge_leads=0, eedge_posonly=0;
		int eedge_type=0;
		char eedgeconf[40], eedgetypename[40];
		//strength (1) and coefficients (2,3,4) for each sublattice at top T and bottom B edges
		//(2) generally shows the extent of the potential (decay constant, etc)
		//(3) is generally a distance modifying the edge position from the edgemost atom to elsewhere
		double eedge_AT1=0.0, eedge_AT2=0, eedge_AT3=0.0, eedge_AT4=0;
		double eedge_BT1=0.0, eedge_BT2=0, eedge_BT3=0.0, eedge_BT4=0;
		double eedge_AB1=0.0, eedge_AB2=0, eedge_AB3=0.0, eedge_AB4=0;
		double eedge_BB1=0.0, eedge_BB2=0, eedge_BB3=0.0, eedge_BB4=0;
		
			sprintf(eedgeconf, "misc");
			sprintf(eedgetypename, "");
		
		cedgepot_para eedgepara = {};
		
		if(eedge==1)
		{
			for(i=1; i<argc-1; i++)
			{
				if(strcmp("-eeleads", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &eedge_leads);
				}
				if(strcmp("-eeposonly", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &eedge_posonly);
				}
				if(strcmp("-eetype", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &eedge_type);
				}
				//shorthand arguments for multiple params
				//sets that there is no dependence on particlar sublattice or edge
				if(strcmp("-ee1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AT1);
					eedge_BT1 = eedge_AT1;
					eedge_AB1 = eedge_AT1;
					eedge_BB1 = eedge_AT1;
				}
				if(strcmp("-ee2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AT2);
					eedge_BT2 = eedge_AT2;
					eedge_AB2 = eedge_AT2;
					eedge_BB2 = eedge_AT2;
				}
				if(strcmp("-ee3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AT3);
					eedge_BT3 = eedge_AT3;
					eedge_AB3 = eedge_AT3;
					eedge_BB3 = eedge_AT3;
				}
				if(strcmp("-ee4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AT4);
					eedge_BT4 = eedge_AT4;
					eedge_AB4 = eedge_AT4;
					eedge_BB4 = eedge_AT4;
				}
				
				//individual sublattice/edge params
				if(strcmp("-eeAT1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AT1);
				}
				if(strcmp("-eeAT2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AT2);
				}
				if(strcmp("-eeAT3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AT3);
				}
				if(strcmp("-eeAT4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AT4);
				}
				
				if(strcmp("-eeBT1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_BT1);
				}
				if(strcmp("-eeBT2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_BT2);
				}
				if(strcmp("-eeBT3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_BT3);
				}
				if(strcmp("-eeBT4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_BT4);
				}
				
				if(strcmp("-eeAB1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AB1);
				}
				if(strcmp("-eeAB2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AB2);
				}
				if(strcmp("-eeAB3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AB3);
				}
				if(strcmp("-eeAB4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_AB4);
				}
				
				if(strcmp("-eeBB1", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_BB1);
				}
				if(strcmp("-eeBB2", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_BB2);
				}
				if(strcmp("-eeBB3", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_BB3);
				}
				if(strcmp("-eeBB4", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%lf", &eedge_BB4);
				}
				
				if(strcmp("-eename", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%s", eedgeconf);
				}
				
				if(eedge_type==0)
					sprintf(eedgetypename, "cnst");  //const pot within distance of edge
					
				if(eedge_type==1)
					sprintf(eedgetypename, "expdcy");  //power law from edge
				
				if(eedge_type==2)
					sprintf(eedgetypename, "pwrlaw");  //exponential decay from edge
					
				if(eedge_type==3)
					sprintf(eedgetypename, "lindrop");  //power law from edge
				
				if(eedge_type==4)
					sprintf(eedgetypename, "efetov");  //exponential decay from edge	
				
			
			}
		
			eedgepara.type=eedge_type;
			eedgepara.AT1=eedge_AT1;
			eedgepara.AT2=eedge_AT2;
			eedgepara.AT3=eedge_AT3;
			eedgepara.AT4=eedge_AT4; 
			eedgepara.BT1=eedge_BT1;
			eedgepara.BT2=eedge_BT2;
			eedgepara.BT3=eedge_BT3;
			eedgepara.BT4=eedge_BT4;
			eedgepara.AB1=eedge_AB1;
			eedgepara.AB2=eedge_AB2;
			eedgepara.AB3=eedge_AB3;
			eedgepara.AB4=eedge_AB4;
			eedgepara.BB1=eedge_BB1;
			eedgepara.BB2=eedge_BB2;
			eedgepara.BB3=eedge_BB3;
			eedgepara.BB4=eedge_BB4;
			
			sprintf(eedgeinfo, "CEDGE_%s_%s", eedgetypename, eedgeconf);
			
			if(eedge_leads==0)
				sprintf(eedgeinfo, "CEDGEnl_%s_%s", eedgetypename, eedgeconf);
		}
		
	  
	  connectrules = default_connect_rule;
	  
	  //standard hall settings
	   if(ishallbar==1 || ishallbar==2)
            {
            connectrules = default_connect_rule;
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
                
                if(ishallbar==1)
                {
                    sprintf(peritype, "HALLBAR_%d_%d_rw_%d", ntop, nbot, hall_denom);
                }
                if(ishallbar==2)
                {
                    sprintf(peritype, "SHBAR_%d_%d_rw_%d", ntop, nbot, hall_denom);
                }
                

                
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
            
            
            double *leadsxmin, *leadsxmax;
            double leadwidth, endsep, endpos, convunit;
            int typeleads;
	    int currIn, currOut;
	    int *leadOrder;
        
            //settings for (generic) 4-probe configuration
            //these leads are generally trivial self energies, rather than self consistent SGF based
            //this should work for metallic (FM) contacts, STM tips etc with some modification
            //this setup is based on simple metallic leads (typeleads 0), attaching across the entire device width within an x range.
            
            //positions at endpos from start and end, and at endpos +  endsep from each.
            //make sure system is long enough to support this
            //leads are of leadwidth units wide
            if(ishallbar ==3)
            {
                num_leads=4; //default - checks again below!
                connectrules = default_connect_rule;
                isperiodic=0;
                cnxpara.periodic = 0;
                kpts = 1;
		currIn=1;
		currOut=0;
                
	      
                //metallic leads are default for now...
                //x direction divided in 20)
                metal_leads=1;
                
                
                //default lead geometries
                    convunit = (endloc-startloc)/20;
                    
                    leadwidth = convunit;
                    endpos = convunit;
                    endsep = 3*convunit;
                    
                
                
                
                for(i=1; i<argc-1; i++)
                {
                    if(strcmp("-numleads", argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%d", &num_leads);
                    }
                    if(strcmp("-mptypeleads", argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%d", &typeleads);
                    }
                    if(strcmp("-mpleadwidth", argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%lf", &leadwidth);
                    }
                    if(strcmp("-mpendsep", argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%lf", &endsep);
                    }
                    if(strcmp("-mpcurrin", argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%d", &currIn);
                    }
                    if(strcmp("-mpcurrout", argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%d", &currOut);
                    }
                }
                
                if(num_leads ==4)
                {
                        leadsxmin = createDoubleArray(4);
                        leadsxmax = createDoubleArray(4);
                        
                        leadsxmin[0] = endpos;
                        leadsxmax[0] = endpos + leadwidth;
                        leadsxmin[1] = endpos + endsep;
                        leadsxmax[1] = endpos + endsep + leadwidth;
                        
                        leadsxmin[2] = endloc -endpos -endsep -leadwidth;
                        leadsxmax[2] = endloc -endpos -endsep;
                        leadsxmin[3] = endloc - endpos - leadwidth;
                        leadsxmax[3] = endloc - endpos;
                        
                    
                }
                
                leadOrder = createIntArray(num_leads);
		
		leadOrder[0] = currIn;
		leadOrder[num_leads-1] = currOut;
		j=1;
		
		for(i=0; i<num_leads; i++)
		{
			if(i != currIn && i != currOut)
			{
				leadOrder[j] = i;
				j++;
			}
		}
                
//                 for(i=0; i<num_leads; i++)
// 			printf("# lead order %d %d\n", i, leadOrder[i]);
                
                sprintf(peritype, "METAL_%d_PROBES_%d_to_%d_w%.0lf_s%.0lf_%.0lf", num_leads, currIn, currOut, leadwidth, endpos, endsep);
                
            }
            
		gen_hop_params metal_hop_p = {};
		multix_start_params mxsp = {};
		double metal_alpha=0.0, metal_sig=-1.0, metal_beta=1.0, metal_hop=t, metal_default_width=4;;
		
		
		 
                for(i=1; i<argc-1; i++)
                {
                    if(strcmp("-metalalpha", argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%lf", &metal_alpha);
                    }
                    if(strcmp("-metalbeta", argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%lf", &metal_beta);
                    }
                    if(strcmp("-metalsig", argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%lf", &metal_sig);
                    }
                    if(strcmp("-metalhop", argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%lf", &metal_hop);
                    }
                    
                }
		
		
		metal_hop_p.hops = createCompArray(4);
		metal_hop_p.hops[0] = metal_sig;
                metal_hop_p.hops[1] = metal_alpha;
                metal_hop_p.hops[2] = metal_beta;
                metal_hop_p.hops[3] = metal_hop;
	    
            if(metal_leads == 1)
            {
                starting_cell_mode=4;
            
               
                
                metal_hop_p.hops[0] = metal_sig;
                metal_hop_p.hops[1] = metal_alpha;
                metal_hop_p.hops[2] = metal_beta;
                metal_hop_p.hops[3] = metal_hop;
                
                mxsp.num_leads=4;
                mxsp.startx=leadsxmin;
                mxsp.endx=leadsxmax;
                
            }
	  


	  
		
		char leadconf[40], additional_lead_info[40];
		multiple_para *mleadps[num_leads];
		custom_start_params cstart_p ={};
		cstart_p.leadtype = createIntArray(num_leads);
		sprintf(leadconf, "CUSTOM");
                
                if(ispatched !=0)
                {
                    sprintf(leadconf, "PGF");
                }

		sprintf(additional_lead_info, "");

		int nleft, nright, nfull;
		int counttop, countbot, countleft, countright, countfull;
		double metaldim2=2.0;
		if(ishallbar==4)  //CUSTOM LEADS MODE
		{
			currIn=1;
			currOut=0;
			connectrules = default_connect_rule;
			for(i=1; i<argc-1; i++)
			{
				//uses main code numleads check
// 				if(strcmp("-numleads", argv[i]) == 0)
// 				{
// 					sscanf(argv[i+1], "%d", &num_leads);
// 				}
				if(strcmp("-mpcurrin", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &currIn);
				}
				if(strcmp("-mpcurrout", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%d", &currOut);
				}
				if(strcmp("-leadconfname", argv[i]) == 0)
				{
					sscanf(argv[i+1], "%s", leadconf);
				}
			}
			
			counttop=0; countbot=0; countleft=0; countright=0, countfull=0;
			ntop=(int)(num_leads-2 - (num_leads-2)/2 );
			nbot=(num_leads-2)/2;
			nleft=1;
			nright=1;
			
			leadOrder = createIntArray(num_leads);
		
			leadOrder[0] = currIn;
			leadOrder[num_leads-1] = currOut;
			j=1;
			
			for(i=0; i<num_leads; i++)
			{
				if(i != currIn && i != currOut)
				{
					leadOrder[j] = i;
					j++;
				}
			}
			
			for(j=0; j<num_leads; j++)
			{
				mleadps[j] = (multiple_para *)malloc(sizeof(multiple_para));
 				(mleadps[j]->name) = (char *)malloc(sizeof(char) * (41));
				sprintf(mleadps[j]->name, "RIBBON"); //default is ribbon
				cstart_p.leadtype[j]=0;
				

								
				if(j==0)
					(mleadps[j]->def_pos)=0;
				if(j==1)
					(mleadps[j]->def_pos)=1;
				if(j>1 && j<2+ ntop)
					(mleadps[j]->def_pos)=2;
				if(j>= 2+ (int)(num_leads-2 - (num_leads-2)/2 ))
					(mleadps[j]->def_pos)=3;
				
				sprintf(temp_in_string, "-lead%dside", j);
				for(i=1; i<argc-1; i++)
				{
					if(strcmp(temp_in_string, argv[i]) == 0)
					{
						sscanf(argv[i+1], "%d", &(mleadps[j]->def_pos));
					}
				}
				
				if((mleadps[j]->def_pos)==0)
					countleft++;
				if((mleadps[j]->def_pos)==1)
					countright++;
				if((mleadps[j]->def_pos)==2)
					counttop++;
				if((mleadps[j]->def_pos)==3)
					countbot++;
				if((mleadps[j]->def_pos)==4)
					countfull++;
			}
			ntop=counttop; nbot=countbot; nleft=countleft; nright=countright, nfull = countfull;
			counttop=0; countbot=0; countleft=0; countright=0, countfull=0;
			
			
			//individual lead types and default settings
			for(j=0; j<num_leads; j++)
			{
				sprintf(temp_in_string, "-lead%dtype", j);
				
				//set lead types and default settings
				for(i=1; i<argc-1; i++)
				{
					if(strcmp(temp_in_string, argv[i]) == 0)
					{
						if(strcmp("RIBBON", argv[i+1]) == 0)
						{
							//lead j is a ribbon
							sprintf(mleadps[j]->name, "RIBBON");
							cstart_p.leadtype[j]=0;
						}
						
						if(strcmp("METALX", argv[i+1]) == 0)
						{
							//lead j is a metallic strip with finite xdim
							sprintf(mleadps[j]->name, "METALX");
							cstart_p.leadtype[j]=1;
						}
						
						if(strcmp("METALXY", argv[i+1]) == 0)
						{
							//lead j is a metallic strip with finite x&y dims
							sprintf(mleadps[j]->name, "METALXY");
							cstart_p.leadtype[j]=1;
						}
						
						if(strcmp("STM", argv[i+1]) == 0)
						{
							//lead j is a simple STM tip
							sprintf(mleadps[j]->name, "STM");
							cstart_p.leadtype[j]=2;
						}
						
						if(strcmp("PATCHED", argv[i+1]) == 0)
						{
							//lead j is a patched boundary
							sprintf(mleadps[j]->name, "PATCHED");
							cstart_p.leadtype[j]=3;
						}
						
					}
				}
				
				
				//set other lead properties
				if(strcmp("RIBBON", mleadps[j]->name) == 0)
				{
					(mleadps[j]->indiv_lead_para) = (rib_lead_para *)malloc(sizeof(rib_lead_para));
					
					(mleadps[j]->indiv_lead_fn) = &singleRibbonLead;
					(mleadps[j]->indiv_gen_fn) = &genSingleRibbonLead;
					((rib_lead_para *)(mleadps[j]->indiv_lead_para))->hopfn = hopfn;
					((rib_lead_para *)(mleadps[j]->indiv_lead_para))->hoppara = &hoppara;
					((rib_lead_para *)(mleadps[j]->indiv_lead_para))->homo_A_pot = 0.0;
					((rib_lead_para *)(mleadps[j]->indiv_lead_para))->homo_B_pot = 0.0;
					
					if (xleadsavgpot == 1)
					{
						((rib_lead_para *)(mleadps[j]->indiv_lead_para))->homo_A_pot = suba_conc * suba_pot;
						((rib_lead_para *)(mleadps[j]->indiv_lead_para))->homo_B_pot = subb_conc * subb_pot;
					}

					
					//default sizes and geos
						if((mleadps[j]->def_pos)==0)
						{
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->width=(length1 / (nleft)) - 2*(1 - (int)(nright/(nleft*nleft)) );
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->geo=geo;
						}
						if((mleadps[j]->def_pos)==1)
						{
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->width=(length1 / (nright)) - 2*(1 - (int)(nright/(nright*nright)) );
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->geo=geo;
							
						}
						if((mleadps[j]->def_pos)==2)
						{
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->width= (int) 2*length2 / (3*ntop);
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->geo=1-geo;
						}
						if((mleadps[j]->def_pos)==3)
						{
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->width= (int) 2*length2 / (3*nbot);
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->geo=1-geo;
						}
						
					//command line size and geo options
						for(i=1; i<argc-1; i++)
						{
							sprintf(temp_in_string, "-lead%dsize", j);
							if(strcmp(temp_in_string, argv[i]) == 0)
							{
								sscanf(argv[i+1], "%d", &((rib_lead_para *)(mleadps[j]->indiv_lead_para))->width);
							}
							sprintf(temp_in_string, "-lead%dgeo", j);
							if(strcmp(temp_in_string, argv[i]) == 0)
							{
								sscanf(argv[i+1], "%d", &(((rib_lead_para *)(mleadps[j]->indiv_lead_para))->geo));
							}
						}
					
					//default positioning
						if((mleadps[j]->def_pos)==0)
						{
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord = HallPositioning(length1, nleft, countleft, 0, 1, ((rib_lead_para *)(mleadps[j]->indiv_lead_para))->width);
							countleft++;
						}
						if((mleadps[j]->def_pos)==1)
						{
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord = HallPositioning(length1, nright, countright, 0, 1, ((rib_lead_para *)(mleadps[j]->indiv_lead_para))->width);
							countright++;
							
						}
						if((mleadps[j]->def_pos)==2)
						{
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord = HallPositioning(length2, ntop, counttop, buffer_rows, 2, ((rib_lead_para *)(mleadps[j]->indiv_lead_para))->width);
							counttop++;
						}
						if((mleadps[j]->def_pos)==3)
						{
							((rib_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord = HallPositioning(length2, nbot, countbot, buffer_rows, 2, ((rib_lead_para *)(mleadps[j]->indiv_lead_para))->width);
							countbot++;
						}
						
					//command line positioning options
						for(i=1; i<argc-1; i++)
						{
							sprintf(temp_in_string, "-lead%dpos", j);
							if(strcmp(temp_in_string, argv[i]) == 0)
							{
								sscanf(argv[i+1], "%d", &((rib_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord);
							}
							
						}
						
					((rib_lead_para *)(mleadps[j]->indiv_lead_para))->def_pos=(mleadps[j]->def_pos);	
					
					printf("#lead %d: %s\n", j, mleadps[j]->name);
					printf("#lead %d, defpos %d, width %d, start %d\n", j, (mleadps[j]->def_pos),  ((rib_lead_para *)(mleadps[j]->indiv_lead_para))->width, ((rib_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord);
				
				}
				
				
				
				if(strcmp("METALX", mleadps[j]->name) == 0 || strcmp("METALXY", mleadps[j]->name) == 0)
				{
					(mleadps[j]->indiv_lead_para) = (metal_lead_para *)malloc(sizeof(metal_lead_para));
					
					(mleadps[j]->indiv_lead_fn) = &singleSimplestMetalLead;
					(mleadps[j]->indiv_gen_fn) = &genSingleMetalLead;
					
												
					//default hopping paramaters
						((metal_lead_para *)(mleadps[j]->indiv_lead_para))->hoppara = (gen_hop_params *)malloc(sizeof(gen_hop_params));
						((gen_hop_params*)((metal_lead_para *)(mleadps[j]->indiv_lead_para))->hoppara)->hops = createCompArray(4);
						((gen_hop_params*)((metal_lead_para *)(mleadps[j]->indiv_lead_para))->hoppara)->hops[0] = metal_sig;
						((gen_hop_params*)((metal_lead_para *)(mleadps[j]->indiv_lead_para))->hoppara)->hops[1] = metal_alpha;
						((gen_hop_params*)((metal_lead_para *)(mleadps[j]->indiv_lead_para))->hoppara)->hops[2] = metal_beta;
						((gen_hop_params*)((metal_lead_para *)(mleadps[j]->indiv_lead_para))->hoppara)->hops[3] = metal_hop;
					
						
						
					//default sizes
						if((mleadps[j]->def_pos)==0)
						{
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width2=(length1 / (nleft)) - 2*(1 - (int)(nright/(nleft*nleft)) );
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width = metal_default_width;
						}
						if((mleadps[j]->def_pos)==1)
						{
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width2=(length1 / (nright)) - 2*(1 - (int)(nright/(nright*nright)) );
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width = metal_default_width;
						}
						if((mleadps[j]->def_pos)==2)
						{
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width= (int) 2*length2 / (3*ntop);
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width2 = metal_default_width;
						}
						if((mleadps[j]->def_pos)==3)
						{
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width= (int) 2*length2 / (3*nbot);
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width2 = metal_default_width;
						}

						
						//special settings for metalX!!
						if(strcmp("METALX", mleadps[j]->name) == 0)
						{
							(mleadps[j]->def_pos)=4;
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width= ((int) 2*length2 / (5*num_leads))*1.0;
// 							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width2= ((int) 2*length2 / (3*num_leads))*1.0;
						}
						
						
						
					//command line size option	
						for(i=1; i<argc-1; i++)
						{
							sprintf(temp_in_string, "-lead%dsize", j);
							if(strcmp(temp_in_string, argv[i]) == 0)
							{
								sscanf(argv[i+1], "%lf", &((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width);
							}
							sprintf(temp_in_string, "-lead%dsize2", j);
							if(strcmp(temp_in_string, argv[i]) == 0)
							{
								sscanf(argv[i+1], "%lf", &((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width2);
							}
							
						}
						
					//default positioning (defpos5 gives absolute positioning, specified from command line)
						if((mleadps[j]->def_pos)==0)
						{
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord2 = (double) HallPositioning(length1, nleft, countleft, 0, 1, (int) ( ((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width2));
							countleft++;
						}
						if((mleadps[j]->def_pos)==1)
						{
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord2 = (double) HallPositioning(length1, nright, countright, 0, 1, (int) ( ((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width2));
							countright++;
							
						}
						if((mleadps[j]->def_pos)==2)
						{
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord = (double) HallPositioning(length2, ntop, counttop, buffer_rows, 2, (int) ( ((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width));
							counttop++;
						}
						if((mleadps[j]->def_pos)==3)
						{
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord = (double) HallPositioning(length2, nbot, countbot, buffer_rows, 2, (int) ( ((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width));
							countbot++;
						}
						
						if((mleadps[j]->def_pos)==4)
						{
							((metal_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord = (double) HallPositioning(length2, num_leads, countfull, buffer_rows, 2, (int) ( ((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width));
							countfull++;
						}
						
					//command line positioning options
						for(i=1; i<argc-1; i++)
						{
							sprintf(temp_in_string, "-lead%dpos", j);
							if(strcmp(temp_in_string, argv[i]) == 0)
							{
								sscanf(argv[i+1], "%lf", &((metal_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord);
							}
							sprintf(temp_in_string, "-lead%dpos2", j);
							if(strcmp(temp_in_string, argv[i]) == 0)
							{
								sscanf(argv[i+1], "%lf", &((metal_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord2);
							}
							
						}
						
					((metal_lead_para *)(mleadps[j]->indiv_lead_para))->def_pos=(mleadps[j]->def_pos);
					
					printf("#lead %d: %s\n", j, mleadps[j]->name);
					printf("#lead %d, defpos %d, width %lf %lf, start %lf %lf\n", j, (mleadps[j]->def_pos),  ((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width, ((metal_lead_para *)(mleadps[j]->indiv_lead_para))->width2, ((metal_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord,  ((metal_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord2);
				
				}
				
						
				
			
				
			
				
			}
// 			printf("%d	%d\n", ntop, nbot);
			for(j=0; j<num_leads; j++)
			{
				

			}
			
			sprintf(peritype, "%s_%d_LEADS_%d_to_%d_%s", leadconf, num_leads, currIn, currOut, additional_lead_info);
			
		}
		
	  
	  //loop info
	      double loop_step, loop_min_temp, loopmax, loopmin;
	      
	      if(strcmp("E", loop_type) == 0)
	      {
		Bfield = Bmin;
		loopmin = Emin;
		loopmax = Emax;
		sprintf(loopinfo, "Eloop_%+.2lf_to_%+.2lf_Bfixed_%+.3lf", Emin, Emax, Bfield);
		sprintf(loopmininfo, "E_%+.3lf_B_%+.3lf", Emin, Bfield);
	      }

	      if(strcmp("B", loop_type) == 0)
	      {
		realE = Emin;
		loopmin = Bmin;
		loopmax = Bmax;
		sprintf(loopinfo, "Bloop_%+.2lf_to_%+.2lf_Efixed_%+.3lf", Bmin, Bmax, realE);
		sprintf(loopmininfo, "E_%+.3lf_B_%+.3lf", realE, Bmin);


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
		else if (vgtype ==3)
		{
		  sprintf(vgtname, "j2efetov");
		}
	
		
		sprintf(loopinfo, "VG_%s_loop_%+.2lf_to_%+.2lf_Bfixed_%+.3lf_Efixed_%+.3lf", vgtname, VGmin, VGmax, Bfield, realE);
		sprintf(loopmininfo, "VG_%s_%+.2lf_E_%+.3lf_B_%+.3lf", vgtname, VGmin, realE, Bfield);

		
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
	
	    FILE *output, *fulloutput, *tsoutput;
	    char filename[300], filename3[300], filename_temp[300], fullfilename[300], tsfilename[300];
	    char checkname[300], direcname[300], conffile[350], strucfile[300], lstrucfile[350], disorderfile[300], pstrucfile[350];
	    char bandname1[300], bandname2[300], bandname3[300], mapname[400], maindirec[100], gsname[300];
            char command[600], inputlist[300];
	    char altconffile[100];
	    int usealtconf=0;
	    FILE *bandfile;
	    FILE *mapfile;
	    
	//Create directory and filenaming convention
	    sprintf(maindirec, "..");
	    
	     //check for command line arguments which vary these
		    for(i=1; i<argc-1; i++)
		    {
		      if(strcmp("-maindirec", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s", maindirec);
		      }
		      
		      //specifies that we should use an existing configuration
		      //useful for multiple runs on a single disorder configuration.
		      //must be in the same folder etc beforehand
		      //(best to use the "input" file from the initial run, together with this flag, and changes to runtime loop)
		      if(strcmp("-altconf", argv[i]) == 0)
		      {
			  sscanf(argv[i+1], "%s", altconffile);
			  usealtconf=1;
		      }
		      
		    }
	    
	    
	    
	    sprintf(direcname, "%s/res/%s_%s_%.0e/%s%d_%s", maindirec, systemtype, peritype, eta, geotype, length1, sysinfo);
	    
	    if(potdis==1)
	    {
	      sprintf(direcname, "%s/%s", direcname, disinfo);
	    }
	     if(cedge==1)
	    {
	      sprintf(direcname, "%s/%s", direcname, cedgeinfo);
	    }
	    if(dedge==1)
	    {
	      sprintf(direcname, "%s/%s", direcname, dedgeinfo);
	    }
	    if(eedge==1)
	    {
	      sprintf(direcname, "%s/%s", direcname, eedgeinfo);
	    }
	    if(abs_pots == 1)
            {
                sprintf(direcname, "%s/%s", direcname, "abs_pot");
            }
            if(abs_pots == 2)
            {
                sprintf(direcname, "%s/%s", direcname, "abs_pot2");
            }
	    
	    
	    sprintf(command, "mkdir -p %s", direcname);
	    system(command);
	    printf("# directory: %s\n", direcname);
	    
	    
	    sprintf(filename_temp, "%s_%s.conf%02d", loopinfo, job_name, conf_num); 

	    sprintf(strucfile, "%s/%s.struct", direcname, filename_temp);
	    sprintf(inputlist, "%s/%s.inputs", direcname, filename_temp);

	    sprintf(disorderfile, "%s/%s.disprof", direcname, filename_temp);

	    sprintf(filename, "%s/.%s.part%02d", direcname, filename_temp, this_proc);
		sprintf(fullfilename, "%s/.%s.full.part%02d", direcname, filename_temp, this_proc);
		sprintf(tsfilename, "%s/.%s.ts.part%02d", direcname, filename_temp, this_proc);
	    
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
	    
	FILE *listinputs;
	
	
	//Create main output file
	    output =fopen(filename, "w");
	    fclose(output);
	    
	    if (ishallbar==3 || ishallbar==4)
	    {
		    fulloutput = fopen(fullfilename, "w");
		    fclose(fulloutput);
		    tsoutput = fopen(tsfilename, "w");
		    fclose(tsoutput);
		    
	    }

	    
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
            System.patched=0;
            
            //printf("HERE! %d\n", Nrem);
	    
	    RectRedux *LeadCells[num_leads];

	    //timing
	    clock_t time;
	    time = clock();
	    
	    int can_start_yet=-1;
	    

	    

		if (usealtconf==0)
		{
			//generate, or read in, disorder configuration
			sprintf(conffile, "%s/.%s", direcname, filename_temp);
			//printf("%s\n", conffile);
		}
		if (usealtconf==1)
		{
			//generate, or read in, disorder configuration
			sprintf(conffile, "%s/.%s", direcname, altconffile);
			//printf("%s\n", conffile);

		}
			
			
		if(this_proc==0)
		{
			
			if(usealtconf==0)
			{
				(SysFunction) ( &System, SysPara, output_type, strucfile);

				//additional, general disorder routine(s) here.
				if(potdis == 1)
				{
					potentialDisorder (&System, &dispara, pdmap, disorderfile );
				}
				
				
				//additional edge related potentials
				if(cedge == 1)
				{
					customEdgePots (&System, &cedgepara);
				}
				if(dedge == 1)
				{
					customEdgePots (&System, &dedgepara);
				}
				if(eedge == 1)
				{
					customEdgePots (&System, &eedgepara);
				}
				
				//export the disorder configuration for the other processes calculating the same configuration
				if(procs>0)
				{
					exportRectConf(&System, conffile);
				}
			}
			if(usealtconf==1)
			{
				importRectConf(&System, length1, length2, conffile);
			}

				check = fopen(checkname, "w");
				fprintf(check, "%d", 0);
				fclose(check);
				
				
				listinputs=fopen(inputlist, "w");
				for(i=1; i<argc; i++)
				{
					fprintf(listinputs, "%s ", argv[i]);
				}
				fprintf(listinputs, "\n");
				fclose(listinputs);
		}




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
	    
            
                        
// 		BASE SYSTEM AND POTENTIALS FULLY GENERATED AND EXPORTED
//              FROM NOW ON LEADS, PATCHES ETC CAN BE ADDED
//              EVERYTHING FROM HERE ON RUNS EVERY TIME.... (is not loaded from config files)
            
            
            //double *cap_pots = NULL;

            if(abs_pots == 1)
            {
                cap_potential(&System, abs_pots_width);
            }
            if(abs_pots == 2)
            {
                cap_potential2(&System, abs_pots_width);
            }
            else
            {
                System.cap_pots = NULL;
            }
            
            

            



	  	  int halloutput;
		  
	    
		  
	  if(magsetup==0)
	  {
	    gauge = 0;
	  }
	  if(magsetup==1)
	  {
	    gauge = 1;
	  }
	
	
//LEAD AND PATCH GENERATION

            char PGFlib[100];
            cellDivision cellinfo; // definition moved to here so that 'group' sites can be set in Patchify
            cellinfo.group_dim = 0;

            int patchcenter=0, pct=0, pct2=0;     // tries to align patches by their centre (use similar odd/even lengths as device for best results)
            patch_para patchp = {};
            patchp.numpatches = numpatches;
            char PGFdir[100];
            int patchv=0; //auto position patchs along y (ac direction), rather than x
            
            sprintf(PGFdir, "%s/PGF/graphene/", maindirec);
            
            
            //PATCH SETTINGS -default is no additional patches, and main System region patched
            if (ispatched != 0)
            {
                if (ishallbar != 4 )
                {
                    printf("Patched mode only works for custom leads mode (hallbar = 4)");
                    exit(1);
                }
                if (nngm != 1 )
                {
                    printf("Patched mode only works for NNTB model");
                    exit(1);
                }
                
                 if (geo != 0 )
                {
                    printf("Patched mode only works for ZZ-edged devices at present!");
                    exit(1);
                }
                
                patchp.pcoords = createNonSquareIntMatrix(numpatches, 2);
                patchp.pl1 = createIntArray(numpatches);
                patchp.pl2 = createIntArray(numpatches);
                
                
                
//sscanf(argv[i+1], "%d", &((rib_lead_para *)(mleadps[j]->indiv_lead_para))->start_coord);
                
                
                for(i=1; i<argc-1; i++)
                {
                    if(strcmp("-usePGFlib" , argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%d", &patchp.usePGFlib);
                    }
                    
                    patchp.PGFlibloc = PGFdir;
                    if(strcmp("-PGFlibloc" , argv[i]) == 0)
                    {
                        sscanf(argv[i+1], "%s", patchp.PGFlibloc);
                    }
                }
                
                
                for (j=0; j<numpatches; j++)
                {
                     for(i=1; i<argc-1; i++)
                    {
                        sprintf(temp_in_string, "-patchv");
                        if(strcmp(temp_in_string, argv[i]) == 0)
                        {
                            sscanf(argv[i+1], "%d", &patchv);
                        }
                    }
                    
                    
                    patchp.pl1[j] = 10;
                    patchp.pl2[j] = 10;
                    if(patchv == 0)
                    {
                        patchp.pcoords[j][0] = - (length1+length2)*(j+1);
                        patchp.pcoords[j][1] = (length1+length2)*(j+1);
                    }
                    if(patchv == 1)
                    {
                        patchp.pcoords[j][0] =  (int) ((length1+length2)*(j+1) /sqrt(3));
                        patchp.pcoords[j][1] =  (int) ((length1+length2)*(j+1) /sqrt(3));
                    }
                    if(patchv == 2)
                    {
                        patchp.pcoords[j][0] =  -(int) ((length1+length2)*(j+1) /sqrt(3)) -length1/2;
                        patchp.pcoords[j][1] =  -(int) ((length1+length2)*(j+1) /sqrt(3)) -length1/2;
                    }
                    
                    
                    for(i=1; i<argc-1; i++)
                    {
                    
                        sprintf(temp_in_string, "-patch%da1", j);
                        if(strcmp(temp_in_string, argv[i]) == 0)
                        {
                            sscanf(argv[i+1], "%d", &patchp.pcoords[j][0]);
                        }
                        sprintf(temp_in_string, "-patch%da2", j);
                        if(strcmp(temp_in_string, argv[i]) == 0)
                        {
                            sscanf(argv[i+1], "%d", &patchp.pcoords[j][1]);
                        }
                        sprintf(temp_in_string, "-patch%dl1", j);
                        if(strcmp(temp_in_string, argv[i]) == 0)
                        {
                            sscanf(argv[i+1], "%d", &patchp.pl1[j]);
                        }
                        sprintf(temp_in_string, "-patch%dl2", j);
                        if(strcmp(temp_in_string, argv[i]) == 0)
                        {
                            sscanf(argv[i+1], "%d", &patchp.pl2[j]);
                        }
                        
                        if(strcmp("-usePGFlib" , argv[i]) == 0)
                        {
                            sscanf(argv[i+1], "%d", &patchp.usePGFlib);
                        }
                        
                        patchp.PGFlibloc = PGFdir;
                        if(strcmp("-PGFlibloc" , argv[i]) == 0)
                        {
                            sscanf(argv[i+1], "%s", patchp.PGFlibloc);
                        }
                        
                        if(strcmp("-patchcenter" , argv[i]) == 0 || strcmp("-patchcentre" , argv[i]) == 0)
                        {
                            sscanf(argv[i+1], "%d", &patchcenter);
                            
                            
                           
                        }
                        
                        
                    }
                    if(patchcenter == 1)
                    {
                        if(patchv == 0)
                        {
                            pct = (int)  ((length1 - patchp.pl1[j])/4);
                            pct2 = (int)  ((length1 - patchp.pl1[j])/2 -  (length1 - patchp.pl1[j])/4) ;
                            
                            patchp.pcoords[j][0] -=  pct;
                            patchp.pcoords[j][1] -= pct2;
                        }
                        
                        if(patchv == 1 || patchv == 2)
                        {
                            pct = (int)  ((length2 - patchp.pl2[j])/2);
                            pct2 = (int)  ((length2 - patchp.pl2[j])/2) ;
                            
                            patchp.pcoords[j][0] +=  pct;
                            patchp.pcoords[j][1] -= pct2;
                        }
                        
//                          printf("# %d, %d\n", patchp.pcoords[j][0], patchp.pcoords[j][1]);
                        
                    }
                    
                   
                    
                    
                }
                
                 
                    if(patchp.usePGFlib == 1)
                    {
                        sprintf(command, "mkdir -p %s", PGFdir);
                        system(command);
                        printf("# PGF directory: %s\n", patchp.PGFlibloc);
                        
                        
                        
                    }
                
                if(output_type == 1)
                {
                        sprintf(pstrucfile, "%s.patch", strucfile);
                }
                    
                //call Patchify...
                Patchify (&System, &patchp, &cellinfo, output_type, pstrucfile);
                System.patched=1;
                System.patch_params = &patchp;
                pos = System.pos;
                //exit(1);

                //filename -- change leadconf from CUSTOM to PGF (this can be overridden by cmd line using -leadconfname

                
                
            }
	  



		  
	  if(ishallbar == 1 || ishallbar == 2)
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
		
                if(ismagnetic != 0)
                {
                    if(magsetup==0)
                    {
                    gauge = 3;
                    }
                    if(magsetup==1)
                    {
                    gauge = 2;
                    }
                    if(magsetup==2)
                    {
                    gauge = 0;
                    }
                }
		
		
		pos = System.pos;
		
	  }

	  
	  
	  
	    if(ishallbar ==0 && bandsonly == 0)
	    {
	      genLeads(&System, LeadCells, num_leads, gen_leads_mode, &leadp);
	    }
	  
		//adds (averaged) sublattice dependent potentials to left and right leads (leads 0 and 1)
		if(xleadsavgpot==1 && strcmp("SUBLATTICEPOT", systemtype) == 0 && ishallbar != 4)
		{
			genSublatticeLeadPots(LeadCells, &sublp);
		}
	  
                //similar for strain
                if(leadstrain==1 && strcmp("SYMSTRAIN", systemtype) == 0 && ishallbar != 4)
		{
			customLeadStrain(LeadCells, &symsp);
		}
		
		
	  
	  
	  if(ishallbar == 4)
	  {
		leadp.multiple = mleadps;
		leadp.leadsfn = &multipleCustomLeads;
		genCustomLeads (&System, LeadCells, num_leads, &leadp);
		
	  }
	  
	  
	  //lead edge potentials! 
		//this is only implemented for ribbon and custom leads (ishallbar == 0 || 4)
		if(cedge==1 && cedge_leads == 1)
		{
			for(i=0; i<num_leads; i++)
			{
				if(ishallbar==0)
					customEdgePots (LeadCells[i], &cedgepara);
				
				if(ishallbar==4)
				{
					if(strcmp("RIBBON", (mleadps[i]->name)) == 0 && ( (mleadps[i]->def_pos) == 0 || (mleadps[i]->def_pos) == 1 ) )
					{
						customEdgePots (LeadCells[i], &cedgepara);
					}
				}
				
			}				
		}
		
		if(dedge==1 && dedge_leads == 1)
		{
			for(i=0; i<num_leads; i++)
			{
				if(ishallbar==0)
					customEdgePots (LeadCells[i], &dedgepara);
				
				if(ishallbar==4)
				{
					if(strcmp("RIBBON", (mleadps[i]->name)) == 0 && ( (mleadps[i]->def_pos) == 0 || (mleadps[i]->def_pos) == 1 ) )
					{
						customEdgePots (LeadCells[i], &dedgepara);
					}
				}
				
			}				
		}
		
		if(eedge==1 && eedge_leads == 1)
		{
			for(i=0; i<num_leads; i++)
			{
				if(ishallbar==0)
					customEdgePots (LeadCells[i], &eedgepara);
				
				if(ishallbar==4)
				{
					if(strcmp("RIBBON", (mleadps[i]->name)) == 0 && ( (mleadps[i]->def_pos) == 0 || (mleadps[i]->def_pos) == 1 ) )
					{
						customEdgePots (LeadCells[i], &eedgepara);
					}
				}
				
			}				
		}
		
		
// 		if(abs_pots == 1)
//                 {
//                     cap_potential(&System, abs_pots_width);
//                 }
//                 else
//                 {
//                     System.cap_pots = NULL;
//                 }
		
		
        //lead edge absorbing potentials
		
                for(i=0; i<num_leads; i++)
                {
                    if(abs_pots == 1)
                    {
                        if(ishallbar==0)
                        {
                                cap_potential(LeadCells[i], abs_pots_width);
                        }
                        
                        if(ishallbar==4)
                        {
                                if(strcmp("RIBBON", (mleadps[i]->name)) == 0 && ( (mleadps[i]->def_pos) == 0 || (mleadps[i]->def_pos) == 1 ) )
                                {
                                        cap_potential(LeadCells[i], abs_pots_width);
                                }
                        }
                    }
                    else
                    {
                        if( bandsonly == 0)
                            (LeadCells[i]->cap_pots) = NULL;
                    }
                        
                }				
		
		
	  
//Device Connectivity

	  cnxProfile cnxp;
	 
 	    cnxp.max_neigh=max_neigh;
	    device_connectivity (&System, connectrules, default_connection_params, &cnxp);

		  time = clock() - time;
		  printf("#made connection profile in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		  time = clock();
   	   //printConnectivity (&System, &cnxp);

		
	  gen_start_params start_p ={};
	 
	    start_p.rule = defdouble_connect_rule;
	    start_p.rule_params = default_connection_params;
	    start_p.num_leads = num_leads;
	    start_p.Leads = LeadCells;
	  
	  void *starting_ps;
          starting_ps = &start_p;
          
          
          
	if(metal_leads ==1)
              starting_ps = &mxsp;
	  
	if(ishallbar == 4)
	{	
		starting_ps = &cstart_p;
		cstart_p.rule = defdouble_connect_rule;
		cstart_p.rule_params = default_connection_params;
		cstart_p.num_leads = num_leads;
		cstart_p.Leads = LeadCells;
		starting_cell_mode = 5;
	}
		
	 // genStartingCell(&System, &cellinfo, 2, NULL);
 
	 if(bandsonly == 0)
		genStartingCell(&System, &cellinfo, starting_cell_mode, starting_ps);
	 
 	//exit(1);


		if(output_type == 1 && bandsonly == 0)
		{
			sprintf(lstrucfile, "%s.leads", strucfile);
			printOutLeadStrucs(&System, LeadCells, &cellinfo, lstrucfile);
			
		}
		  
		  
	  
		time = clock() - time;
		printf("#generated starting cell in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();
		
	if(bandsonly == 0)
	  cellSplitter(&System, &cnxp, &cellinfo);
	  
		time = clock() - time;
		printf("#split cells in %f seconds\n", ((float)time)/CLOCKS_PER_SEC);
		time = clock();

        //exit(1);

		
		
	  //hopping parameters and gauge info

	//LOAD HOPPING PARAMS

		hoppara.num_neigh = num_neigh;
		hoppara.hops=createCompArray(hoppara.num_neigh);
		hoppara.NN_lowdis=createDoubleArray(hoppara.num_neigh);
		hoppara.NN_highdis=createDoubleArray(hoppara.num_neigh);
		hoppara.NN_shifts=createDoubleArray(hoppara.num_neigh);
		hoppara.NN_zmin=createDoubleArray(hoppara.num_neigh);
		hoppara.NN_zmax=createDoubleArray(hoppara.num_neigh);
		
		
		
		//read a hopping profile
		//taken to nngm-th order
		//nngm must be less than the maximum for this to work
 		if(num_neigh<=(hop_to_load->num_neigh))
		{
// 			hoppara.max_neigh = (hop_to_load->max_neigh)[nngm];
			for(i=0; i<num_neigh; i++)
			{
				hoppara.hops[i]=(hop_to_load->hops)[i];	
				hoppara.NN_lowdis[i] = (hop_to_load->NN_lowdis)[i];
				hoppara.NN_highdis[i] = (hop_to_load->NN_highdis)[i];
				hoppara.NN_shifts[i] = (hop_to_load->NN_shifts)[i];
				hoppara.NN_zmin[i] = (hop_to_load->NN_zmin)[i];
				hoppara.NN_zmax[i] = (hop_to_load->NN_zmax)[i];
				
			}
		}
		else
			exit(1);
// 	  
// 	  hoppara.hops[0]=t0;
// 	  hoppara.NN_lowdis[0] = NNlowdis;
// 	  hoppara.NN_highdis[0] = NNhighdis;
// 
// 	  
// 	  if(nngm>1)
// 	  {
// 	    hoppara.hops[1]=t1;
// 	    hoppara.NN_lowdis[1]=NNlowdis1;
// 	    hoppara.NN_highdis[1]=NNhighdis1;
// 	    hoppara.NN_shifts[1]= onsitec1;
// 	  }
// 	    
// 	  if(nngm>2)
// 	  {
// 	    hoppara.hops[2]=t2;
// 	    hoppara.NN_lowdis[2]=NNlowdis2;
// 	    hoppara.NN_highdis[2]=NNhighdis2;
// 	    hoppara.NN_shifts[2]= onsitec2;
// 	    
// 	  }
// 	  hoppara.t0=t0;
		
// 		hoppara.num_neigh = hop_to_load->num_neigh;
		
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
	  
	  //this is the special case of no field in left and right leads for a hall bar
	  //field is turned on at the end of the buffer region
	  if(magsetup == 2 && ishallbar == 1 && gauge == 0)
	  {
	    res[0] = 1;
	    reslimits[0][0] = pos[2*length1*(buffer_rows)][0];
	    reslimits[0][1] = pos[2*length1*(length2-buffer_rows-1)][0];
	    
	  }
	  
	  if(gauge == 3)
	  {

	    
	    reslimits[1][0] = pos[2*length1*length2 + 2*(hallp.toppc[0]*hallp.toppw[0] +  hallp.toppc[1]*hallp.toppw[1] + hallp.botpc[0]*hallp.botpw[0])-1 ][1];
	    reslimits[1][3] = pos[2*length1*length2 + 2*hallp.toppc[0]*hallp.toppw[0]  -1][1];

	    

	    
	  }

	
	  hoppara.restrics = res;
	  hoppara.limits = reslimits;
	  
	 

	  //these are the default values for most cases, i.e. leads are continuation of system
          //different functions are needed for different systems.
            leadp.hopfn = hopfn;
            leadp.hoppara = &hoppara;
            leadp.leadsfn = &multipleLeads;
            
            if(metal_leads == 1)
            {
                leadp.hoppara = &metal_hop_p;
                leadp.leadsfn = &multipleSimplestMetalLeads;
            }
            
            if(ishallbar == 4)
	    {
		    leadp.leadsfn = &multipleCustomLeads;
	    }

	  
	  trans_params tpara = {};
	  tpara.num_leads=num_leads;
	  tpara.TRsym=ismagnetic;
          tpara.filename=gsname;
	  double **transmissions = createNonSquareDoubleMatrix(num_leads, num_leads);
	  double **ttildas=createNonSquareDoubleMatrix(num_leads, num_leads);
	  double **mtrans = createNonSquareDoubleMatrix(num_leads-1, num_leads-1);
	  double *vecx = createDoubleArray(num_leads-1);
	  double *vecb = createDoubleArray(num_leads-1);
	  

	  tpara.transmissions = transmissions;
          tpara.ispatched = ispatched;
          if(ispatched != 0)
          {
              tpara.patchpara = &patchp;
          }
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
	 if(strcmp("VG", loop_type) == 0 )
	 {
	   engdeppots = createDoubleArray(Ntot);
	   origpots = createDoubleArray(Ntot);
	   for(i=0; i<Ntot; i++)
	   {
	     origpots[i] = System.site_pots[i];
	   }
	   
	   if( bandsonly == 0)
	   {
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
	 }

	 
	 //BAND STRUCTURES
	  
	 if(makebands == 1)
	 {
	   bands = createDoubleArray(Nrem);
	   weights = createNonSquareDoubleMatrix(Nrem, length2);
	   projections = createNonSquareDoubleMatrix(Nrem, Nrem);
	   bandmode=0;
	   kxstep = (kxmax-kxmin)/(kxpoints -1);
	   if(kxpoints ==1)
		   kxstep=100.0;
	   
	   
	    if(ismagnetic == 1)
	    {
	      hoppara.Btes=Bmin;
	    }
	    		double *engdeppots2;

	    //generate engdeppots if required
	      if(strcmp("VG", loop_type) == 0)
	      {
		engdeppots2 = createDoubleArray(Ntot);
		gate_induced_pot (vgtype, &System, engdeppots2, VGmin, edge_cut, subs_thick, subs_epsr);
	      }

	    
	    
	    //set total onsite potentials
	    if(strcmp("VG", loop_type) == 0)
	    {
	      for(i=0; i<Ntot; i++)
	      {
		System.site_pots[i] = origpots[i] + engdeppots2[i];
	      }
	    }
	    
	    
	      for(kxl=kxmin; kxl< kxmax; kxl += kxstep)
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
	    
	    free(engdeppots2);
	 }
	  
	  
	  
	  
	    int do_linear_decomp=1;
	  
            
          // MAIN LOOP!  
	  
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
		    if(ishallbar==0)
		    {
			for(j=0; j<num_leads; j++)
			{
				gate_induced_pot (vgtype, LeadCells[j], lead_engdeppots[j], VG, edge_cut, subs_thick, subs_epsr);
			}
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
		  mapnow = mapmode;
		}
		  
	      }
	      
	      
	      if(mapnow > 0)
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
	      sprintf(gsname, "%s/E_%+.4lf_B_%+.3lf_%s.conf%02d", direcname, realE, Bfield, job_name, conf_num);
	     
	    
	      hoppara.Btes=Bfield;
	      
	      //needs generalising for multiple transmission cases
	      kavg=0.0;
	      
	      for(k=0; k<kpts; k++)
	      {
		hoppara.kpar = kmin + k*kstep;

		genTransmissions(realE+eta*I, &System, LeadCells, &cnxp, &cellinfo, hopfn, &hoppara, &leadp, &tpara, mapnow, ldoses, currents);
		kavg += (transmissions[0][1]/kpts);
	      }

// 	      printDMatrix(transmissions, 4);

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
	      double probe_currs[num_leads];

	      
              
	      if(mapnow > 0)
	      {
		     sprintf(mapname, "%s/E_%+.4lf_B_%+.3lf_%s.conf%02d.ldos", direcname, realE, Bfield, job_name, conf_num);
		     
			if(strcmp("VG", loop_type) == 0)
			{
				sprintf(mapname, "%s/VG_%s_%.2lf_E_%+.2lf_B_%+.3lf_%s.conf%02d.ldos", direcname, vgtname, VG, realE, Bfield, job_name, conf_num);
			}

		      mapfile = fopen(mapname, "w");
		      
		      for(i=0; i<Ntot; i++)
		      {
			fprintf(mapfile, "%lf	%lf	%e\n", (pos)[i][0], (pos)[i][1], ldoses[i]);
		      }
		      fclose(mapfile);
		      

	      
		      for(j=0; j<num_leads; j++)
		      {
			
			sprintf(mapname, "%s/E_%+.4lf_B_%+.3lf_%s.conf%02d.cmaps_l%d", direcname, realE, Bfield, job_name, conf_num, j);

			if(strcmp("VG", loop_type) == 0)
			{
				sprintf(mapname, "%s/VG_%s_%.2lf_E_%+.2lf_B_%+.3lf_%s.conf%02d.cmaps_l%d", direcname, vgtname, VG, realE, Bfield, job_name, conf_num, j);
			}

			  mapfile = fopen(mapname, "w");
			  for(i=0; i<Ntot; i++)
			  {
			    fprintf(mapfile, "%lf	%lf	%e %e\n", (pos)[i][0], (pos)[i][1], currents[j][i][0], currents[j][i][1] );
			  }
			  fclose(mapfile);
			  
		      }
		      
		      
		  
			   

		      
		    
	      }
	      

	      
	      int iprime, jprime; 
	      
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
//                                                  printDMatrix(transmissions, 6);
//                                                 printDMatrix(ttildas, 6);
//                                                 printDMatrix(mtrans, 5);

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
                            if(mapnow > 0)
                            {
                                probe_pots[0] = 1.0;
                                probe_pots[1] = 0.0;
                                probe_pots[2] = vecx[1];
                                probe_pots[3] = vecx[2];
                                probe_pots[4] = vecx[3];
                                probe_pots[5] = vecx[4];
                                
        // 			printf("#curr: %.10e\n", vecx[0]);

                                sprintf(mapname, "%s/E_%+.4lf_B_%+.3lf_%s.conf%02d.cmaps_multi", direcname, realE, Bfield, job_name, conf_num);
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
                                
                     
                    
                    //spin hall bar transmissions and resistances
                    if (ishallbar==2)
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
                        
                        mtrans[2][2] = 1 / ttildas[2][2];
                        
                        for(i=0; i<num_leads-1; i++)
                        {
                            if(i!=2)
                            {
                                iprime = i;
                                if(i == 4)
                                {
                                    iprime = 5;
                                }
                                mtrans[2][i] = - ttildas[2][iprime] / ttildas[2][2];
                                mtrans[i][2] =  ttildas[iprime][2] / ttildas[2][2];
                                
                                for(j=0; j< num_leads-1; j++)
                                {
                                    if(j!=2)
                                    {
                                        jprime = j;
                                        if(j == 4)
                                        {
                                            jprime = 5;
                                        }
                                        mtrans[i][j] = ttildas[iprime][jprime] - (ttildas[iprime][2] * ttildas[2][jprime])/ttildas[2][2] ;
                                    }
                                }
                            }
                        }
                        vecb[0] = 0.0; vecb[1] = 0.0;  vecb[2] = 1.0; vecb[3] = 0.0; vecb[4] = 0.0; 
                        vecx[0] = 0.0; vecx[1] = 0.0;  vecx[2] = 0.0; vecx[3] = 0.0; vecx[4] = 0.0; 
                        
//                         printDMatrix(mtrans, 5);
                        LinEqnDouble (mtrans, vecb, vecx, num_leads-1);
                        
                        
                        
                        if(strcmp("E", loop_type) == 0)
                        {
                            fprintf(output, "%lf	%.12e\n", realE, (vecx[3]-vecx[4])/vecx[2]);		

                        }
                        if(strcmp("B", loop_type) == 0)
                        {
                            fprintf(output, "%lf	%.12e\n", Bfield, (vecx[3]-vecx[4])/vecx[2]);		

                        }
                        
                        
                        
                            //make a composite map for multiprobe systems
                            if(mapnow > 0)
                            {
                                probe_pots[0] = vecx[0];
                                probe_pots[1] = vecx[1];
                                probe_pots[2] = 1.0;
                                probe_pots[3] = vecx[3];
                                probe_pots[4] = 0.0;
                                probe_pots[5] = vecx[4];
                                

                                sprintf(mapname, "%s/E_%+.4lf_B_%+.3lf_%s.conf%02d.cmaps_multi", direcname, realE, Bfield, job_name, conf_num);
                                mapfile = fopen(mapname, "w");
                                for(i=0; i<Ntot; i++)
                                { 
                                xc=0.0; yc=0.0;
                                k=4;

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
                                    
                                }
                                
                                
                        
                        
                                
                                    fprintf(mapfile, "%lf	%lf	%.12e %.12e\n", (pos)[i][0], (pos)[i][1], xc, yc);
                                }
                                fclose(mapfile);
                            }
                        
                        
                        
                        
                        
                    }
                        

                     
                    //non local resistances in general multi probe measurements
                    //currIn and currOut define the probes current is driven between
                    //probes are reordered 0-N-1 (above) so that 0=currIn, N-2=currOut, 1 .. N-1 define the nonlocal potential difference
                    if (ishallbar==3 || ishallbar==4)
                    {	
			//fulloutput prints potentials and currents of each lead
			fulloutput =fopen(fullfilename, "a");
			
			//tsoutput prints all elements of the transmission matrix
			tsoutput =fopen(tsfilename, "a");
			
//                         if(num_leads==4)
// 			{
				EmptyDoubleMatrix(ttildas, num_leads, num_leads);
				EmptyDoubleMatrix(mtrans, num_leads-1, num_leads-1);
				
				for(i=0; i<num_leads; i++)
				{
					for(j=0; j<num_leads; j++)
					{
						if(i!=j)
						{
							ttildas[i][i] += transmissions[leadOrder[j]][leadOrder[i]] ; 
							ttildas[i][j] = -transmissions[leadOrder[i]][leadOrder[j]] ;
						}
					}
				}
				
				mtrans[0][0] = 1 / ttildas[0][0];
				vecb[0] = 1.0;
				
				for(i=1; i<num_leads-1; i++)
				{
					mtrans[0][i] = - ttildas[0][i] / ttildas[0][0];
					mtrans[i][0] =  ttildas[i][0] / ttildas[0][0];
					
					for(j=1; j< num_leads-1; j++)
					{
						mtrans[i][j] = ttildas[i][j] - (ttildas[i][0] * ttildas[0][j])/ttildas[0][0] ;
					}
					
					vecb[i] = 0.0;
					vecx[i] = 0.0;
					
					
				}

				
// 				for(i=0; i<num_leads-1; i++)
// 				{
// 					for(j=0; j<num_leads-1; j++)
// 					{
// 						if(mtrans[i][j] < 1.0E-10)
// 							mtrans[i][j] = 0.0;
// 					}
// 				}
				
			
				LinEqnDouble (mtrans, vecb, vecx, num_leads-1);
				
				probe_pots[leadOrder[0]] = 1.0;
				probe_pots[leadOrder[num_leads-1]] = 0.0;
				probe_currs[leadOrder[0]] = vecx[0];
				probe_currs[leadOrder[num_leads-1]] = -vecx[0];
				
				for(i=1; i<num_leads-1; i++)
				{
					probe_pots[leadOrder[i]] = vecx[i];
					probe_currs[leadOrder[i]] = 0.0;
				}
					
				
				if(strcmp("E", loop_type) == 0)
				{
					fprintf(output, "%lf	%.12e\n", realE, (vecx[2]-vecx[1])/vecx[0]);	
					fprintf(fulloutput, "%lf\t", realE);
					fprintf(tsoutput, "%lf\t", realE);

					for(i=0; i<num_leads; i++)
					{
						fprintf(fulloutput, "%.12e\t", probe_pots[i]);
					}
					for(i=0; i<num_leads; i++)
					{
						fprintf(fulloutput, "%.12e\t", probe_currs[i]);
					}
					fprintf(fulloutput, "\n");
					
					for(i=0; i<num_leads; i++)
					{
						for(j=0; j<num_leads; j++)
						{	
							fprintf(tsoutput, "%.12e\t", transmissions[i][j]);
						}
					}
					fprintf(tsoutput, "\n");
					
				}
				if(strcmp("B", loop_type) == 0)
				{
					fprintf(output, "%lf	%.12e\n", Bfield, (vecx[2]-vecx[1])/vecx[0]);	
					fprintf(fulloutput, "%lf\t", Bfield);
					fprintf(tsoutput, "%lf\t", Bfield);
					for(i=0; i<num_leads; i++)
					{
						fprintf(fulloutput, "%.12e\t", probe_pots[i]);
					}
					for(i=0; i<num_leads; i++)
					{
						fprintf(fulloutput, "%.12e\t", probe_currs[i]);
					}
					fprintf(fulloutput, "\n");
					
					for(i=0; i<num_leads; i++)
					{
						for(j=0; j<num_leads; j++)
						{	
							fprintf(tsoutput, "%.12e\t", transmissions[i][j]);
						}
					}
					fprintf(tsoutput, "\n");
					
				}
				
				if(strcmp("VG", loop_type) == 0)
				{
					fprintf(output, "%lf	%.12e\n", VG, (vecx[2]-vecx[1])/vecx[0]);	
					fprintf(fulloutput, "%lf\t", VG);
					fprintf(tsoutput, "%lf\t", VG);
					for(i=0; i<num_leads; i++)
					{
						fprintf(fulloutput, "%.12e\t", probe_pots[i]);
					}
					for(i=0; i<num_leads; i++)
					{
						fprintf(fulloutput, "%.12e\t", probe_currs[i]);
					}
					fprintf(fulloutput, "\n");
					
					for(i=0; i<num_leads; i++)
					{
						for(j=0; j<num_leads; j++)
						{	
							fprintf(tsoutput, "%.12e\t", transmissions[i][j]);
						}
					}
					fprintf(tsoutput, "\n");
					
				}
				
				
				
				
				//make a composite map for multiprobe systems
				if(mapnow > 0)
				{
					sprintf(mapname, "%s/E_%+.4lf_B_%+.3lf_%s.conf%02d.cmaps_multi", direcname, realE, Bfield, job_name, conf_num);
					mapfile = fopen(mapname, "w");
					for(i=0; i<Ntot; i++)
					{ 
						xc=0.0; yc=0.0;
						k=leadOrder[num_leads-1];

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
						
						}

						fprintf(mapfile, "%lf	%lf	%.12e %.12e\n", (pos)[i][0], (pos)[i][1], xc, yc);
					}
					fclose(mapfile);
				}
				
				
				
				
// 			}
                    }
                    
	      
	      fclose(output);

	      if (ishallbar==3 || ishallbar==4)
		{
			fclose(fulloutput);
			fclose(tsoutput);
		}
	      
                time = clock() - time;
		printf("#iteration %d complete in %f seconds\n", en, ((float)time)/CLOCKS_PER_SEC);
		time = clock();
	      
	      
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
		    
		    if(ishallbar==4 || ishallbar==3)
		    {
			sprintf(command, "cat %s.full.part* | sort -n > %s/%s.full.dat", filename, direcname, filename_temp);
			system(command);
			
			sprintf(command, "cat %s.ts.part* | sort -n > %s/%s.ts.dat", filename, direcname, filename_temp);
			system(command);
		    }
		    
		    
		    
		    if(conf_num != 0)
		    {
		      sprintf(command, "rm %s/.%s.*", direcname, filename_temp);
		      system(command);
		    }
	
		  
		}
		
}





