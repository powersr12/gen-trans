
//General transport code
//Built from sections of ribbon-asymm code and others
//Accepts devices with RectRedux geometries, but calls specified routines for defining cells and setting hopping values
//This file is for developing routines - finished products should go in separate files to be included in compilation


//command line version
#include "test.h"

#define eta 1.0e-5



//remove unneeded

main(int argc, char *argv[])
{
	int i, j,k,l;
	
	double realE=0.0;
	double radius=2.0;
	int Ntot, Ntot2, Ntot3, Nrem, Nrem3, Ntemp, Ntemp2, length1=6, length1a=length1, length2 = length1*3 , geo=0;
	int botcells=0, topcells=0, midcells=1, procs=1, this_proc=0, remainder, Epts_temp=11, en, numfin, num_rows=10, Epoints=11, num_hole_rows;
 	double Emax=-0.9, Emin=0.9, Estep=0.04, Emin_temp=0.0, pos_err=0.0, rad_err=0.0, rad_dis=0.0, dis_pot=0.0, anglein, angle;
	double _Complex En;
	int conf_num; 
	char job_name[64];
	double suba_conc, suba_pot, subb_conc, subb_pot;
	int output_type=1;   //=0 no structure output, =1 just atoms by sublattice, =2 atoms and connectors
	int buffer_rows, mappings;		   
			      	    double dos, cond, cond2;

	
      if(argc ==17)
      {
	      sscanf(argv[1], "%s", job_name);
	      sscanf(argv[2], "%d", &geo);
	      sscanf(argv[3], "%d", &length1);
	      sscanf(argv[4], "%d", &length2);
	      sscanf(argv[5], "%d", &buffer_rows);
	      
	      sscanf(argv[6], "%lf", &suba_conc);
	      sscanf(argv[7], "%lf", &suba_pot);
	      
	      sscanf(argv[8], "%lf", &subb_conc);
	      sscanf(argv[9], "%lf", &subb_pot);
	      
	      sscanf(argv[10], "%lf", &Emin);
	      sscanf(argv[11], "%lf", &Emax);
	      sscanf(argv[12], "%d", &Epoints);
	      sscanf(argv[13], "%d", &mappings); 	//how frequently to make maps (0 = never, 1 = always, 2 = every second energy etc, -n every nth energy but only for conf 0)

	      sscanf(argv[14], "%d", &conf_num);
	      sscanf(argv[15], "%d", &procs);
	      sscanf(argv[16], "%d", &this_proc);
      }
      
      else
      {
	      exit(1);
      }
      

	
      FILE *output, *output2;
     char filename[160], filename2[160], filename3[160], filename_temp[128], checkname[128], finalname[128], direcname[128], conffile[128], strucfile[128];
      
      char command[256];
      sprintf(direcname, "res/%s", job_name);
      sprintf(command, "mkdir -p %s", direcname);
      system(command);
      


      sprintf(filename_temp, "%s.conf%02d", job_name, conf_num); 
      sprintf(strucfile, "%s/%s.struct", direcname, filename_temp);

      sprintf(filename, "%s/.%s.part%02d", direcname, filename_temp, this_proc);
      output =fopen(filename, "w");
      fclose(output);

      
      
      //status files for multiprocessor runs
      
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
      
      
      int number, number2;
      int rows =1;
      int numholes, numholes2, numholes_clean;
           
      double offset= sqrt(3) *length1a /4;

    
      
//------------------------------ START OF THE ACTUAL INTERESTING COMPUTATIONAL NON-I/O BIT -----------------------------------//
      
           	//generate System Geometry

      
	    RectRedux System = {};
	    System.geo=geo;
	    System.length=length1;
	    System.length2=length2;
	    System.Nrem=&Nrem;

	    //generate, or read in, disorder configuration
		sprintf(conffile, "%s/.%s", direcname, filename_temp);
		
		if(this_proc==0)
		{
		  genSublatticeDevice (&System, buffer_rows, suba_conc, suba_pot, subb_conc, subb_pot, conf_num, 1, strucfile);
	   
		  //export the disorder configuration for the other processes calculating the same configuration
		  if(procs>0)
		  {
		    exportRectConf(&System, conffile);
		  }
		}
		
		  
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


	
 	      

	      double *ldoses = createDoubleArray(*(System.Nrem));
	      double **conds = createNonSquareDoubleMatrix(*(System.Nrem), 3);

	      int *indices = createIntArray(2*System.length*System.length2);
	      int **neigh = createNonSquareIntMatrix(2*System.length*System.length2, 3);
	      
	      for(en=0; en < Epts_temp; en++)
	      {

		      realE =  Emin_temp + en*Estep;
		      sprintf(filename3, "%s/%s.en_%.3lf", direcname, filename_temp, realE);

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
		      
		      printf("%lf	%e\n", realE,  genConduc5(realE + eta*I, &System, -1.0));
			
		      //cond2 =genConduc5(realE + eta*I, &System, -1.0);
		      
			//    printf("%lf	%e	%e\n", realE, genConduc5(realE + eta*I, &System, -1.0),genConduc4(realE + eta*I, &System, ldoses, conds, indices, neigh, -1.0) );
			    
		      output =fopen(filename, "a");
		      fprintf(output, "%lf	%e\n", realE, cond2);
		      fclose(output);

		      //printf("%lf	%lf\n", realE, cond2);


 	      }
	      
	      
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














