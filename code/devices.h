#pragma once


//moved here for compilation reasons...
typedef struct {
	int num_cells;		//number of cells the system has been partitioned into
	int cell1dim;
	int *cells_site_order;  //ordering of the sites in new cell structure
	int *starting_index;	//index of first site in each cell within cells_site_order
	int *cell_dims;	//dimension of each cell in cells_site_order
	int *sites_by_cell;	//what cell each site is in
	
	int group_dim;		//a group has sites that should be in the same cell
	int group_cell;
	int *group_sites;		//i.e. due to self energy terms
	
	int num_leads;
	int *lead_dims;
	int *lead_sites;
}cellDivision;



typedef struct {
	int geo;		//=0 for zigzag, =1 for armchair
	int length;
	int length2;
	int **siteinfo;		//in lattice? which sublattice?
	double **pos;
        double **pert_pos;
	double *site_pots;
        double *cap_pots;
	int **chaininfo;
	int *Nrem;
	int *Ntot;
}RectRedux;

typedef struct {
        int numpatches;
        int *pNrem;
        int **pcoords;
        int *pl1;
        int *pl2;
        int usePGFlib;
        char *PGFlibloc;
}patch_para;

typedef struct {
	char *name;
	void *indiv_lead_fn;
	void *indiv_lead_para;
	void *indiv_gen_fn;
	int def_pos; //0->left, 1->right, 2->top, 3->bottom
}multiple_para;

typedef struct {
	double **shift_vecs;
	void *hopfn;
	void *leadsfn;
	void *hoppara;
	multiple_para **multiple;
	void *additional_params;
}lead_para;

typedef struct {
	double *shift_vec;
	void *hopfn;
	void *hoppara;
	int geo;
	int width;
	int start_coord;
	int def_pos; //0->left, 1->right, 2->top, 3->bottom
	double homo_A_pot;
	double homo_B_pot;
}rib_lead_para;

typedef struct {
	double *shift_vec;
	void *hopfn;
	void *hoppara;
	int geo;
	int width;
	int start_coord;
	int def_pos; //0->left, 1->right, 2->top, 3->bottom
	void *bilayer_para;
}blg_lead_para;

typedef struct {
	void *hopfn;
	void *hoppara;
	int geo;
	double width;	//x-dim
	double width2;	//y-dim
	double start_coord;
	double start_coord2;
	int def_pos; //0->left, 1->right, 2->top, 3->bottom, 4->fullwidth, 5-> absolute positioning
}metal_lead_para;
	
	
	
typedef struct {
	int buffer_rows;
	double a_conc;
	double a_pot;
	double b_conc;
	double b_pot;
	int seed;
}subl_para;

typedef struct {
	int cellswide;	//how many blocks of GB wide the AGNR is 
	int GBpos; 	//after what standard unit cell is the GB dimer included
	int anddis;	//include anderson disorder near GB?
	double andD;	//size of region containing anderson disorder
	double andW;	//strength of anderson disorder
	int vacdis;
	double vacW;
	double vacP;
	double ycut;
	int seed;
}simple558_para;

typedef struct {
	double AA_mass;
	double AA_pot;
	double AB_mass;
	double AB_pot;
	double BA_mass;
	double BA_pot;
	int lM;
	double x_offset;
	double y_offset;
	int seed;
}submoire_para;


typedef struct {
	int type_shift; //=0 AA. =1 AB(Bernal), =2(Custom, uses *shiftvec and  shift_angle (ANGLE NOT IMPLEMENTED YET!)
	double *shift_vec;
	double shift_angle;
	double zsep;
	double *subpots;
	int skews;
}bilayer_para;

typedef struct {
	int buffer_rows;
	double a_conc1;   
	double a_pot1;
	double b_conc1;
	double b_pot1;
	double a_conc2;
	double a_pot2;
	double b_conc2;
	double b_pot2;
	int xory;
	double int_pos;
	double int_width;
	int seed;
}subint_para;

typedef struct {
	int buffer_rows;
	double a_conc1;   
	double a_pot1;
	double b_conc1;
	double b_pot1;
	double a_conc2;
	double a_pot2;
	double b_conc2;
	double b_pot2;
	double a_conc3;
	double a_pot3;
	double b_conc3;
	double b_pot3;
	int xory;
	double int_pos1;
	double int_pos2;
	double int_width1;
	double int_width2;
	int seed;
}sub2int_para;

typedef struct {
	int buffer_rows;
	int AD_length;
	int AD_length2;
	char *latgeo;
	char *dotgeo;
	double AD_rad;
	double AD_rad2;
	int lat_width;
	int lat_length;
	int seed;
	int isperiodic;
	double radfluc;
	double xyfluc;
}adot_para;

typedef struct {
	int buffer_rows;
	int SD_length;
	int SD_length2;
	char *latgeo;
	char *dotgeo;
	double SD_rad;
	double SD_rad2;
	int lat_width;
	int lat_length;
	int seed;
	int isperiodic;
	double radfluc;
	double xyfluc;
        double a_conc;
	double a_pot;
	double b_conc;
	double b_pot;        
}sldot_para;

typedef struct {
	int buffer_rows;
	char *latgeo;
	char *bubgeo;
        int cell_length;
	int cell_length2;
	double bub_rad;
        double bub_height;
	int lat_width;
	int lat_length;
	int seed;
	int isperiodic;
	double radfluc;
	double xyfluc;
        double zfluc;
        double ifluc;
}bubble_para;


typedef struct {
	int buffer_rows;
	int sruns;
	double smax;
	double minper;
	int vacruns;
	double vacprob;
	int seed;
}edgedis_para;


typedef struct {
  int num_top_probes;
  int num_bot_probes;
  int *toppx;
  int *toppw;
  int *toppc;
  int *botpx;
  int *botpw;
  int *botpc;
}hallbpara;


typedef void (generatefn) (RectRedux *, void *, int , char *);
typedef void (leadgenfn) (RectRedux *, RectRedux *, int, void *);


void Patchify ( RectRedux *System, patch_para *ppara, cellDivision *cellinfo, int struc_out, char *filename);



void genLeads (RectRedux *SiteArray, RectRedux **Leads, int numleads, int mode, lead_para *params);
void HallBarify (RectRedux *System, RectRedux **Leads, hallbpara *hall_para, lead_para *params, int struc_out, char *filename);

int HallPositioning(int length2, int num_side_probes, int this_probe, int buffer_rows, int geo_renorm, int width);

void genCustomLeads (RectRedux *SiteArray, RectRedux **Leads, int numleads, lead_para *params);
void genSingleRibbonLead (RectRedux *SiteArray, RectRedux *Lead, int lead_num, void *params);
void genSingleMetalLead (RectRedux *SiteArray, RectRedux *Lead, int lead_num, void *params);

void simple558GB (RectRedux *SiteArray, void *p, int struc_out, char *filename);


void simpleRibbonGeo (RectRedux *SiteArray, void *p, int struc_out, char *filename);

void genSublatticeDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename);
void genSublatticeInterface(RectRedux *SiteArray, void *p, int struc_out, char *filename);
void genSublatticeTwoInts(RectRedux *SiteArray, void *p, int struc_out, char *filename);
void genSublatticeLeadPots(RectRedux **Leads, void *p);
void genSublatticeMoire(RectRedux *SiteArray, void *p, int struc_out, char *filename);
void genSublatticeDots(RectRedux *SiteArray, void *p, int struc_out, char *filename);

void genBubbleDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename);
void genAntidotDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename);

void genEdgeDisorderedDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename);


void exportRectConf(RectRedux *System, char *filename);
void importRectConf(RectRedux *System, int length, int length2, char *filename);


void simpleBilayerGeo (RectRedux *SiteArray, void *p, int struc_out, char *filename);
void genSingleBLGLead (RectRedux *SiteArray, RectRedux *Lead, int lead_num, void *params);


