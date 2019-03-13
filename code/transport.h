#pragma once



typedef double _Complex (hoppingfunc) (RectRedux *,  RectRedux *, int, int, double *, void *p);
typedef void (leadfunction) (double _Complex, RectRedux *, RectRedux **, cellDivision *, void *, double _Complex **); 
typedef void (singleleadfunction) (int, double _Complex, RectRedux *, RectRedux **, cellDivision *, void *, double _Complex **); 



// typedef struct {
// 	double t0;		
// 	int isperiodic;
// 	double kpar;
// 	double NN_lowdis;
// 	double NN_highdis;
// }simpleTB_params;

typedef struct {
	int num_neigh;
	int *max_neigh;
	double _Complex *hops;		
	int isperiodic;
	double kpar;
	double *NN_lowdis;
	double *NN_highdis;
	double *NN_shifts;
	double *NN_zmin;
	double *NN_zmax;
	int gauge;
	double Btes;
	int *restrics;	//are there limits on x, y with field
	double **limits; //what are the limits on x, y with field
}gen_hop_params;


typedef struct {
	int num_leads;
	int TRsym;
	double **transmissions;
        char *filename;
        int ispatched;
        void *patchpara;
}trans_params;


double genConduc4(double _Complex En, RectRedux *DeviceCell, double *ldoses, double **conds, int *indices, int **neigh, double hopping);
double genConduc5(double _Complex En, RectRedux *DeviceCell, double hopping);

void genDeviceGF(double _Complex En, RectRedux *DeviceCell, cnxProfile *cnxp, 
		      cellDivision *cellinfo, hoppingfunc *hoppingfn, void *hoppingparams, int mode, int mode2, 
		      double _Complex **Gon, double _Complex *Gdiags, double _Complex **Goff, double _Complex **Sigma);

double _Complex simpleTB(RectRedux *aDeviceCell, RectRedux *bDeviceCell, int a, int b, double *bshifts, void *hoppingparams);
double _Complex peierlsTB(RectRedux *aDeviceCell, RectRedux *bDeviceCell, int a, int b, double *bshifts, void *hoppingparams);
double _Complex strainedTB(RectRedux *aDeviceCell, RectRedux *bDeviceCell, int a, int b, double *bshifts, void *hoppingparams);


double _Complex graphenePeierlsPhase(double x1, double y1, double x2, double y2, int gauge, double BTesla, int *res, double **limits);


void multipleLeads (double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, lead_para *params, double _Complex **Sigma);

void multipleSimplestMetalLeads (double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, lead_para *params, double _Complex **Sigma);

void simple2leads (double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, lead_para *params, double _Complex **Sigma);
void lead_prep(double _Complex En, RectRedux *LeadCell, int leadindex, lead_para *params, double _Complex **ginv, double _Complex **V12, double _Complex **V21);
void lead_prep2(double _Complex En, RectRedux *LeadCell, int leadindex, rib_lead_para *params, double _Complex **ginv, double _Complex **V12, double _Complex **V21);

void genTransmissions(double _Complex En, RectRedux *DeviceCell, RectRedux **Leads, cnxProfile *cnxp, 
		      cellDivision *cellinfo, hoppingfunc *hoppingfn, void *hoppingparams,
		      lead_para *leadsparams, trans_params *tpara, int mode, double *ldoses, double ***currents);



void genKXbandproj(RectRedux *DeviceCell,  hoppingfunc *hoppingfn, void *hoppingparams, int mode,
		      double kx, double *bands, double **projs, double **weights);


void gate_induced_pot ( int vgtype, RectRedux *DeviceCell, double *engdeppots, double gate_voltage, double edge_cut_off, double subs_thick, double subs_epsr);


void multipleCustomLeads (double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, lead_para *params, double _Complex **Sigma);

void singleRibbonLead (int leadnum, double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, void *params, double _Complex **Sigma);

void singleSimplestMetalLead (int leadnum, double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, void *params, double _Complex **Sigma);