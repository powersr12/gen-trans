#pragma once
#include "devices.h"
#include "connect.h"

typedef double _Complex (hoppingfunc) (RectRedux *,  RectRedux *, int, int, double *, void *p);
typedef void (leadfunction) (double _Complex, RectRedux *, RectRedux **, cellDivision *, void *, double _Complex **); 




// typedef struct {
// 	double t0;		
// 	int isperiodic;
// 	double kpar;
// 	double NN_lowdis;
// 	double NN_highdis;
// }simpleTB_params;

typedef struct {
	int num_neigh;
	double _Complex *hops;		
	int isperiodic;
	double kpar;
	double *NN_lowdis;
	double *NN_highdis;
	int gauge;
	double Btes;
	int *restrics;	//are there limits on x, y with field
	double **limits; //what are the limits on x, y with field
}gen_hop_params;


typedef struct {
	int num_leads;
	int TRsym;
	double **transmissions;
}trans_params;


double genConduc4(double _Complex En, RectRedux *DeviceCell, double *ldoses, double **conds, int *indices, int **neigh, double hopping);
double genConduc5(double _Complex En, RectRedux *DeviceCell, double hopping);

void genDeviceGF(double _Complex En, RectRedux *DeviceCell, cnxProfile *cnxp, 
		      cellDivision *cellinfo, hoppingfunc *hoppingfn, void *hoppingparams, int mode, 
		      double _Complex **Gon, double _Complex **Goff, double _Complex **Sigma);

double _Complex simpleTB(RectRedux *aDeviceCell, RectRedux *bDeviceCell, int a, int b, double *bshifts, void *hoppingparams);
double _Complex peierlsTB(RectRedux *aDeviceCell, RectRedux *bDeviceCell, int a, int b, double *bshifts, void *hoppingparams);

double _Complex graphenePeierlsPhase(double x1, double y1, double x2, double y2, int gauge, double BTesla, int *res, double **limits);



void simple2leads (double _Complex En, RectRedux *DeviceCell, RectRedux **LeadCells, cellDivision *cellinfo, lead_para *params, double _Complex **Sigma);
void lead_prep(double _Complex En, RectRedux *LeadCell, int leadindex, lead_para *params, double _Complex **ginv, double _Complex **V12, double _Complex **V21);
void genTransmissions(double _Complex En, RectRedux *DeviceCell, RectRedux **Leads, cnxProfile *cnxp, 
		      cellDivision *cellinfo, hoppingfunc *hoppingfn, void *hoppingparams,
		      lead_para *leadsparams, trans_params *tpara);