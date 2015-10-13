#pragma once
#include "devices.h"
#include "connect.h"

typedef double _Complex (hoppingfunc) (RectRedux *,  int, int, void *p);

typedef struct {
	double t0;		
	int isperiodic;
	double kpar;
	double NN_lowdis;
	double NN_highdis;
}simpleTB_params;


double genConduc4(double _Complex En, RectRedux *DeviceCell, double *ldoses, double **conds, int *indices, int **neigh, double hopping);
double genConduc5(double _Complex En, RectRedux *DeviceCell, double hopping);

void genDeviceGF(double _Complex En, RectRedux *DeviceCell, cnxProfile *cnxp, 
		      cellDivision *cellinfo, hoppingfunc *hoppingfn, void *hoppingparams, int mode, 
		      double _Complex **Gon, double _Complex **Goff, double _Complex **Sigma);

double _Complex simpleTB(RectRedux *DeviceCell, int a, int b, void *hoppingparams);
