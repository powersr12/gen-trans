#pragma once

#include <stdio.h>
#include "../libs/matrices.h"

#include "useful.h"



typedef struct {
	int geo;		//=0 for zigzag, =1 for armchair
	int length;
	int length2;
	int **siteinfo;		//in lattice? which sublattice?
	double **pos;
	double *site_pots;
	int **chaininfo;
	int *Nrem;
	int *Ntot;
}RectRedux;

typedef struct {
	double **shift_vecs;
	void *hopfn;
	void *leadsfn;
	void *hoppara;
}lead_para;

void genLeads (RectRedux *SiteArray, RectRedux **Leads, int numleads, int mode, lead_para *params);

void **simpleRibbonGeo (RectRedux *SiteArray, void *p, int struc_out, char *filename);

void genSublatticeDevice (RectRedux *SiteArray, int buffer_rows, double a_conc, double a_pot, double b_conc, double b_pot, int seed, int struc_out, char *filename);
void genAntidotDevice (RectRedux *SiteArray, int buffer_rows, int AD_length, double AD_rad, int lat_width, int lat_length, int seed, int struc_out, char *filename);


void exportRectConf(RectRedux *System, char *filename);
void importRectConf(RectRedux *System, int length, int length2, char *filename);