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
}RectRedux;



void genSublatticeDevice (RectRedux *SiteArray, int buffer_rows, double a_conc, double a_pot, double b_conc, double b_pot, int seed, int struc_out, char *filename);

void exportRectConf(RectRedux *System, char *filename);
void importRectConf(RectRedux *System, int length, int length2, char *filename);