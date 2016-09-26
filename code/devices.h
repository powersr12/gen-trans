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

typedef struct {
	int buffer_rows;
	double a_conc;
	double a_pot;
	double b_conc;
	double b_pot;
	int seed;
}subl_para;

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

void genLeads (RectRedux *SiteArray, RectRedux **Leads, int numleads, int mode, lead_para *params);
void HallBarify (RectRedux *System, RectRedux **Leads, hallbpara *hall_para, lead_para *params, int struc_out, char *filename);

void simpleRibbonGeo (RectRedux *SiteArray, void *p, int struc_out, char *filename);

void genSublatticeDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename);
void genSublatticeInterface(RectRedux *SiteArray, void *p, int struc_out, char *filename);
void genSublatticeTwoInts(RectRedux *SiteArray, void *p, int struc_out, char *filename);

void genAntidotDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename);

void genEdgeDisorderedDevice(RectRedux *SiteArray, void *p, int struc_out, char *filename);


void exportRectConf(RectRedux *System, char *filename);
void importRectConf(RectRedux *System, int length, int length2, char *filename);