#pragma once

typedef int (cnxRulesFn) (RectRedux *, void *, int, int);
typedef int (cnxRulesFn2) (RectRedux *, RectRedux *, void *, int, int);

// 
//says whether or not two atoms are connected, NOT the hopping value
typedef struct {
	int max_neigh;		
	int *site_cnxnum;
	int **site_cnx;
}cnxProfile;



typedef struct {
	double conn_sep_thresh_min;
	double conn_sep_thresh_max;
	int periodic;
}graph_conn_para;

typedef struct {
	double intra_thresh_min;
	double intra_thresh_max;
	double inter_thresh_min;
	double inter_thresh_max;
	double zthresh1;	//for deciding if two atoms in same layer
	double zthresh2;	//max z sep for neighbouring layers.
	int periodic;
}blg_conn_para;


typedef struct {
  cnxRulesFn2 *rule;
  void *rule_params;
  int num_leads;
  RectRedux **Leads;
}gen_start_params;
  
typedef struct {
  int num_leads;
  double *startx;
  double *endx;
}multix_start_params;

typedef struct {
  cnxRulesFn2 *rule;
  void *rule_params;
  int num_leads;
  RectRedux **Leads;
  int *leadtype;
}custom_start_params;

void device_connectivity (RectRedux *DeviceCell, cnxRulesFn *rule, void *rule_params, cnxProfile *cnxp);
void printConnectivity (RectRedux *DeviceCell, cnxProfile *cnxp);

void genStartingCell (RectRedux *DeviceCell, cellDivision *cellinfo, int config, void *start_params);
void cellSplitter (RectRedux *DeviceCell, cnxProfile *cnxp, cellDivision *cellinfo);
void printOutLeadStrucs(RectRedux *DeviceCell, RectRedux **Leads, cellDivision *cellinfo, char *output);



int zzacnn (RectRedux *DeviceCell, void *rule_params, int a, int b);
int zzacnnk (RectRedux *DeviceCell, void *rule_params, int a, int b);

int graph_conn_sep (RectRedux *DeviceCell, void *rule_params, int a, int b);
int graph_conn_sep2 (RectRedux *DeviceCell, RectRedux *LeadCell, void *rule_params, int a, int b);
int blg_conn_sep2 (RectRedux *DeviceCell, RectRedux *LeadCell, void *rule_params, int a, int b);
int blg_conn_sep (RectRedux *DeviceCell, void *rule_params, int a, int b);

int mlg_conn_sep2 (RectRedux *DeviceCell, RectRedux *LeadCell, void *rule_params, int a, int b);
int mlg_conn_sep (RectRedux *DeviceCell, void *rule_params, int a, int b);

