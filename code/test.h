
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "mapping.h"
#include "transport.h"

#include "../libs/matrices.h"
#include "../libs/hubbard.h"
#include "../libs/graphene_gf.h"
#include "../libs/suscept.h"
#include "../libs/static.h"


typedef int (cnxRulesFn) (RectRedux *, void *, int, int);

//says whether or not two atoms are connected, NOT the hopping value
typedef struct {
	int max_neigh;		
	int *site_cnxnum;
	int **site_cnx;
}cnxProfile;

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


void device_connectivity (RectRedux *DeviceCell, cnxRulesFn *rule, void *rule_params, cnxProfile *cnxp);
void printConnectivity (RectRedux *DeviceCell, cnxProfile *cnxp);

void genStartingCell (RectRedux *DeviceCell, cellDivision *cellinfo, int config, void *start_params);
void cellSplitter (RectRedux *DeviceCell, cnxProfile *cnxp, cellDivision *cellinfo);

int zzacnn (RectRedux *DeviceCell, void *rule_params, int i, int j);
