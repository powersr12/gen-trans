
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


typedef int (*cnxRulesFn) (RectRedux *, void *, int, int);

//says whether or not two atoms are connected, NOT the hopping value
typedef struct {
	int max_neigh;		//=0 for zigzag, =1 for armchair
	int *site_cnxnum;
	int **site_cnx;
}cnxProfile;


void device_connectivity (RectRedux *DeviceCell, cnxRulesFn *rule, void *rule_params, cnxProfile *cnxp);

int zzacnn (RectRedux *DeviceCell, void *rule_params, int i, int j);
