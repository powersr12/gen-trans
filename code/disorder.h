#pragma once
#include "devices.h"


void potentialDisorder (RectRedux *DeviceCell, void *p, int disprof_out, char *filename );

typedef struct {
	double conc;
	double delta;
	double xi;
	int seed;
}potDis_para;