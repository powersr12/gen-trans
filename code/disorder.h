#pragma once


void potentialDisorder (RectRedux *DeviceCell, void *p, int disprof_out, char *filename );
void customEdgePots (RectRedux *DeviceCell, void *p);
void cap_potential(RectRedux *DeviceCell, double width);
void cap_potential2(RectRedux *DeviceCell, double width);



typedef struct {
	double conc;
	double delta;
	double xi;
	int seed;
}potDis_para;

typedef struct {
	int type;
	double AT1;
	double AT2;
	double AT3;
	double AT4;
	double BT1;
	double BT2;
	double BT3;
	double BT4;
	double AB1;
	double AB2;
	double AB3;
	double AB4;
	double BB1;
	double BB2;
	double BB3;
	double BB4;
}cedgepot_para;