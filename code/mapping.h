#pragma once
#include "devices.h"

void LDOSMapOut(RectRedux *DeviceCell, double *ldoses, char *out);
void CurrentMapOut(RectRedux *DeviceCell, double **conds, int *indices, int **neigh, char *out);

