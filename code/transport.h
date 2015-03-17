#pragma once
#include "devices.h"

double genConduc4(double _Complex En, RectRedux *DeviceCell, double *ldoses, double **conds, int *indices, int **neigh, double hopping);
double genConduc5(double _Complex En, RectRedux *DeviceCell, double hopping);

