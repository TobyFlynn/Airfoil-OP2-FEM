#include "airfoil_data.h"

#include <memory>

AirfoilData::AirfoilData(int numCells) {
  nodeX_data = (double*)malloc(3 * numCells * sizeof(double));
  nodeY_data = (double*)malloc(3 * numCells * sizeof(double));
  x_data = (double *)malloc(15 * numCells * sizeof(double));
  y_data = (double *)malloc(15 * numCells * sizeof(double));
  xr_data = (double *)malloc(15 * numCells * sizeof(double));
  yr_data = (double *)malloc(15 * numCells * sizeof(double));
  xs_data = (double *)malloc(15 * numCells * sizeof(double));
  ys_data = (double *)malloc(15 * numCells * sizeof(double));
  rx_data = (double *)malloc(15 * numCells * sizeof(double));
  ry_data = (double *)malloc(15 * numCells * sizeof(double));
  sx_data = (double *)malloc(15 * numCells * sizeof(double));
  sy_data = (double *)malloc(15 * numCells * sizeof(double));
  nx_data = (double *)malloc(3 * 5 * numCells * sizeof(double));
  ny_data = (double *)malloc(3 * 5 * numCells * sizeof(double));
  fscale_data = (double *)malloc(3 * 5 * numCells * sizeof(double));
  for(int i = 0; i < 4; i++) {
    Q_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    rhs_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    F_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    G_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    dFdr_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    dFds_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    dGdr_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    dGds_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    workingQ_data[i]  = (double *)malloc(15 * numCells * sizeof(double));
    exteriorQ_data[i] = (double *)malloc(3 * 5 * numCells * sizeof(double));
    flux_data[i] = (double *)malloc(3 * 5 * numCells * sizeof(double));
    rk1_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    rk2_data[i] = (double *)malloc(15 * numCells * sizeof(double));
    rk3_data[i] = (double *)malloc(15 * numCells * sizeof(double));
  }
}

AirfoilData::~AirfoilData() {
  free(nodeX_data);
  free(nodeY_data);
  free(x_data);
  free(y_data);
  free(xr_data);
  free(yr_data);
  free(xs_data);
  free(ys_data);
  free(rx_data);
  free(ry_data);
  free(sx_data);
  free(sy_data);
  free(nx_data);
  free(ny_data);
  free(fscale_data);
  for(int i = 0; i < 4; i++) {
    free(Q_data[i]);
    free(rhs_data[i]);
    free(F_data[i]);
    free(G_data[i]);
    free(dFdr_data[i]);
    free(dFds_data[i]);
    free(dGdr_data[i]);
    free(dGds_data[i]);
    free(workingQ_data[i]);
    free(exteriorQ_data[i]);
    free(flux_data[i]);
    free(rk1_data[i]);
    free(rk2_data[i]);
    free(rk3_data[i]);
  }
}
