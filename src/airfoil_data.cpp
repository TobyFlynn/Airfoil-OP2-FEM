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
  q_data  = (double *)malloc(4 * 15 * numCells * sizeof(double));
  F_data  = (double *)malloc(4 * 15 * numCells * sizeof(double));
  G_data  = (double *)malloc(4 * 15 * numCells * sizeof(double));
  dFdr_data  = (double *)malloc(4 * 15 * numCells * sizeof(double));
  dFds_data  = (double *)malloc(4 * 15 * numCells * sizeof(double));
  dGdr_data  = (double *)malloc(4 * 15 * numCells * sizeof(double));
  dGds_data  = (double *)malloc(4 * 15 * numCells * sizeof(double));
  workingQ_data  = (double *)malloc(4 * 15 * numCells * sizeof(double));
  exteriorQ_data = (double *)malloc(4 * 3 * 5 * numCells * sizeof(double));
  flux_data = (double *)malloc(4 * 3 * 5 * numCells * sizeof(double));
  rk1_data = (double *)malloc(4 * 15 * numCells * sizeof(double));
  rk2_data = (double *)malloc(4 * 15 * numCells * sizeof(double));
  rk3_data = (double *)malloc(4 * 15 * numCells * sizeof(double));
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
  free(q_data);
  free(F_data);
  free(G_data);
  free(dFdr_data);
  free(dFds_data);
  free(dGdr_data);
  free(dGds_data);
  free(workingQ_data);
  free(exteriorQ_data);
  free(flux_data);
  free(rk1_data);
  free(rk2_data);
  free(rk3_data);
}
