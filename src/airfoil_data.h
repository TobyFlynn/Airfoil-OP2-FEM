#ifndef __AIRFOIL_DATA_H
#define __AIRFOIL_DATA_H

class AirfoilData {
public:
  AirfoilData(int numCells);
  ~AirfoilData();
  double *nodeX_data;
  double *nodeY_data;
  double *x_data;
  double *y_data;
  double *xr_data;
  double *yr_data;
  double *xs_data;
  double *ys_data;
  double *rx_data;
  double *ry_data;
  double *sx_data;
  double *sy_data;
  double *nx_data;
  double *ny_data;
  double *fscale_data;
  double *q_data;
  double *F_data;
  double *G_data;
  double *dFdr_data;
  double *dFds_data;
  double *dGdr_data;
  double *dGds_data;
  double *workingQ_data;
  double *exteriorQ_data;
  double *flux_data;
  double *rk1_data;
  double *rk2_data;
  double *rk3_data;
};

#endif
