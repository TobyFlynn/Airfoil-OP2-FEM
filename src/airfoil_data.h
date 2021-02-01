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
  double *Q_data[4];
  double *F_data[4];
  double *G_data[4];
  double *dFdr_data[4];
  double *dFds_data[4];
  double *dGdr_data[4];
  double *dGds_data[4];
  double *workingQ_data[4];
  double *exteriorQ_data[4];
  double *flux_data[4];
  double *rk1_data[4];
  double *rk2_data[4];
  double *rk3_data[4];
  double *rhs_data[4];
};

#endif
