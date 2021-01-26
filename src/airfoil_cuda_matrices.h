#ifndef __AIRFOIL_CUDA_MATRICES_H
#define __AIRFOIL_CUDA_MATRICES_H

#include "cublas_v2.h"

extern double ones[15];
extern double r[15];
extern double s[15];
extern double Dr[15 * 15];
extern double Ds[15 * 15];

void init_grid_matrices(cublasHandle_t handle, const int numCells,
                        const double *node_coords, const int *cell2nodes,
                        double *x_d, double *y_d, double *xr_d, double *xs_d,
                        double *yr_d, double *ys_d);

void internal_fluxes_matrices(cublasHandle_t handle, const int numCells,
                              const double *F_d, const double *G_d,
                              double *dFdr_d, double *dFds_d, double *dGdr_d,
                              double *dGds_d);

#endif
