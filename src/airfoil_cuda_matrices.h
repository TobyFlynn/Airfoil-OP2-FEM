#ifndef __AIRFOIL_CUDA_MATRICES_H
#define __AIRFOIL_CUDA_MATRICES_H

#define NUMBER_OF_MATRIX_STREAMS 8

#include "cublas_v2.h"

extern double ones[15];
extern double r[15];
extern double s[15];
extern double Dr[15 * 15];
extern double Ds[15 * 15];
extern double Drw[15 * 15];
extern double Dsw[15 * 15];
extern double LIFT[15 * 15];

void init_grid_matrices(cublasHandle_t handle, const int numCells,
                        const double *node_coords, const int *cell2nodes,
                        double *x_d, double *y_d, double *xr_d, double *xs_d,
                        double *yr_d, double *ys_d);

void internal_fluxes_matrices(cublasHandle_t handle, const int numCells,
                              const double *f0_d, const double *f1_d,
                              const double *f2_d, const double *f3_d,
                              const double *g0_d, const double *g1_d,
                              const double *g2_d, const double *g3_d,
                              double *dFdr0_d, double *dFdr1_d, double *dFdr2_d,
                              double *dFdr3_d, double *dFds0_d, double *dFds1_d,
                              double *dFds2_d, double *dFds3_d, double *dGdr0_d,
                              double *dGdr1_d, double *dGdr2_d, double *dGdr3_d,
                              double *dGds0_d, double *dGds1_d, double *dGds2_d,
                              double *dGds3_d);

void face_fluxes_matrices(cublasHandle_t handle, const int numCells,
                          const double *flux0_d, const double *flux1_d,
                          const double *flux2_d, const double *flux3_d,
                          double *qRHS0_d, double *qRHS1_d, double *qRHS2_d,
                          double *qRHS3_d);
#endif
