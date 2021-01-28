#include "airfoil_cuda_matrices.h"

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
                              double *dGds3_d) {
  double *Drw_d;
  cudaMalloc((void**)&Drw_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Drw_d, Drw, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *Dsw_d;
  cudaMalloc((void**)&Dsw_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Dsw_d, Dsw, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Drw_d, 15, f0_d, 15, &beta, dFdr0_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Drw_d, 15, f1_d, 15, &beta, dFdr1_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Drw_d, 15, f2_d, 15, &beta, dFdr2_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Drw_d, 15, f3_d, 15, &beta, dFdr3_d, 15);

  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Drw_d, 15, g0_d, 15, &beta, dGdr0_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Drw_d, 15, g1_d, 15, &beta, dGdr1_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Drw_d, 15, g2_d, 15, &beta, dGdr2_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Drw_d, 15, g3_d, 15, &beta, dGdr3_d, 15);

  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Dsw_d, 15, f0_d, 15, &beta, dFds0_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Dsw_d, 15, f1_d, 15, &beta, dFds1_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Dsw_d, 15, f2_d, 15, &beta, dFds2_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Dsw_d, 15, f3_d, 15, &beta, dFds3_d, 15);

  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Dsw_d, 15, g0_d, 15, &beta, dGds0_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Dsw_d, 15, g1_d, 15, &beta, dGds1_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Dsw_d, 15, g2_d, 15, &beta, dGds2_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, Dsw_d, 15, g3_d, 15, &beta, dGds3_d, 15);

  cudaDeviceSynchronize();
  cudaFree(Drw_d);
  cudaFree(Dsw_d);
}
