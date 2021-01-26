#include "airfoil_cuda_matrices.h"

void internal_fluxes_matrices(cublasHandle_t handle, const int numCells,
                              const double *F_d, const double *G_d,
                              double *dFdr_d, double *dFds_d, double *dGdr_d,
                              double *dGds_d) {
  double *ones_d;
  cudaMalloc((void**)&ones_d, 15 * sizeof(double));
  cudaMemcpy(ones_d, ones, 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *r_d;
  cudaMalloc((void**)&r_d, 15 * sizeof(double));
  cudaMemcpy(r_d, r, 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *s_d;
  cudaMalloc((void**)&s_d, 15 * sizeof(double));
  cudaMemcpy(s_d, s, 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *temp_d;
  cudaMalloc((void**)&temp_d, numCells * 15 * sizeof(double));

  double *Dr_d;
  cudaMalloc((void**)&Dr_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Dr_d, Dr, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *Ds_d;
  cudaMalloc((void**)&Ds_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Ds_d, Ds, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  for(int c = 0; c < numCells; c++) {
    const double *F = F_d + c * 4 * 15;
    const double *G = G_d + c * 4 * 15;
    double *dFdr = dFdr_d + c * 4 * 15;
    double *dFds = dFds_d + c * 4 * 15;
    double *dGdr = dGdr_d + c * 4 * 15;
    double *dGds = dGds_d + c * 4 * 15;

    for(int i = 0; i < 4; i++) {
      cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Drw, 15, &F[i], 4, 0.0, dFdr, 1);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dsw, 15, &F[i], 4, 0.0, dFds, 1);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Drw, 15, &G[i], 4, 0.0, dGdr, 1);
      cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dsw, 15, &G[i], 4, 0.0, dGds, 1);
    }
  }

  cudaFree(ones_d);
  cudaFree(r_d);
  cudaFree(temp_d);
  cudaFree(Dr_d);
  cudaFree(Ds_d);
}
