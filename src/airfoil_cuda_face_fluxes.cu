#include "airfoil_cuda_matrices.h"

void face_fluxes_matrices(cublasHandle_t handle, const int numCells,
                          const double *flux0_d, const double *flux1_d,
                          const double *flux2_d, const double *flux3_d,
                          double *qRHS0_d, double *qRHS1_d, double *qRHS2_d,
                          double *qRHS3_d) {
  double *LIFT_d;
  cudaMalloc((void**)&LIFT_d, 15 * 15 * sizeof(double));
  cudaMemcpy(LIFT_d, LIFT, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double alpha = -1.0;
  double beta = 1.0;

  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, LIFT_d, 15, flux0_d, 15, &beta, qRHS0_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, LIFT_d, 15, flux1_d, 15, &beta, qRHS1_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, LIFT_d, 15, flux2_d, 15, &beta, qRHS2_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, LIFT_d, 15, flux3_d, 15, &beta, qRHS3_d, 15);

  cudaDeviceSynchronize();
  cudaFree(LIFT_d);
}
