#include "airfoil_cuda_matrices.h"

void face_fluxes_matrices(cublasHandle_t handle, const int numCells,
                          const double *flux_d, double *qRHS_d) {
  cudaDeviceSynchronize();
  double *LIFT_d;
  cudaMalloc((void**)&LIFT_d, 15 * 15 * sizeof(double));
  cudaMemcpy(LIFT_d, LIFT, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  for(int c = 0; c < numCells; c++) {
    const double *flux = flux_d + c * 4 * 3 * 5;
    double *qRHS = qRHS_d + c * 4 * 15;

    for(int i = 0; i < 4; i++) {
      double alpha = -1.0;
      double beta = 1.0;
      cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, LIFT_d, 15, flux + i, 4, &beta, qRHS + i, 4);
    }
  }

  cudaFree(LIFT_d);
  cudaDeviceSynchronize();
}
