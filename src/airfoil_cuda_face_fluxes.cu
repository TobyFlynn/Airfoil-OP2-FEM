#include "airfoil_cuda_matrices.h"

void face_fluxes_matrices(cublasHandle_t handle, const int numCells,
                          const double *flux0_d, const double *flux1_d,
                          const double *flux2_d, const double *flux3_d,
                          double *qRHS0_d, double *qRHS1_d, double *qRHS2_d,
                          double *qRHS3_d) {
  double *LIFT_d;
  cudaMalloc((void**)&LIFT_d, 15 * 15 * sizeof(double));
  cudaMemcpy(LIFT_d, LIFT, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();
  // cudaStream_t streams[NUMBER_OF_MATRIX_STREAMS];
  // for(int i = 0; i < NUMBER_OF_MATRIX_STREAMS; i++) {
  //   cudaStreamCreate(&streams[i]);
  // }
  for(int c = 0; c < numCells; c++) {
    // cublasSetStream(handle, streams[c % NUMBER_OF_MATRIX_STREAMS]);
    const double *flux0 = flux0_d + c * 3 * 5;
    const double *flux1 = flux1_d + c * 3 * 5;
    const double *flux2 = flux2_d + c * 3 * 5;
    const double *flux3 = flux3_d + c * 3 * 5;
    double *qRHS0 = qRHS0_d + c * 15;
    double *qRHS1 = qRHS1_d + c * 15;
    double *qRHS2 = qRHS2_d + c * 15;
    double *qRHS3 = qRHS3_d + c * 15;

    double alpha = -1.0;
    double beta = 1.0;

    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, LIFT_d, 15, flux0, 1, &beta, qRHS0, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, LIFT_d, 15, flux1, 1, &beta, qRHS1, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, LIFT_d, 15, flux2, 1, &beta, qRHS2, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, LIFT_d, 15, flux3, 1, &beta, qRHS3, 1);

    // for(int i = 0; i < 4; i++) {
    //   double alpha = -1.0;
    //   double beta = 1.0;
    //   cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, LIFT_d, 15, flux + i, 4, &beta, qRHS + i, 4);
    // }
  }

  // for(int i = 0; i < NUMBER_OF_MATRIX_STREAMS; i++) {
  //   cudaStreamDestroy(streams[i]);
  // }

  cudaFree(LIFT_d);
  cudaDeviceSynchronize();
}
