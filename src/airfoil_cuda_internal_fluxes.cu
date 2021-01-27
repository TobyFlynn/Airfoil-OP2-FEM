#include "airfoil_cuda_matrices.h"

void internal_fluxes_matrices(cublasHandle_t handle, const int numCells,
                              const double *F_d, const double *G_d,
                              double *dFdr_d, double *dFds_d, double *dGdr_d,
                              double *dGds_d) {
  double *Drw_d;
  cudaMalloc((void**)&Drw_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Drw_d, Drw, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *Dsw_d;
  cudaMalloc((void**)&Dsw_d, 15 * 15 * sizeof(double));
  cudaMemcpy(Dsw_d, Dsw, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  cudaDeviceSynchronize();
  cudaStream_t streams[NUMBER_OF_MATRIX_STREAMS];
  for(int i = 0; i < NUMBER_OF_MATRIX_STREAMS; i++) {
    cudaStreamCreate(&streams[i]);
  }
  for(int c = 0; c < numCells; c++) {
    cublasSetStream(handle, streams[c % NUMBER_OF_MATRIX_STREAMS]);
    const double *F = F_d + c * 4 * 15;
    const double *G = G_d + c * 4 * 15;
    double *dFdr = dFdr_d + c * 4 * 15;
    double *dFds = dFds_d + c * 4 * 15;
    double *dGdr = dGdr_d + c * 4 * 15;
    double *dGds = dGds_d + c * 4 * 15;

    for(int i = 0; i < 4; i++) {
      double alpha = 1.0;
      double beta = 0.0;
      // CUBLAS_OP_T because cublas is column major but constants are stored row major
      cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, F + i, 4, &beta, dFdr + i, 4);
      cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, F + i, 4, &beta, dFds + i, 4);
      cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, G + i, 4, &beta, dGdr + i, 4);
      cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, G + i, 4, &beta, dGds + i, 4);
    }
  }

  for(int i = 0; i < NUMBER_OF_MATRIX_STREAMS; i++) {
    cudaStreamDestroy(streams[i]);
  }

  cudaFree(Drw_d);
  cudaFree(Dsw_d);
  cudaDeviceSynchronize();
}
