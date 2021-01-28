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

  cudaDeviceSynchronize();
  // cudaStream_t streams[NUMBER_OF_MATRIX_STREAMS];
  // for(int i = 0; i < NUMBER_OF_MATRIX_STREAMS; i++) {
  //   cudaStreamCreate(&streams[i]);
  // }
  for(int c = 0; c < numCells; c++) {
    const double *f0 = f0_d + c * 15;
    const double *f1 = f1_d + c * 15;
    const double *f2 = f2_d + c * 15;
    const double *f3 = f3_d + c * 15;
    const double *g0 = g0_d + c * 15;
    const double *g1 = g1_d + c * 15;
    const double *g2 = g2_d + c * 15;
    const double *g3 = g3_d + c * 15;
    double *dFdr0 = dFdr0_d + c * 15;
    double *dFdr1 = dFdr1_d + c * 15;
    double *dFdr2 = dFdr2_d + c * 15;
    double *dFdr3 = dFdr3_d + c * 15;
    double *dFds0 = dFds0_d + c * 15;
    double *dFds1 = dFds1_d + c * 15;
    double *dFds2 = dFds2_d + c * 15;
    double *dFds3 = dFds3_d + c * 15;
    double *dGdr0 = dGdr0_d + c * 15;
    double *dGdr1 = dGdr1_d + c * 15;
    double *dGdr2 = dGdr2_d + c * 15;
    double *dGdr3 = dGdr3_d + c * 15;
    double *dGds0 = dGds0_d + c * 15;
    double *dGds1 = dGds1_d + c * 15;
    double *dGds2 = dGds2_d + c * 15;
    double *dGds3 = dGds3_d + c * 15;

    double alpha = 1.0;
    double beta = 0.0;

    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, f0, 1, &beta, dFdr0, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, f1, 1, &beta, dFdr1, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, f2, 1, &beta, dFdr2, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, f3, 1, &beta, dFdr3, 1);

    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, f0, 1, &beta, dFds0, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, f1, 1, &beta, dFds1, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, f2, 1, &beta, dFds2, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, f3, 1, &beta, dFds3, 1);

    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, g0, 1, &beta, dGdr0, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, g1, 1, &beta, dGdr1, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, g2, 1, &beta, dGdr2, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, g3, 1, &beta, dGdr3, 1);

    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, g0, 1, &beta, dGds0, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, g1, 1, &beta, dGds1, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, g2, 1, &beta, dGds2, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, g3, 1, &beta, dGds3, 1);

    // cublasSetStream(handle, streams[(2*c) % NUMBER_OF_MATRIX_STREAMS]);
    // for(int i = 0; i < 4; i++) {
    //   double alpha = 1.0;
    //   double beta = 0.0;
    //   // CUBLAS_OP_T because cublas is column major but constants are stored row major
    //   cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, F + i, 4, &beta, dFdr + i, 4);
    //   cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, F + i, 4, &beta, dFds + i, 4);
    //   // cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, G + i, 4, &beta, dGdr + i, 4);
    //   // cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, G + i, 4, &beta, dGds + i, 4);
    // }
    //
    // cublasSetStream(handle, streams[(2*c + 1) % NUMBER_OF_MATRIX_STREAMS]);
    // for(int i = 0; i < 4; i++) {
    //   double alpha = 1.0;
    //   double beta = 0.0;
    //   // CUBLAS_OP_T because cublas is column major but constants are stored row major
    //   // cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, F + i, 4, &beta, dFdr + i, 4);
    //   // cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, F + i, 4, &beta, dFds + i, 4);
    //   cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Drw_d, 15, G + i, 4, &beta, dGdr + i, 4);
    //   cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dsw_d, 15, G + i, 4, &beta, dGds + i, 4);
    // }
  }

  // for(int i = 0; i < NUMBER_OF_MATRIX_STREAMS; i++) {
  //   cudaStreamDestroy(streams[i]);
  // }

  cudaFree(Drw_d);
  cudaFree(Dsw_d);
  cudaDeviceSynchronize();
}
