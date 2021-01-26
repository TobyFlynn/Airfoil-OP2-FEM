#include "airfoil_cuda_matrices.h"

void init_grid_matrices(cublasHandle_t handle, const int numCells,
                        const double *node_coords, const int *cell2nodes,
                        double *x_d, double *y_d, double *xr_d, double *xs_d,
                        double *yr_d, double *ys_d) {
  cudaDeviceSynchronize();
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
    // Get nodes for this cell (on host)
    const double *n0 = &node_coords[2 * cell2nodes[3 * c]];
    const double *n1 = &node_coords[2 * cell2nodes[3 * c + 1]];
    const double *n2 = &node_coords[2 * cell2nodes[3 * c + 2]];

    double *temp = temp_d + c * 15;
    double *x = x_d + c * 15;
    double *y = y_d + c * 15;
    double *xr = xr_d + c * 15;
    double *xs = xs_d + c * 15;
    double *yr = yr_d + c * 15;
    double *ys = ys_d + c * 15;

    double alpha = 1.0;
    cublasDcopy(handle, 15, ones_d, 1, x, 1);
    cublasDaxpy(handle, 15, &alpha, r_d, 1, x, 1);
    alpha = 0.5 * n1[0];
    cublasDscal(handle, 15, &alpha, x, 1);
    cublasDcopy(handle, 15, ones_d, 1, temp, 1);
    alpha = 1.0;
    cublasDaxpy(handle, 15, &alpha, s_d, 1, temp, 1);
    alpha = 0.5 * n2[0];
    cublasDaxpy(handle, 15, &alpha, temp, 1, x, 1);
    cublasDcopy(handle, 15, s_d, 1, temp, 1);
    alpha = 1.0;
    cublasDaxpy(handle, 15, &alpha, r_d, 1, temp, 1);
    alpha = -0.5 * n0[0];
    cublasDaxpy(handle, 15, &alpha, temp, 1, x, 1);

    cublasDcopy(handle, 15, ones_d, 1, y, 1);
    alpha = 1.0;
    cublasDaxpy(handle, 15, &alpha, r_d, 1, y, 1);
    alpha = 0.5 * n1[1];
    cublasDscal(handle, 15, &alpha, y, 1);
    cublasDcopy(handle, 15, ones_d, 1, temp, 1);
    alpha = 1.0;
    cublasDaxpy(handle, 15, &alpha, s_d, 1, temp, 1);
    alpha = 0.5 * n2[1];
    cublasDaxpy(handle, 15, &alpha, temp, 1, y, 1);
    cublasDcopy(handle, 15, s_d, 1, temp, 1);
    alpha = 1.0;
    cublasDaxpy(handle, 15, &alpha, r_d, 1, temp, 1);
    alpha = -0.5 * n0[1];
    cublasDaxpy(handle, 15, &alpha, temp, 1, y, 1);

    // CUBLAS_OP_T because cublas is column major but constants are stored row major
    // xr = Dr * x
    alpha = 1.0;
    double beta = 0.0;
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dr_d, 15, x, 1, &beta, xr, 1);
    // xs = Ds * x
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Ds_d, 15, x, 1, &beta, xs, 1);
    // yr = Dr * y
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Dr_d, 15, y, 1, &beta, yr, 1);
    // ys = Ds * y
    cublasDgemv(handle, CUBLAS_OP_T, 15, 15, &alpha, Ds_d, 15, y, 1, &beta, ys, 1);
  }

  cudaFree(ones_d);
  cudaFree(r_d);
  cudaFree(temp_d);
  cudaFree(Dr_d);
  cudaFree(Ds_d);
  cudaDeviceSynchronize();
}
