// #include <cblas.h>

inline void init_grid(const double *n0, const double *n1, const double *n2,
                      double *nodeX, double *nodeY, double *x, double *y,
                      double *xr, double *yr, double *xs, double *ys,
                      double *rx, double *ry, double *sx, double *sy,
                      double *nx, double *ny, double *fscale) {
  // Calculate the solution point coordinates

  nodeX[0] = n0[0];
  nodeX[1] = n1[0];
  nodeX[2] = n2[0];
  nodeY[0] = n0[1];
  nodeY[1] = n1[1];
  nodeY[2] = n2[1];
  // x = 0.5(-n0_x (r + s) + n1_x(1 + r) + n2_x(1 + s))
  // x = 0.5*n1_x * (1 + r) + 0.5*n2_x * (1 + s) - 0.5*n0_x * (r + s)

  // x <- 0.5*n1_x * (1 + r)
  // cblas_dcopy(15, ones, 1, x, 1);
  // cblas_daxpy(15, 1.0, r, 1, x, 1);
  // cblas_dscal(15, 0.5 * n1[0], x, 1);
  // // temp <- 1 + s
  // double temp[15];
  // cblas_dcopy(15, ones, 1, temp, 1);
  // cblas_daxpy(15, 1.0, s, 1, temp, 1);
  // // x <- 0.5*n2_x * temp + x (i.e. 0.5*n1_x * (1 + r) + 0.5*n2_x * (1 + s))
  // cblas_daxpy(15, 0.5 * n2[0], temp, 1, x, 1);
  // // temp <- r + s
  // cblas_dcopy(15, s, 1, temp, 1);
  // cblas_daxpy(15, 1.0, r, 1, temp, 1);
  // // x <- -0.5 * n0_x * temp + x (i.e. final target)
  // cblas_daxpy(15, -0.5 * n0[0], temp, 1, x, 1);
  //
  // // y = 0.5(-n0_y (r + s) + n1_y(1 + r) + n2_y(1 + s))
  // // y = 0.5*n1_y * (1 + r) + 0.5*n2_y * (1 + s) - 0.5*n0_y * (r + s)
  //
  // // y <- 0.5*n1_y * (1 + r)
  // cblas_dcopy(15, ones, 1, y, 1);
  // cblas_daxpy(15, 1.0, r, 1, y, 1);
  // cblas_dscal(15, 0.5 * n1[1], y, 1);
  // // temp <- 1 + s
  // cblas_dcopy(15, ones, 1, temp, 1);
  // cblas_daxpy(15, 1.0, s, 1, temp, 1);
  // // y <- 0.5*n2_y * temp + y (i.e. 0.5*n1_y * (1 + r) + 0.5*n2_y * (1 + s))
  // cblas_daxpy(15, 0.5 * n2[1], temp, 1, y, 1);
  // // temp <- r + s
  // cblas_dcopy(15, s, 1, temp, 1);
  // cblas_daxpy(15, 1.0, r, 1, temp, 1);
  // // x <- -0.5 * n0_y * temp + x (i.e. final target)
  // cblas_daxpy(15, -0.5 * n0[1], temp, 1, y, 1);
  //
  // // Calculate geometric factors
  //
  // // xr = Dr * x
  // cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dr, 15, x, 1, 0.0, xr, 1);
  //
  // // xs = Ds * x
  // cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Ds, 15, x, 1, 0.0, xs, 1);
  //
  // // yr = Dr * y
  // cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dr, 15, y, 1, 0.0, yr, 1);
  //
  // // ys = Ds * y
  // cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Ds, 15, y, 1, 0.0, ys, 1);

  // J = -xs.*yr + xr.*ys
  double J[15];
  for(int i = 0; i < 15; i++) {
    J[i] = -xs[i] * yr[i] + xr[i] * ys[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < 15; i++) {
    rx[i] = ys[i] / J[i];
    sx[i] = -yr[i] / J[i];
    ry[i] = -xs[i] / J[i];
    sy[i] = xr[i] / J[i];
  }

  // Calculate normals

  // Face 0
  for(int i = 0; i < 5; i++) {
    nx[i] = yr[FMASK[i]];
    ny[i] = -xr[FMASK[i]];
  }
  // Face 1
  for(int i = 0; i < 5; i++) {
    nx[5 + i] = ys[FMASK[5 + i]] - yr[FMASK[5 + i]];
    ny[5 + i] = xr[FMASK[5 + i]] - xs[FMASK[5 + i]];
  }
  // Face 2
  for(int i = 0; i < 5; i++) {
    nx[2 * 5 + i] = -ys[FMASK[2 * 5 + i]];
    ny[2 * 5 + i] = xs[FMASK[2 * 5 + i]];
  }

  // Normalise
  for(int i = 0; i < 3 * 5; i++) {
    double sJ = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    nx[i] = nx[i] / sJ;
    ny[i] = ny[i] / sJ;
    fscale[i] = sJ / J[FMASK[i]];
  }
}
