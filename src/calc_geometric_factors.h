#include <cblas.h>

inline void calc_geometric_factors(const double *x, const double *y,
                                   double *xr, double *yr, double *xs,
                                   double *ys, double *J, double *rx,
                                   double *ry, double *sx, double *sy) {
  // xr = Dr * x
  // cblas_dcopy(NUM_SOLUTION_PTS, zeros, 1, xr, 1); - Apparently not needed if beta set to 0.0
  cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Dr, NUM_SOLUTION_PTS, x, 1, 0.0, xr, 1);

  // xs = Ds * x
  cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Ds, NUM_SOLUTION_PTS, x, 1, 0.0, xs, 1);

  // yr = Dr * y
  cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Dr, NUM_SOLUTION_PTS, y, 1, 0.0, yr, 1);

  // ys = Ds * y
  cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Ds, NUM_SOLUTION_PTS, y, 1, 0.0, ys, 1);

  // J = -xs.*yr + xr.*ys
  for(int i = 0; i < NUM_SOLUTION_PTS; i++) {
    J[i] = -xs[i] * yr[i] + xr[i] * ys[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < NUM_SOLUTION_PTS; i++) {
    rx[i] = ys[i] / J[i];
    sx[i] = -yr[i] / J[i];
    ry[i] = -xs[i] / J[i];
    sy[i] = xr[i] / J[i];
  }
}
