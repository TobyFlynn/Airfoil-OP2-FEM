#include <cblas.h>

inline void calc_solution_pt_coords(const double *n0, const double *n1, const double *n2,
                                    double *x, double *y) {
  // x = 0.5(-n0_x (r + s) + n1_x(1 + r) + n2_x(1 + s))
  // x = 0.5*n1_x * (1 + r) + 0.5*n2_x * (1 + s) - 0.5*n0_x * (r + s)

  // x <- 0.5*n1_x * (1 + r)
  cblas_dcopy(NUM_SOLUTION_PTS, ones, 1, x, 1);
  cblas_daxpy(NUM_SOLUTION_PTS, 1.0, solution_pts_r, 1, x, 1);
  cblas_dscal(NUM_SOLUTION_PTS, 0.5 * n1[0], x, 1);
  // temp <- 1 + s
  double temp[NUM_SOLUTION_PTS];
  cblas_dcopy(NUM_SOLUTION_PTS, ones, 1, temp, 1);
  cblas_daxpy(NUM_SOLUTION_PTS, 1.0, solution_pts_s, 1, temp, 1);
  // x <- 0.5*n2_x * temp + x (i.e. 0.5*n1_x * (1 + r) + 0.5*n2_x * (1 + s))
  cblas_daxpy(NUM_SOLUTION_PTS, 0.5 * n2[0], temp, 1, x, 1);
  // temp <- r + s
  cblas_dcopy(NUM_SOLUTION_PTS, solution_pts_s, 1, temp, 1);
  cblas_daxpy(NUM_SOLUTION_PTS, 1.0, solution_pts_r, 1, temp, 1);
  // x <- -0.5 * n0_x * temp + x (i.e. final target)
  cblas_daxpy(NUM_SOLUTION_PTS, -0.5 * n0[0], temp, 1, x, 1);

  // y = 0.5(-n0_y (r + s) + n1_y(1 + r) + n2_y(1 + s))
  // y = 0.5*n1_y * (1 + r) + 0.5*n2_y * (1 + s) - 0.5*n0_y * (r + s)

  // y <- 0.5*n1_y * (1 + r)
  cblas_dcopy(NUM_SOLUTION_PTS, ones, 1, y, 1);
  cblas_daxpy(NUM_SOLUTION_PTS, 1.0, solution_pts_r, 1, y, 1);
  cblas_dscal(NUM_SOLUTION_PTS, 0.5 * n1[1], y, 1);
  // temp <- 1 + s
  cblas_dcopy(NUM_SOLUTION_PTS, ones, 1, temp, 1);
  cblas_daxpy(NUM_SOLUTION_PTS, 1.0, solution_pts_s, 1, temp, 1);
  // y <- 0.5*n2_y * temp + y (i.e. 0.5*n1_y * (1 + r) + 0.5*n2_y * (1 + s))
  cblas_daxpy(NUM_SOLUTION_PTS, 0.5 * n2[1], temp, 1, y, 1);
  // temp <- r + s
  cblas_dcopy(NUM_SOLUTION_PTS, solution_pts_s, 1, temp, 1);
  cblas_daxpy(NUM_SOLUTION_PTS, 1.0, solution_pts_r, 1, temp, 1);
  // x <- -0.5 * n0_y * temp + x (i.e. final target)
  cblas_daxpy(NUM_SOLUTION_PTS, -0.5 * n0[1], temp, 1, y, 1);
}
