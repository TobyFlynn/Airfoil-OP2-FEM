#include <cmath>

// fmask0 = 0, 1, 2, 3, 4
// fmask1 = 4, 8, 11, 13, 14
// fmask2 = 0, 5, 9, 12, 14

inline void calc_normals(const double *xr, const double *yr, const double *xs,
                         const double *ys, double *nx, double *ny, double *sJ) {
  // Face 0
  nx[0] = yr[0];
  nx[1] = yr[1];
  nx[2] = yr[2];
  nx[3] = yr[3];
  nx[4] = yr[4];
  ny[0] = -xr[0];
  ny[1] = -xr[1];
  ny[2] = -xr[2];
  ny[3] = -xr[3];
  ny[4] = -xr[4];
  // Face 1
  nx[5] = ys[4] - yr[4];
  nx[6] = ys[8] - yr[8];
  nx[7] = ys[11] - yr[11];
  nx[8] = ys[13] - yr[13];
  nx[9] = ys[14] - yr[14];
  ny[5] = xr[4] - xs[4];
  ny[6] = xr[8] - xs[8];
  ny[7] = xr[11] - xs[11];
  ny[8] = xr[13] - xs[13];
  ny[9] = xr[14] - xs[14];
  // Face 2
  nx[10] = -ys[0];
  nx[11] = -ys[5];
  nx[12] = -ys[9];
  nx[13] = -ys[12];
  nx[14] = -ys[14];
  ny[10] = xs[0];
  ny[11] = xs[5];
  ny[12] = xs[9];
  ny[13] = xs[12];
  ny[14] = xs[14];

  // Normalise
  for(int i = 0; i < 3 * NUM_FACE_PTS; i++) {
    sJ[i] = std::sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    nx[i] = nx[i] / sJ[i];
    ny[i] = ny[i] / sJ[i];
  }
}
