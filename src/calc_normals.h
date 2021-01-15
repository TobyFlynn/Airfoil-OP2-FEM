#include <cmath>

// fmask0 = 0, 1, 2, 3, 4
// fmask1 = 4, 8, 11, 13, 14
// fmask2 = 0, 5, 9, 12, 14

inline void calc_normals(const double *xr, const double *yr, const double *xs,
                         const double *ys, double *nx, double *ny, double *sJ) {
  // Face 0
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    nx[i] = yr[fmask0[i]];
    ny[i] = -xr[fmask0[i]];
  }
  // Face 1
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    nx[NUM_FACE_PTS + i] = ys[fmask1[i]] - yr[fmask1[i]];
    ny[NUM_FACE_PTS + i] = xr[fmask1[i]] - xs[fmask1[i]];
  }
  // Face 2
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    nx[2 * NUM_FACE_PTS + i] = -ys[fmask2[i]];
    ny[2 * NUM_FACE_PTS + i] = xs[fmask2[i]];
  }

  // Normalise
  for(int i = 0; i < 3 * NUM_FACE_PTS; i++) {
    sJ[i] = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    nx[i] = nx[i] / sJ[i];
    ny[i] = ny[i] / sJ[i];
  }
}
