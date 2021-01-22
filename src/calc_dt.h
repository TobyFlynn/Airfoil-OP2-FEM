#include <cmath>
#include <iostream>

inline void calc_dt(const double *q, const double *fscale, double *dt1) {
  double dt1_arr[3 * NUM_FACE_PTS];
  for(int i = 0; i < 3 * NUM_FACE_PTS; i++) {
    double rho = q[FMASK[i] * 4];
    double u = q[FMASK[i] * 4 + 1] / rho;
    double v = q[FMASK[i] * 4 + 2] / rho;
    double p = (gam - 1.0) * (q[FMASK[i] * 4 + 3] - rho * (u * u + v * v) * 0.5);
    double c = sqrt(abs(gam * p / rho));
    dt1_arr[i] = ((ORDER + 1) * (ORDER + 1)) * 0.5 * fscale[FMASK[i]] *(sqrt(u * u + v * v) + c);
  }

  // Find local max
  double max = dt1_arr[0];
  for(int i = 1; i < 3 * NUM_FACE_PTS; i++) {
    if(dt1_arr[i] > max) {
      max = dt1_arr[i];
    }
  }
  *dt1 = max;
}
