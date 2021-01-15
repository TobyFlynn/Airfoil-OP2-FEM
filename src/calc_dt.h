#include <cmath>
#include <iostream>

inline void calc_dt(const double *q, const double *J, const double *sJ,
                    double *dt1) {
  double dt1_arr[3 * NUM_FACE_PTS];
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    double fscale = sJ[i] / J[fmask0[i]];
    double rho = q[fmask0[i] * 4];
    double u = q[fmask0[i] * 4 + 1] / rho;
    double v = q[fmask0[i] * 4 + 2] / rho;
    double p = (gam - 1.0) * (q[fmask0[i] * 4 + 3] - rho * (u * u + v * v) * 0.5);
    double c = sqrt(abs(gam * p / rho));
    dt1_arr[i] = ((ORDER + 1) * (ORDER + 1)) * 0.5 * fscale *(sqrt(u * u + v * v) + c);
  }
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    double fscale = sJ[NUM_FACE_PTS + i] / J[fmask1[i]];
    double rho = q[fmask1[i] * 4];
    double u = q[fmask1[i] * 4 + 1] / rho;
    double v = q[fmask1[i] * 4 + 2] / rho;
    double p = (gam - 1.0) * (q[fmask1[i] * 4 + 3] - rho * (u * u + v * v) * 0.5);
    double c = sqrt(abs(gam * p / rho));
    dt1_arr[NUM_FACE_PTS + i] = ((ORDER + 1) * (ORDER + 1)) * 0.5 * fscale *(sqrt(u * u + v * v) + c);
  }
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    double fscale = sJ[2 * NUM_FACE_PTS + i] / J[fmask2[i]];
    double rho = q[fmask2[i] * 4];
    double u = q[fmask2[i] * 4 + 1] / rho;
    double v = q[fmask2[i] * 4 + 2] / rho;
    double p = (gam - 1.0) * (q[fmask2[i] * 4 + 3] - rho * (u * u + v * v) * 0.5);
    double c = sqrt(abs(gam * p / rho));
    dt1_arr[2 * NUM_FACE_PTS + i] = ((ORDER + 1) * (ORDER + 1)) * 0.5 * fscale *(sqrt(u * u + v * v) + c);
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
