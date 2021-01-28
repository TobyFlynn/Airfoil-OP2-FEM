inline void calc_dt(const double *q0, const double *q1, const double *q2,
                    const double *q3, const double *fscale, double *dt1) {
  double dt1_arr[3 * 5];
  for(int i = 0; i < 3 * 5; i++) {
    double rho = q0[FMASK[i]];
    double u = q1[FMASK[i]] / rho;
    double v = q2[FMASK[i]] / rho;
    double p = (gam - 1.0) * (q3[FMASK[i]] - rho * (u * u + v * v) * 0.5);
    double c = sqrt(fabs(gam * p / rho));
    dt1_arr[i] = ((4 + 1) * (4 + 1)) * 0.5 * fscale[FMASK[i]] *(sqrt(u * u + v * v) + c);
  }

  // Find local max
  double max = *dt1;
  for(int i = 0; i < 3 * 5; i++) {
    if(dt1_arr[i] > max) {
      max = dt1_arr[i];
    }
  }
  *dt1 = max;
}
