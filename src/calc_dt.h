inline void calc_dt(const double *q, const double *fscale, double *dt1) {
  double dt1_arr[3 * 5];
  for(int i = 0; i < 3 * 5; i++) {
    double rho = q[FMASK[i] * 4];
    double u = q[FMASK[i] * 4 + 1] / rho;
    double v = q[FMASK[i] * 4 + 2] / rho;
    double p = (gam - 1.0) * (q[FMASK[i] * 4 + 3] - rho * (u * u + v * v) * 0.5);
    double c = sqrt(abs(gam * p / rho));
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
