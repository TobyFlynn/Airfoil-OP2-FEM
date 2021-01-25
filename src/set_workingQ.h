inline void set_workingQ(const double *dt, const int *stage, const double *q,
                         const double *k1, const double *k2, double *workingQ) {
  // for(int i = 0; i < 4 * 15; i++) {
  //   workingQ[i] = q[i] + (*dt) * (*rk_frac) * rk[i];
  // }
  if(stage == 0) {
    for(int i = 0; i < 4 * 15; i++) {
      workingQ[i] = q[i] + (*dt) * k1[i];
    }
  } else {
    for(int i = 0; i < 4 * 15; i++) {
      workingQ[i] = q[i] + (*dt) * (k1[i] / 4.0 + k2[i] / 4.0);
    }
  }
}
