inline void set_workingQ(const double *dt, const int *stage, const double *q0,
                         const double *q1, const double *q2, const double *q3,
                         const double *k10, const double *k11,
                         const double *k12, const double *k13,
                         const double *k20, const double *k21,
                         const double *k22, const double *k23,
                         double *workingQ0, double *workingQ1,
                         double *workingQ2, double *workingQ3) {
  if(*stage == 0) {
    for(int i = 0; i < 15; i++) {
      workingQ0[i] = q0[i] + (*dt) * k10[i];
      workingQ1[i] = q1[i] + (*dt) * k11[i];
      workingQ2[i] = q2[i] + (*dt) * k12[i];
      workingQ3[i] = q3[i] + (*dt) * k13[i];
    }
  } else {
    for(int i = 0; i < 15; i++) {
      workingQ0[i] = q0[i] + (*dt) * (k10[i] / 4.0 + k20[i] / 4.0);
      workingQ1[i] = q1[i] + (*dt) * (k11[i] / 4.0 + k21[i] / 4.0);
      workingQ2[i] = q2[i] + (*dt) * (k12[i] / 4.0 + k22[i] / 4.0);
      workingQ3[i] = q3[i] + (*dt) * (k13[i] / 4.0 + k23[i] / 4.0);
    }
  }
}
