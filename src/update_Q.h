inline void update_Q(const double *dt, double *q0, double *q1, double *q2,
                     double *q3, const double *rk10, const double *rk11,
                     const double *rk12, const double *rk13, const double *rk20,
                     const double *rk21, const double *rk22, const double *rk23,
                     const double *rk30, const double *rk31, const double *rk32,
                     const double *rk33, double *workingQ0, double *workingQ1,
                     double *workingQ2, double *workingQ3) {
  for(int i = 0; i < 15; i++) {
    q0[i] = q0[i] + (*dt) * (rk10[i]/ 6.0 + rk20[i] / 6.0 + 2.0 * rk30[i] / 3.0);
    workingQ0[i] = q0[i];
    q1[i] = q1[i] + (*dt) * (rk11[i]/ 6.0 + rk21[i] / 6.0 + 2.0 * rk31[i] / 3.0);
    workingQ1[i] = q1[i];
    q2[i] = q2[i] + (*dt) * (rk12[i]/ 6.0 + rk22[i] / 6.0 + 2.0 * rk32[i] / 3.0);
    workingQ2[i] = q2[i];
    q3[i] = q3[i] + (*dt) * (rk13[i]/ 6.0 + rk23[i] / 6.0 + 2.0 * rk33[i] / 3.0);
    workingQ3[i] = q3[i];
  }
}
