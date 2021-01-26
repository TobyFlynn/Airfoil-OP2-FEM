inline void update_Q(const double *dt, double *q, const double *rk1,
                     const double *rk2, const double *rk3, //const double *rk4,
                     double *workingQ) {
  for(int i = 0; i < 4 * 15; i++) {
    q[i] = q[i] + (*dt) * (rk1[i]/ 6.0 + rk2[i] / 6.0 + 2.0 * rk3[i] / 3.0);
    workingQ[i] = q[i];
  }
}
