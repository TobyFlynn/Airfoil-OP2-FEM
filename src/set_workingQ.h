inline void set_workingQ(const double *dt, const double *q, const double *rk,
                         const double *rk_frac, double *workingQ) {
  for(int i = 0; i < 4 * NUM_SOLUTION_PTS; i++) {
    workingQ[i] = q[i] + (*dt) * (*rk_frac) * rk[i];
  }
}
