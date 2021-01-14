inline void update_Q(double *q, const double *rk1, const double *rk2,
                     const double *rk3, const double *rk4, double *workingQ) {
  for(int i = 0; i < 4 * NUM_SOLUTION_PTS; i++) {
    q[i] = q[i] + (1.0/6.0) * dt * (rk1[i] + 2 * rk2[i] + 2 * rk3[i] + rk4[i]);
    workingQ[i] = q[i];
  }
}
