inline void backwards_euler_update_Q(const double *dt, const double *q0,
                     const double *q1, const double *q2, const double *q3,
                     double *rhs0, double *rhs1, double *rhs2, double *rhs3) {
  for(int i = 0; i < 15; i++) {
    rhs0[i] = q0[i] - (*dt) * rhs0[i];
    rhs1[i] = q1[i] - (*dt) * rhs1[i];
    rhs2[i] = q2[i] - (*dt) * rhs2[i];
    rhs3[i] = q3[i] - (*dt) * rhs3[i];
  }
}
