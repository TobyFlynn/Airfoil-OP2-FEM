inline void set_ic(double *q0, double *q1, double *q2, double *q3,
                   double *workingQ0, double *workingQ1, double *workingQ2,
                   double *workingQ3, double *exQ0, double *exQ1, double *exQ2,
                   double *exQ3) {
  for(int i = 0; i < 15; i++) {
    q0[i] = bc_r;
    q1[i] = bc_r * bc_u;
    q2[i] = bc_r * bc_v;
    q3[i] = bc_e;
    workingQ0[i] = q0[i];
    workingQ1[i] = q1[i];
    workingQ2[i] = q2[i];
    workingQ3[i] = q3[i];
    exQ0[i] = 0.0;
    exQ1[i] = 0.0;
    exQ2[i] = 0.0;
    exQ3[i] = 0.0;
  }
}
