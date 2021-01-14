inline void set_ic(double *q) {
  for(int i = 0; i < NUM_SOLUTION_PTS; i++) {
    q[i * 4]     = bc_r;
    q[i * 4 + 1] = bc_r * bc_u;
    q[i * 4 + 2] = 0.0;
    q[i * 4 + 3] = bc_e;
  }
}
