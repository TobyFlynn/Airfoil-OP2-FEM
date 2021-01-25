inline void neighbour_zero(double *exQ) {
  for(int i = 0; i < 4 * 3 * 5; i++) {
    exQ[i] = 0.0;
  }
}
