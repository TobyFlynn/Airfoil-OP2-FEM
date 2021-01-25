inline void flux_zero(double *flux) {
  for(int i = 0; i < 4 * 3 * 5; i++) {
    flux[i] = 0.0;
  }
}
