inline void neighbour_zero(double *exQ) {
  for(int i = 0; i < 4 * 3 * NUM_FACE_PTS; i++) {
    exQ[i] = 0.0;
  }
}
