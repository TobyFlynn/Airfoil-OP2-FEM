#include <cblas.h>

inline void filterQ(double *rhsQ) {
  for(int i = 0; i < 4; i++) {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, filter, NUM_SOLUTION_PTS, rhsQ + i, 4, 0.0, rhsQ, 4);
  }
}
