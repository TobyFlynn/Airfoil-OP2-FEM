#include <cmath>

inline void get_bedge_q(const int *bedgeNum, const double *n0, const double *n1,
                        const double *nodeX, const double *nodeY,
                        const double *x, const double *y, double *exteriorQ) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = NUM_FACE_PTS;
  else if(*bedgeNum == 2) exInd = 2 * NUM_FACE_PTS;

  for(int i = 0; i < NUM_FACE_PTS; i++) {
    exteriorQ[exInd + i] += 4.5e-5;
  }
}
