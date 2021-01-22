#include <cmath>

inline void set_ic(const double *x, const double *y, double *q,
                   double *workingQ) {
  double t = 0.0;
  for(int i = 0; i < NUM_SOLUTION_PTS; i++) {
    q[i]        = exp(-40.0 * (x[i] - 0.5) * (x[i] - 0.5)) * exp(-40.0 * (y[i] - 0.5) * (y[i] - 0.5));
    workingQ[i] = q[i];
  }
}
