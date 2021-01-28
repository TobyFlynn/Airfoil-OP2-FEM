#include "fluxes.h"

inline void internal_fluxes(const double *q0, const double *q1,
                            const double *q2, const double *q3, double *f0,
                            double *f1, double *f2, double *f3, double *g0,
                            double *g1, double *g2, double *g3) {
  for(int i = 0; i < 15; i++) {
    double rho, u, v, p;
    euler_flux(q0[i], q1[i], q2[i], q3[i], &f0[i], &f1[i], &f2[i], &f3[i],
               &g0[i], &g1[i], &g2[i], &g3[i], &rho, &u, &v, &p);
  }
}
