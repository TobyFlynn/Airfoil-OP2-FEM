#include "fluxes.h"

inline void internal_fluxes(const double *q, double *F, double *G) {
  for(int i = 0; i < 15; i++) {
    double rho, u, v, p;
    euler_flux(&q[i * 4], &F[i * 4], &G[i * 4], &rho, &u, &v, &p);
  }
}
