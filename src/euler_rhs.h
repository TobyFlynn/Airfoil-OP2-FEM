#include <cblas.h>
#include <cmath>

#include "fluxes.h"

inline void euler_rhs(const double *q, double *flux,
                      const double *rx, const double *ry, const double *sx,
                      const double *sy, double *qRHS) {
  double F[4 * 15];
  double G[4 * 15];
  for(int i = 0; i < 15; i++) {
    double rho, u, v, p;
    euler_flux(&q[i * 4], &F[i * 4], &G[i * 4], &rho, &u, &v, &p);
  }

  // Compute weak derivatives
  for(int i = 0; i < 4; i++) {
    double dFdr[15];
    double dFds[15];
    double dGdr[15];
    double dGds[15];

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Drw, 15, &F[i], 4, 0.0, dFdr, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dsw, 15, &F[i], 4, 0.0, dFds, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Drw, 15, &G[i], 4, 0.0, dGdr, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dsw, 15, &G[i], 4, 0.0, dGds, 1);

    for(int j = 0; j < 15; j++) {
      qRHS[i + j * 4] = (rx[j] * dFdr[j] + sx[j] * dFds[j]) + (ry[j] * dGdr[j] + sy[j] * dGds[j]);
    }
  }

  for(int i = 0; i < 4; i++) {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, -1.0, LIFT, 15, flux + i, 4, 1.0, qRHS + i, 4);
  }

  for(int i = 0; i < 4 * 3 * 5; i++) {
    flux[i] = 0.0;
  }
}
