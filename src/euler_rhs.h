// #include <cblas.h>
// #include <algorithm>
// #include <cmath>

#include "fluxes.h"


inline void euler_rhs(const double *q, double *exteriorQ,
                      const double *rx, const double *ry, const double *sx,
                      const double *sy, const double *fscale, const double *nx,
                      const double *ny, const double *dFdr, const double *dFds,
                      const double *dGdr, const double *dGds, double *qRHS) {
  // Compute weak derivatives
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 15; j++) {
      qRHS[i + j * 4] = (rx[j] * dFdr[i + j * 4] + sx[j] * dFds[i + j * 4]) + (ry[j] * dGdr[i + j * 4] + sy[j] * dGds[i + j * 4]);
    }
  }

  // Compute primative variables and flux functions for interior points on edges
  double mQ[4 * 3 * 5];
  double mF[4 * 3 * 5];
  double mG[4 * 3 * 5];
  double mRho[3 * 5];
  double mU[3 * 5];
  double mV[3 * 5];
  double mP[3 * 5];

  for(int i = 0; i < 3 * 5; i++) {
    int ind = FMASK[i] * 4;
    mQ[i * 4]     = q[ind];
    mQ[i * 4 + 1] = q[ind + 1];
    mQ[i * 4 + 2] = q[ind + 2];
    mQ[i * 4 + 3] = q[ind + 3];

    euler_flux(&mQ[i * 4], &mF[i * 4], &mG[i * 4], &mRho[i], &mU[i], &mV[i], &mP[i]);
  }

  // Compute primative variables and flux functions for exterior points on edges
  double pF[4 * 3 * 5];
  double pG[4 * 3 * 5];
  double pRho[3 * 5];
  double pU[3 * 5];
  double pV[3 * 5];
  double pP[3 * 5];
  for(int i = 0; i < 3 * 5; i++) {
    euler_flux(&exteriorQ[i * 4], &pF[i * 4], &pG[i * 4], &pRho[i], &pU[i], &pV[i], &pP[i]);
  }

  double flux[4 * 3 * 5];
  // lax_friedrichs(flux, nx, ny, fscale, q, exteriorQ);
  roe(flux, nx, ny, fscale, q, exteriorQ);

  // for(int i = 0; i < 4; i++) {
  //   cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, -1.0, LIFT, 15, &flux[i], 4, 1.0, qRHS + i, 4);
  // }

  for(int i = 0; i < 4 * 3 * 5; i++) {
    exteriorQ[i] = 0.0;
  }
}
