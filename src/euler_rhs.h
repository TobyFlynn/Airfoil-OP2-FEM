#include "fluxes.h"


inline void euler_rhs(const double *q0, const double *q1, const double *q2,
                      const double *q3, double *exteriorQ0, double *exteriorQ1,
                      double *exteriorQ2, double *exteriorQ3, const double *rx,
                      const double *ry, const double *sx, const double *sy,
                      const double *fscale, const double *nx, const double *ny,
                      const double *dFdr0, const double *dFdr1,
                      const double *dFdr2, const double *dFdr3,
                      const double *dFds0, const double *dFds1,
                      const double *dFds2, const double *dFds3,
                      const double *dGdr0, const double *dGdr1,
                      const double *dGdr2, const double *dGdr3,
                      const double *dGds0, const double *dGds1,
                      const double *dGds2, const double *dGds3, double *flux0,
                      double *flux1, double *flux2, double *flux3,
                      double *qRHS0, double *qRHS1, double *qRHS2,
                      double *qRHS3) {
  // Compute weak derivatives
  for(int j = 0; j < 15; j++) {
    qRHS0[j] = (rx[j] * dFdr0[j] + sx[j] * dFds0[j]) + (ry[j] * dGdr0[j] + sy[j] * dGds0[j]);
  }

  for(int j = 0; j < 15; j++) {
    qRHS1[j] = (rx[j] * dFdr1[j] + sx[j] * dFds1[j]) + (ry[j] * dGdr1[j] + sy[j] * dGds1[j]);
  }

  for(int j = 0; j < 15; j++) {
    qRHS2[j] = (rx[j] * dFdr2[j] + sx[j] * dFds2[j]) + (ry[j] * dGdr2[j] + sy[j] * dGds2[j]);
  }

  for(int j = 0; j < 15; j++) {
    qRHS3[j] = (rx[j] * dFdr3[j] + sx[j] * dFds3[j]) + (ry[j] * dGdr3[j] + sy[j] * dGds3[j]);
  }

  // lax_friedrichs(flux0, flux1, flux2, flux3, nx, ny, fscale, q0, q1, q2, q3, exteriorQ0,
  //     exteriorQ1, exteriorQ2, exteriorQ3);
  roe(flux0, flux1, flux2, flux3, nx, ny, fscale, q0, q1, q2, q3, exteriorQ0,
      exteriorQ1, exteriorQ2, exteriorQ3);

  for(int i = 0; i < 3 * 5; i++) {
    exteriorQ0[i] = 0.0;
    exteriorQ1[i] = 0.0;
    exteriorQ2[i] = 0.0;
    exteriorQ3[i] = 0.0;
  }
}
