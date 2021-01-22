#include <cblas.h>
#include <algorithm>
#include <cmath>

using namespace std;

inline void advec_rhs(const double *q, const double *exteriorQ,
                      const double *rx, const double *ry, const double *sx,
                      const double *sy, const double *fscale, const double *nx,
                      const double *ny, double *qRHS) {
  // Calculate fluxes
  double F[NUM_SOLUTION_PTS];
  double G[NUM_SOLUTION_PTS];
  for(int i = 0; i < NUM_SOLUTION_PTS; i++) {
    // Calculate flux
    F[i] = velX * q[i];
    G[i] = velY * q[i];
  }

  // Compute weak derivatives
  double dFdr[NUM_SOLUTION_PTS];
  double dFds[NUM_SOLUTION_PTS];
  double dGdr[NUM_SOLUTION_PTS];
  double dGds[NUM_SOLUTION_PTS];

  cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Drw, NUM_SOLUTION_PTS, F, 1, 0.0, dFdr, 1);
  cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Dsw, NUM_SOLUTION_PTS, F, 1, 0.0, dFds, 1);
  cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Drw, NUM_SOLUTION_PTS, G, 1, 0.0, dGdr, 1);
  cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Dsw, NUM_SOLUTION_PTS, G, 1, 0.0, dGds, 1);

  for(int i = 0; i < NUM_SOLUTION_PTS; i++) {
    qRHS[i] = rx[i] * dFdr[i] + sx[i] * dFds[i] + ry[i] * dGdr[i] + sy[i] * dGds[i];
  }

  // Compute primative variables and flux functions for interior points on edges
  double mQ[3 * NUM_FACE_PTS];
  double mF[3 * NUM_FACE_PTS];
  double mG[3 * NUM_FACE_PTS];

  for(int i = 0; i < 3 * NUM_FACE_PTS; i++) {
    int ind = FMASK[i];
    mQ[i] = q[ind];

    // Calculate flux
    mF[i] = velX * q[ind];
    mG[i] = velY * q[ind];
  }

  // Compute primative variables and flux functions for exterior points on edges
  double pF[3 * NUM_FACE_PTS];
  double pG[3 * NUM_FACE_PTS];
  for(int i = 0; i < 3 * NUM_FACE_PTS; i++) {
    // Calculate flux
    pF[i] = velX * exteriorQ[i];
    pG[i] = velY * exteriorQ[i];
  }

  // Lift fluxes
  double nflux[3 * NUM_FACE_PTS];
  double a = sqrt(velX * velX + velY * velY);
  for(int j = 0; j < 3 * NUM_FACE_PTS; j++) {
    nflux[j] = nx[j] * (mF[j] + pF[j]) + ny[j] * (mG[j] + pG[j]) + a * (mQ[j] - exteriorQ[j]);
    nflux[j] *= 0.5 * fscale[j];
  }

  cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, -1.0, LIFT, 15, nflux, 1, 1.0, qRHS, 1);
}
