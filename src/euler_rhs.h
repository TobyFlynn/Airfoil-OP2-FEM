#include <cblas.h>
#include <algorithm>
#include <cmath>

using namespace std;

inline void euler_rhs(const double *q, const double *exteriorQ,
                      const double *rx, const double *ry, const double *sx,
                      const double *sy, const double *J, const double *sJ,
                      const double *nx, const double *ny, double *qRHS) {
  double F[4 * NUM_SOLUTION_PTS];
  double G[4 * NUM_SOLUTION_PTS];
  for(int i = 0; i < NUM_SOLUTION_PTS; i++) {
    double rho  = q[i * 4 + 0];
    double rhou = q[i * 4 + 1];
    double rhov = q[i * 4 + 2];
    double ener = q[i * 4 + 3];
    // Calculate primative variables
    double u = rhou / rho;
    double v = rhov / rho;
    double p = (gam - 1) * (ener - 0.5 * (rhou * u + rhov * v));
    // Calculate flux
    F[i * 4 + 0] = rhou;
    F[i * 4 + 1] = rhou * u + p;
    F[i * 4 + 2] = rhov * u;
    F[i * 4 + 3] = u * (ener + p);
    G[i * 4 + 0] = rhov;
    G[i * 4 + 1] = rhou * v;
    G[i * 4 + 2] = rhov * v + p;
    G[i * 4 + 3] = v * (ener + p);
  }

  // Compute weak derivatives
  for(int i = 0; i < 4; i++) {
    double dFdr[NUM_SOLUTION_PTS];
    double dFds[NUM_SOLUTION_PTS];
    double dGdr[NUM_SOLUTION_PTS];
    double dGds[NUM_SOLUTION_PTS];

    cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Drw, NUM_SOLUTION_PTS, &F[i], 4, 0.0, dFdr, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Dsw, NUM_SOLUTION_PTS, &F[i], 4, 0.0, dFds, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Drw, NUM_SOLUTION_PTS, &G[i], 4, 0.0, dGdr, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, NUM_SOLUTION_PTS, NUM_SOLUTION_PTS, 1.0, Dsw, NUM_SOLUTION_PTS, &G[i], 4, 0.0, dGds, 1);

    for(int j = 0; j < NUM_SOLUTION_PTS; j++) {
      qRHS[i + j * 4] = (rx[j] * dFdr[j] + sx[j] * dFds[j]) + (ry[j] * dGdr[j] + sy[j] * dGds[j]);
    }
  }

  // Compute primative variables and flux functions for interior points on edges
  double mQ[4 * 3 * NUM_FACE_PTS];
  double mF[4 * 3 * NUM_FACE_PTS];
  double mG[4 * 3 * NUM_FACE_PTS];
  double mRho[3 * NUM_FACE_PTS];
  double mU[3 * NUM_FACE_PTS];
  double mV[3 * NUM_FACE_PTS];
  double mP[3 * NUM_FACE_PTS];
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    int ind = fmask0[i];
    mQ[i * 4]     = q[ind * 4];
    mQ[i * 4 + 1] = q[ind * 4 + 1];
    mQ[i * 4 + 2] = q[ind * 4 + 2];
    mQ[i * 4 + 3] = q[ind * 4 + 3];

    mRho[i]     = q[ind * 4 + 0];
    double rhou = q[ind * 4 + 1];
    double rhov = q[ind * 4 + 2];
    double ener = q[ind * 4 + 3];
    // Calculate primative variables
    mU[i] = rhou / mRho[i];
    mV[i] = rhov / mRho[i];
    mP[i] = (gam - 1) * (ener - 0.5 * (rhou * mU[i] + rhov * mV[i]));
    // Calculate flux
    mF[i * 4 + 0] = rhou;
    mF[i * 4 + 1] = rhou * mU[i] + mP[i];
    mF[i * 4 + 2] = rhov * mU[i];
    mF[i * 4 + 3] = mU[i] * (ener + mP[i]);
    mG[i * 4 + 0] = rhov;
    mG[i * 4 + 1] = rhou * mV[i];
    mG[i * 4 + 2] = rhov * mV[i] + mP[i];
    mG[i * 4 + 3] = mV[i] * (ener + mP[i]);
  }

  for(int i = 0; i < NUM_FACE_PTS; i++) {
    int ind = fmask1[i];
    int mInd = NUM_FACE_PTS + i;
    int mQInd = 4 * NUM_FACE_PTS;
    mQ[mQInd + i * 4]     = q[ind * 4];
    mQ[mQInd + i * 4 + 1] = q[ind * 4 + 1];
    mQ[mQInd + i * 4 + 2] = q[ind * 4 + 2];
    mQ[mQInd + i * 4 + 3] = q[ind * 4 + 3];

    mRho[mInd]  = q[ind * 4 + 0];
    double rhou = q[ind * 4 + 1];
    double rhov = q[ind * 4 + 2];
    double ener = q[ind * 4 + 3];
    // Calculate primative variables
    mU[mInd] = rhou / mRho[mInd];
    mV[mInd] = rhov / mRho[mInd];
    mP[mInd] = (gam - 1) * (ener - 0.5 * (rhou * mU[mInd] + rhov * mV[mInd]));
    // Calculate flux
    mF[mInd * 4 + 0] = rhou;
    mF[mInd * 4 + 1] = rhou * mU[mInd] + mP[mInd];
    mF[mInd * 4 + 2] = rhov * mU[mInd];
    mF[mInd * 4 + 3] = mU[mInd] * (ener + mP[mInd]);
    mG[mInd * 4 + 0] = rhov;
    mG[mInd * 4 + 1] = rhou * mV[mInd];
    mG[mInd * 4 + 2] = rhov * mV[mInd] + mP[mInd];
    mG[mInd * 4 + 3] = mV[mInd] * (ener + mP[mInd]);
  }

  for(int i = 0; i < NUM_FACE_PTS; i++) {
    int ind = fmask2[i];
    int mInd = 2 * NUM_FACE_PTS + i;
    int mQInd = 2 * 4 * NUM_FACE_PTS;
    mQ[mQInd + i * 4]     = q[ind * 4];
    mQ[mQInd + i * 4 + 1] = q[ind * 4 + 1];
    mQ[mQInd + i * 4 + 2] = q[ind * 4 + 2];
    mQ[mQInd + i * 4 + 3] = q[ind * 4 + 3];

    mRho[mInd]  = q[ind * 4 + 0];
    double rhou = q[ind * 4 + 1];
    double rhov = q[ind * 4 + 2];
    double ener = q[ind * 4 + 3];
    // Calculate primative variables
    mU[mInd] = rhou / mRho[mInd];
    mV[mInd] = rhov / mRho[mInd];
    mP[mInd] = (gam - 1) * (ener - 0.5 * (rhou * mU[mInd] + rhov * mV[mInd]));
    // Calculate flux
    mF[mInd * 4 + 0] = rhou;
    mF[mInd * 4 + 1] = rhou * mU[mInd] + mP[mInd];
    mF[mInd * 4 + 2] = rhov * mU[mInd];
    mF[mInd * 4 + 3] = mU[mInd] * (ener + mP[mInd]);
    mG[mInd * 4 + 0] = rhov;
    mG[mInd * 4 + 1] = rhou * mV[mInd];
    mG[mInd * 4 + 2] = rhov * mV[mInd] + mP[mInd];
    mG[mInd * 4 + 3] = mV[mInd] * (ener + mP[mInd]);
  }

  // Compute primative variables and flux functions for exterior points on edges
  double pF[4 * 3 * NUM_FACE_PTS];
  double pG[4 * 3 * NUM_FACE_PTS];
  double pRho[3 * NUM_FACE_PTS];
  double pU[3 * NUM_FACE_PTS];
  double pV[3 * NUM_FACE_PTS];
  double pP[3 * NUM_FACE_PTS];
  for(int i = 0; i < 3 * NUM_FACE_PTS; i++) {
    pRho[i]     = exteriorQ[i * 4 + 0];
    double rhou = exteriorQ[i * 4 + 1];
    double rhov = exteriorQ[i * 4 + 2];
    double ener = exteriorQ[i * 4 + 3];
    // Calculate primative variables
    pU[i] = rhou / pRho[i];
    pV[i] = rhov / pRho[i];
    pP[i] = (gam - 1) * (ener - 0.5 * (rhou * pU[i] + rhov * pV[i]));
    // Calculate flux
    pF[i * 4 + 0] = rhou;
    pF[i * 4 + 1] = rhou * pU[i] + pP[i];
    pF[i * 4 + 2] = rhov * pU[i];
    pF[i * 4 + 3] = pU[i] * (ener + pP[i]);
    pG[i * 4 + 0] = rhov;
    pG[i * 4 + 1] = rhou * pV[i];
    pG[i * 4 + 2] = rhov * pV[i] + pP[i];
    pG[i * 4 + 3] = pV[i] * (ener + pP[i]);
  }

  // Compute local Lax-Friedrichs flux
  double lamda[3 * NUM_FACE_PTS];
  for(int i = 0; i < 3 * NUM_FACE_PTS; i++) {
    double m = sqrt(mU[i] * mU[i] + mV[i] * mV[i]) + sqrt(abs(gam * mP[i] / mRho[i]));
    double p = sqrt(pU[i] * pU[i] + pV[i] * pV[i]) + sqrt(abs(gam * pP[i] / pRho[i]));
    lamda[i] = max(m, p);
    // Maybe need the following?
    //lamda[i] = max(lamda[i], 1.0);
  }

  // Calculate fscale (might move to init kernel later)
  // And also combine the 3 fmasks to 1 array
  double fscale[3 * NUM_FACE_PTS];
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    fscale[i] = sJ[i] / J[fmask0[i]];
  }
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    fscale[NUM_FACE_PTS + i] = sJ[NUM_FACE_PTS + i] / J[fmask1[i]];
  }
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    fscale[2 * NUM_FACE_PTS + i] = sJ[2 * NUM_FACE_PTS + i] / J[fmask2[i]];
  }

  // Lift fluxes
  for(int i = 0; i < 4; i++) {
    double nflux[3 * NUM_FACE_PTS];
    for(int j = 0; j < 3 * NUM_FACE_PTS; j++) {
      nflux[j] = nx[j] * (pF[i + 4 * j] + mF[i + 4 * j])
                 + ny[j] * (pG[i + 4 * j] + mG[i + 4 * j])
                 + lamda[j] * (mQ[i + j * 4] - exteriorQ[i + j * 4]);
      nflux[j] *= fscale[j] * 0.5;
    }
    double liftMul[3 * NUM_FACE_PTS];
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, lift, 15, nflux, 1, -1.0, qRHS + i, 4);
  }
}
