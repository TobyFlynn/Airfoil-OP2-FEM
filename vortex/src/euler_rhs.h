#include <cblas.h>
#include <algorithm>
#include <cmath>

using namespace std;

inline void euler_flux(const double *q, double *F, double *G, double *rho,
                       double *u, double *v, double *p) {
  *rho = q[0];
  *u = q[1] / (*rho);
  *v = q[2] / (*rho);
  *p = (gam - 1) * (q[3] - 0.5 * (q[1] * (*u) + q[2] * (*v)));

  F[0] = q[1];
  F[1] = q[1] * (*u) + (*p);
  F[2] = q[2] * (*u);
  F[3] = (*u) * (q[3] + (*p));

  G[0] = q[2];
  G[1] = q[1] * (*v);
  G[2] = q[2] * (*v) + (*p);
  G[3] = (*v) * (q[3] + (*p));
}

inline void euler_rhs(const double *q, const double *exteriorQ,
                      const double *rx, const double *ry, const double *sx,
                      const double *sy, const double *fscale, const double *nx,
                      const double *ny, double *qRHS) {
  double F[4 * NUM_SOLUTION_PTS];
  double G[4 * NUM_SOLUTION_PTS];
  for(int i = 0; i < NUM_SOLUTION_PTS; i++) {
    double rho, u, v, p;
    euler_flux(&q[i * 4], &F[i * 4], &G[i * 4], &rho, &u, &v, &p);
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

  for(int i = 0; i < 3 * NUM_FACE_PTS; i++) {
    int ind = FMASK[i] * 4;
    mQ[i * 4]     = q[ind];
    mQ[i * 4 + 1] = q[ind + 1];
    mQ[i * 4 + 2] = q[ind + 2];
    mQ[i * 4 + 3] = q[ind + 3];

    euler_flux(&mQ[i * 4], &mF[i * 4], &mG[i * 4], &mRho[i], &mU[i], &mV[i], &mP[i]);
  }

  // Compute primative variables and flux functions for exterior points on edges
  double pF[4 * 3 * NUM_FACE_PTS];
  double pG[4 * 3 * NUM_FACE_PTS];
  double pRho[3 * NUM_FACE_PTS];
  double pU[3 * NUM_FACE_PTS];
  double pV[3 * NUM_FACE_PTS];
  double pP[3 * NUM_FACE_PTS];
  for(int i = 0; i < 3 * NUM_FACE_PTS; i++) {
    euler_flux(&exteriorQ[i * 4], &pF[i * 4], &pG[i * 4], &pRho[i], &pU[i], &pV[i], &pP[i]);
  }

  // // Compute local Lax-Friedrichs flux
  // // Max lamda for each face
  // double maxL;
  // double maxLamda[3 * NUM_FACE_PTS];
  // for(int i = 0; i < NUM_FACE_PTS; i++) {
  //   double m = sqrt(mU[i] * mU[i] + mV[i] * mV[i]) + sqrt(abs(gam * mP[i] / mRho[i]));
  //   double p = sqrt(pU[i] * pU[i] + pV[i] * pV[i]) + sqrt(abs(gam * pP[i] / pRho[i]));
  //   double lamda = max(m, p);
  //   if(i == 0 || lamda > maxL) {
  //     maxL = lamda;
  //   }
  // }
  //
  // for(int i = 0; i < NUM_FACE_PTS; i++) {
  //   maxLamda[i] = maxL;
  // }
  //
  // for(int i = NUM_FACE_PTS; i < 2 * NUM_FACE_PTS; i++) {
  //   double m = sqrt(mU[i] * mU[i] + mV[i] * mV[i]) + sqrt(abs(gam * mP[i] / mRho[i]));
  //   double p = sqrt(pU[i] * pU[i] + pV[i] * pV[i]) + sqrt(abs(gam * pP[i] / pRho[i]));
  //   double lamda = max(m, p);
  //   if(i == NUM_FACE_PTS || lamda > maxL) {
  //     maxL = lamda;
  //   }
  // }
  //
  // for(int i = NUM_FACE_PTS; i < 2 * NUM_FACE_PTS; i++) {
  //   maxLamda[i] = maxL;
  // }
  //
  // for(int i = 2 * NUM_FACE_PTS; i < 3 * NUM_FACE_PTS; i++) {
  //   double m = sqrt(mU[i] * mU[i] + mV[i] * mV[i]) + sqrt(abs(gam * mP[i] / mRho[i]));
  //   double p = sqrt(pU[i] * pU[i] + pV[i] * pV[i]) + sqrt(abs(gam * pP[i] / pRho[i]));
  //   double lamda = max(m, p);
  //   if(i == 2 * NUM_FACE_PTS || lamda > maxL) {
  //     maxL = lamda;
  //   }
  // }
  //
  // for(int i = 2 * NUM_FACE_PTS; i < 3 * NUM_FACE_PTS; i++) {
  //   maxLamda[i] = maxL;
  // }
  //
  // // double maxLamda[3 * NUM_FACE_PTS];
  // // for(int i = 0; i < NUM_FACE_PTS; i++) {
  // //   double m = sqrt(mU[i] * mU[i] + mV[i] * mV[i]) + sqrt(abs(gam * mP[i] / mRho[i]));
  // //   double p = sqrt(pU[i] * pU[i] + pV[i] * pV[i]) + sqrt(abs(gam * pP[i] / pRho[i]));
  // //   maxLamda[i] = max(m, p);
  // // }
  //
  // // Lift fluxes
  // for(int i = 0; i < 4; i++) {
  //   double nflux[3 * NUM_FACE_PTS];
  //   for(int j = 0; j < 3 * NUM_FACE_PTS; j++) {
  //     nflux[j] = nx[j] * (pF[i + 4 * j] + mF[i + 4 * j])
  //                + ny[j] * (pG[i + 4 * j] + mG[i + 4 * j])
  //                + maxLamda[j] * (mQ[i + j * 4] - exteriorQ[i + j * 4]);
  //     nflux[j] *= fscale[j] * 0.5;
  //   }
  //   cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, -1.0, LIFT, 15, nflux, 1, 1.0, qRHS + i, 4);
  // }
  // Roe Flux
  double flux[4 * 3 * NUM_FACE_PTS];
  for(int i = 0; i < 3 * NUM_FACE_PTS; i++) {
    double mRoeQ[4];
    mRoeQ[0] = mQ[i * 4];
    mRoeQ[1] = nx[i] * mQ[i * 4 + 1] + ny[i] * mQ[i * 4 + 2];
    mRoeQ[2] = -ny[i] * mQ[i * 4 + 1] + nx[i] * mQ[i * 4 + 2];
    // mRoeQ[1] = mQ[i * 4 + 1];
    // mRoeQ[2] = mQ[i * 4 + 2];
    mRoeQ[3] = mQ[i * 4 + 3];

    double pRoeQ[4];
    pRoeQ[0] = exteriorQ[i * 4];
    pRoeQ[1] = nx[i] * exteriorQ[i * 4 + 1] + ny[i] * exteriorQ[i * 4 + 2];
    pRoeQ[2] = -ny[i] * exteriorQ[i * 4 + 1] + nx[i] * exteriorQ[i * 4 + 2];
    // pRoeQ[1] = exteriorQ[i * 4 + 1];
    // pRoeQ[2] = exteriorQ[i * 4 + 2];
    pRoeQ[3] = exteriorQ[i * 4 + 3];

    double fxMQ[4];
    double fyMQ[4];
    double mRoeRho, mRoeU, mRoeV, mRoeP;
    euler_flux(&mRoeQ[0], &fxMQ[0], &fyMQ[0], &mRoeRho, &mRoeU, &mRoeV, &mRoeP);

    double fxPQ[4];
    double fyPQ[4];
    double pRoeRho, pRoeU, pRoeV, pRoeP;
    euler_flux(&pRoeQ[0], &fxPQ[0], &fyPQ[0], &pRoeRho, &pRoeU, &pRoeV, &pRoeP);

    double mH = (mRoeQ[3] + mRoeP) / mRoeRho;
    double pH = (pRoeQ[3] + pRoeP) / pRoeRho;

    double mRoeRhoSqrt = sqrt(mRoeRho);
    double pRoeRhoSqrt = sqrt(pRoeRho);

    double rho = mRoeRhoSqrt * pRoeRhoSqrt;
    double u = (mRoeRhoSqrt * mRoeU + pRoeRhoSqrt * pRoeU) / (mRoeRhoSqrt + pRoeRhoSqrt);
    double v = (mRoeRhoSqrt * mRoeV + pRoeRhoSqrt * pRoeV) / (mRoeRhoSqrt + pRoeRhoSqrt);
    double H = (mRoeRhoSqrt * mH + pRoeRhoSqrt * pH) / (mRoeRhoSqrt + pRoeRhoSqrt);

    double c2 = (gam - 1) * (H - 0.5 * (u * u + v * v));
    double c = sqrt(c2);

    double dW1 = -0.5 * rho * (pRoeU - mRoeU) / c + 0.5 * (pRoeP - mRoeP) / c2;
    double dW2 = (pRoeRho - mRoeRho) - (pRoeP - mRoeP) / c2;
    double dW3 = rho * (pRoeV - mRoeV);
    double dW4 = 0.5 * rho * (pRoeU - mRoeU) / c + 0.5 * (pRoeP - mRoeP) / c2;

    dW1 = abs(u - c) * dW1;
    dW2 = abs(u) * dW2;
    dW3 = abs(u) * dW3;
    dW4 = abs(u + c) * dW4;

    double fx[4];
    fx[0] = (fxMQ[0] + fxPQ[0]) / 2.0;
    fx[1] = (fxMQ[1] + fxPQ[1]) / 2.0;
    fx[2] = (fxMQ[2] + fxPQ[2]) / 2.0;
    fx[3] = (fxMQ[3] + fxPQ[3]) / 2.0;

    fx[0] = fx[0] - (dW1 + dW2 + dW4) / 2.0;
    fx[1] = fx[1] - (dW1 * (u - c) + dW2 * u + dW4 * (u + c)) / 2.0;
    fx[2] = fx[2] - (dW1 * v + dW2 * v + dW3 + dW4 * v) / 2.0;
    fx[3] = fx[3] - (dW1 * (H - u * c) + dW2 * (u * u + v * v) + dW3 * v + dW4 * (H + u * c)) / 2.0;

    flux[i * 4]     = fx[0];
    flux[i * 4 + 1] = nx[i] * fx[1] - ny[i] * fx[2];
    flux[i * 4 + 2] = ny[i] * fx[1] + nx[i] * fx[2];
    flux[i * 4 + 3] = fx[3];

    flux[i * 4]     *= fscale[i];
    flux[i * 4 + 1] *= fscale[i];
    flux[i * 4 + 2] *= fscale[i];
    flux[i * 4 + 3] *= fscale[i];
  }

  for(int i = 0; i < 4; i++) {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, -1.0, LIFT, 15, &flux[i], 4, 1.0, qRHS + i, 4);
  }
}
