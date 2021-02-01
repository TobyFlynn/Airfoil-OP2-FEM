#ifndef __AIRFOIL_FLUXES_H
#define __AIRFOIL_FLUXES_H

using namespace std;

__device__ void euler_flux(const double q0, const double q1, const double q2,
                           const double q3, double *f0, double *f1, double *f2,
                           double *f3, double *g0, double *g1, double *g2,
                           double *g3, double *rho, double *u, double *v,
                           double *p) {
  double gam = 1.4;
  *rho = q0;
  *u = q1 / (*rho);
  *v = q2 / (*rho);
  *p = (gam - 1.0) * (q3 - 0.5 * (q1 * (*u) + q2 * (*v)));

  *f0 = q1;
  *f1 = q1 * (*u) + (*p);
  *f2 = q2 * (*u);
  *f3 = (*u) * (q3 + (*p));

  *g0 = q2;
  *g1 = q1 * (*v);
  *g2 = q2 * (*v) + (*p);
  *g3 = (*v) * (q3 + (*p));
}

__device__ void lax_friedrichs(double *flux0, double *flux1, double *flux2, double *flux3,
                    const double *nx, const double *ny, const double *fscale,
                    const double *q0, const double *q1, const double *q2,
                    const double *q3, const double *pQ0, const double *pQ1,
                    const double *pQ2, const double *pQ3) {
  double gam = 1.4;
  int FMASK[15] = {0, 1, 2, 3, 4, 4, 8, 11, 13, 14, 14, 12, 9, 5, 0};
  // Compute primative variables and flux functions for interior points on edges
  double mQ0[3 * 5]; double mQ1[3 * 5]; double mQ2[3 * 5]; double mQ3[3 * 5];
  double mF0[3 * 5]; double mF1[3 * 5]; double mF2[3 * 5]; double mF3[3 * 5];
  double mG0[3 * 5]; double mG1[3 * 5]; double mG2[3 * 5]; double mG3[3 * 5];
  double mRho[3 * 5];
  double mU[3 * 5];
  double mV[3 * 5];
  double mP[3 * 5];

  for(int i = 0; i < 3 * 5; i++) {
    int ind = FMASK[i];
    mQ0[i] = q0[ind];
    mQ1[i] = q1[ind];
    mQ2[i] = q2[ind];
    mQ3[i] = q3[ind];

    euler_flux(mQ0[i], mQ1[i], mQ2[i], mQ3[i], &mF0[i], &mF1[i], &mF2[i],
               &mF3[i], &mG0[i], &mG1[i], &mG2[i], &mG3[i], &mRho[i], &mU[i],
               &mV[i], &mP[i]);
  }

  // Compute primative variables and flux functions for exterior points on edges
  double pF0[3 * 5]; double pF1[3 * 5]; double pF2[3 * 5]; double pF3[3 * 5];
  double pG0[3 * 5]; double pG1[3 * 5]; double pG2[3 * 5]; double pG3[3 * 5];
  double pRho[3 * 5];
  double pU[3 * 5];
  double pV[3 * 5];
  double pP[3 * 5];
  for(int i = 0; i < 3 * 5; i++) {
    euler_flux(pQ0[i], pQ1[i], pQ2[i], pQ3[i], &pF0[i], &pF1[i], &pF2[i],
               &pF3[i], &pG0[i], &pG1[i], &pG2[i], &pG3[i], &pRho[i], &pU[i],
               &pV[i], &pP[i]);
  }

  // Compute local Lax-Friedrichs flux
  // Max lamda for each face
  double maxL;
  double maxLamda[3 * 5];
  for(int i = 0; i < 5; i++) {
    double m = sqrt(mU[i] * mU[i] + mV[i] * mV[i]) + sqrt(fabs(gam * mP[i] / mRho[i]));
    double p = sqrt(pU[i] * pU[i] + pV[i] * pV[i]) + sqrt(fabs(gam * pP[i] / pRho[i]));
    double lamda = max(m, p);
    if(i == 0 || lamda > maxL) {
      maxL = lamda;
    }
  }

  for(int i = 0; i < 5; i++) {
    maxLamda[i] = maxL;
  }

  for(int i = 5; i < 2 * 5; i++) {
    double m = sqrt(mU[i] * mU[i] + mV[i] * mV[i]) + sqrt(fabs(gam * mP[i] / mRho[i]));
    double p = sqrt(pU[i] * pU[i] + pV[i] * pV[i]) + sqrt(fabs(gam * pP[i] / pRho[i]));
    double lamda = max(m, p);
    if(i == 5 || lamda > maxL) {
      maxL = lamda;
    }
  }

  for(int i = 5; i < 2 * 5; i++) {
    maxLamda[i] = maxL;
  }

  for(int i = 2 * 5; i < 3 * 5; i++) {
    double m = sqrt(mU[i] * mU[i] + mV[i] * mV[i]) + sqrt(fabs(gam * mP[i] / mRho[i]));
    double p = sqrt(pU[i] * pU[i] + pV[i] * pV[i]) + sqrt(fabs(gam * pP[i] / pRho[i]));
    double lamda = max(m, p);
    if(i == 2 * 5 || lamda > maxL) {
      maxL = lamda;
    }
  }

  for(int i = 2 * 5; i < 3 * 5; i++) {
    maxLamda[i] = maxL;
  }

  for(int j = 0; j < 3 * 5; j++) {
    flux0[j] = nx[j] * (pF0[j] + mF0[j]) + ny[j] * (pG0[j] + mG0[j])
               + maxLamda[j] * (mQ0[j] - pQ0[j]);
    flux0[j] *= fscale[j] * 0.5;
  }

  for(int j = 0; j < 3 * 5; j++) {
    flux1[j] = nx[j] * (pF1[j] + mF1[j]) + ny[j] * (pG1[j] + mG1[j])
               + maxLamda[j] * (mQ1[j] - pQ1[j]);
    flux1[j] *= fscale[j] * 0.5;
  }

  for(int j = 0; j < 3 * 5; j++) {
    flux2[j] = nx[j] * (pF2[j] + mF2[j]) + ny[j] * (pG2[j] + mG2[j])
               + maxLamda[j] * (mQ2[j] - pQ2[j]);
    flux2[j] *= fscale[j] * 0.5;
  }

  for(int j = 0; j < 3 * 5; j++) {
    flux3[j] = nx[j] * (pF3[j] + mF3[j]) + ny[j] * (pG3[j] + mG3[j])
               + maxLamda[j] * (mQ3[j] - pQ3[j]);
    flux3[j] *= fscale[j] * 0.5;
  }
}

__device__ void roe(double *flux0, double *flux1, double *flux2, double *flux3,
                    const double *nx, const double *ny, const double *fscale,
                    const double *q0, const double *q1, const double *q2,
                    const double *q3, const double *pQ0, const double *pQ1,
                    const double *pQ2, const double *pQ3) {
  double gam = 1.4;
  int FMASK[15] = {0, 1, 2, 3, 4, 4, 8, 11, 13, 14, 14, 12, 9, 5, 0};
  // Compute primative variables and flux functions for interior points on edges
  double mQ0[3 * 5]; double mQ1[3 * 5]; double mQ2[3 * 5]; double mQ3[3 * 5];
  double mF0[3 * 5]; double mF1[3 * 5]; double mF2[3 * 5]; double mF3[3 * 5];
  double mG0[3 * 5]; double mG1[3 * 5]; double mG2[3 * 5]; double mG3[3 * 5];
  double mRho[3 * 5];
  double mU[3 * 5];
  double mV[3 * 5];
  double mP[3 * 5];

  for(int i = 0; i < 3 * 5; i++) {
    int ind = FMASK[i];
    mQ0[i] = q0[ind];
    mQ1[i] = q1[ind];
    mQ2[i] = q2[ind];
    mQ3[i] = q3[ind];

    euler_flux(mQ0[i], mQ1[i], mQ2[i], mQ3[i], &mF0[i], &mF1[i], &mF2[i],
               &mF3[i], &mG0[i], &mG1[i], &mG2[i], &mG3[i], &mRho[i], &mU[i],
               &mV[i], &mP[i]);
  }

  // Compute primative variables and flux functions for exterior points on edges
  double pF0[3 * 5]; double pF1[3 * 5]; double pF2[3 * 5]; double pF3[3 * 5];
  double pG0[3 * 5]; double pG1[3 * 5]; double pG2[3 * 5]; double pG3[3 * 5];
  double pRho[3 * 5];
  double pU[3 * 5];
  double pV[3 * 5];
  double pP[3 * 5];
  for(int i = 0; i < 3 * 5; i++) {
    euler_flux(pQ0[i], pQ1[i], pQ2[i], pQ3[i], &pF0[i], &pF1[i], &pF2[i],
               &pF3[i], &pG0[i], &pG1[i], &pG2[i], &pG3[i], &pRho[i], &pU[i],
               &pV[i], &pP[i]);
  }

  for(int i = 0; i < 3 * 5; i++) {
    double mRoeQ0 = mQ0[i];
    double mRoeQ1 = nx[i] * mQ1[i] + ny[i] * mQ2[i];
    double mRoeQ2 = -ny[i] * mQ1[i] + nx[i] * mQ2[i];
    double mRoeQ3 = mQ3[i];

    double pRoeQ0 = pQ0[i];
    double pRoeQ1 = nx[i] * pQ1[i] + ny[i] * pQ2[i];
    double pRoeQ2 = -ny[i] * pQ1[i] + nx[i] * pQ2[i];
    double pRoeQ3 = pQ3[i];

    double fxMQ0, fxMQ1, fxMQ2, fxMQ3;
    double fyMQ0, fyMQ1, fyMQ2, fyMQ3;
    double mRoeRho, mRoeU, mRoeV, mRoeP;
    euler_flux(mRoeQ0, mRoeQ1, mRoeQ2, mRoeQ3, &fxMQ0, &fxMQ1, &fxMQ2, &fxMQ3,
               &fyMQ0, &fyMQ1, &fyMQ2, &fyMQ3, &mRoeRho, &mRoeU, &mRoeV,
               &mRoeP);

    double fxPQ0, fxPQ1, fxPQ2, fxPQ3;
    double fyPQ0, fyPQ1, fyPQ2, fyPQ3;
    double pRoeRho, pRoeU, pRoeV, pRoeP;
    euler_flux(pRoeQ0, pRoeQ1, pRoeQ2, pRoeQ3, &fxPQ0, &fxPQ1, &fxPQ2, &fxPQ3,
               &fyPQ0, &fyPQ1, &fyPQ2, &fyPQ3, &pRoeRho, &pRoeU, &pRoeV,
               &pRoeP);

    double mH = (mRoeQ3 + mRoeP) / mRoeRho;
    double pH = (pRoeQ3 + pRoeP) / pRoeRho;

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

    dW1 = fabs(u - c) * dW1;
    dW2 = fabs(u) * dW2;
    dW3 = fabs(u) * dW3;
    dW4 = fabs(u + c) * dW4;

    double fx[4];
    fx[0] = (fxMQ0 + fxPQ0) / 2.0;
    fx[1] = (fxMQ1 + fxPQ1) / 2.0;
    fx[2] = (fxMQ2 + fxPQ2) / 2.0;
    fx[3] = (fxMQ3 + fxPQ3) / 2.0;

    fx[0] = fx[0] - (dW1 + dW2 + dW4) / 2.0;
    fx[1] = fx[1] - (dW1 * (u - c) + dW2 * u + dW4 * (u + c)) / 2.0;
    fx[2] = fx[2] - (dW1 * v + dW2 * v + dW3 + dW4 * v) / 2.0;
    fx[3] = fx[3] - (dW1 * (H - u * c) + dW2 * (u * u + v * v) + dW3 * v + dW4 * (H + u * c)) / 2.0;

    flux0[i] = fx[0];
    flux1[i] = nx[i] * fx[1] - ny[i] * fx[2];
    flux2[i] = ny[i] * fx[1] + nx[i] * fx[2];
    flux3[i] = fx[3];

    flux0[i] *= fscale[i];
    flux1[i] *= fscale[i];
    flux2[i] *= fscale[i];
    flux3[i] *= fscale[i];
  }
}

#endif
