#include <cmath>

inline void get_bedge_q(const int *bedgeNum, const double *n0, const double *n1,
                        const double *nodeX, const double *nodeY, const double *x,
                        const double *y, const double *t, double *exteriorQ) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = 4 * NUM_FACE_PTS;
  else if(*bedgeNum == 2) exInd = 2 * 4 * NUM_FACE_PTS;

  int *fmask;

  if(*bedgeNum == 0) {
    fmask = FMASK;
  } else if(*bedgeNum == 1) {
    fmask = &FMASK[NUM_FACE_PTS];
  } else {
    fmask = &FMASK[2 * NUM_FACE_PTS];
  }

  for(int i = 0; i < NUM_FACE_PTS; i++) {
    double r = (x[fmask[i]] - *t - vortex_x0) * (x[fmask[i]] - *t - vortex_x0) + (y[fmask[i]] - vortex_y0) * (y[fmask[i]] - vortex_y0);
    r = sqrt(r);
    double u = 1.0 - vortex_beta * exp(1.0 - r * r) * ((y[fmask[i]] - vortex_y0) / (2 * M_PI));
    double v = vortex_beta * exp(1.0 - r * r) * ((x[fmask[i]] - vortex_x0) / (2 * M_PI));
    double rho = 1.0 - ((gam - 1) / (16 * gam * M_PI * M_PI)) * vortex_beta * vortex_beta * exp(2 * (1.0 - r * r));
    rho = pow(rho, 1.0 / (gam - 1));
    double p = pow(rho, gam);
    // u = 0.0;
    // v = 0.0;
    // rho = 1.0;
    // p = pow(rho, gam);
    // double ener = (p / (1.0 - gam)) + 0.5 * (rho * u + rho * v);
    double ener = (p / (gam - 1)) + 0.5 * rho * (u * u + v * v);

    exteriorQ[exInd + i * 4]     += rho;
    exteriorQ[exInd + i * 4 + 1] += rho * u;
    exteriorQ[exInd + i * 4 + 2] += rho * v;
    exteriorQ[exInd + i * 4 + 3] += ener;
  }
}
