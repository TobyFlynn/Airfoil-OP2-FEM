#include <cmath>

#include "fluxes.h"

inline void get_neighbour_q(const int *edgeNum, const double *xL,
                            const double *yL, const double *xR,
                            const double *yR, const double *nxL,
                            const double *nyL, const double *nxR,
                            const double *nyR, const double *fscaleL,
                            const double *fscaleR, const double *qL,
                            const double *qR, double *fluxL, double *fluxR) {
  double exteriorQL[4 * 5];
  double exteriorQR[4 * 5];
  // Work out which edge for each element
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse;

  if(edgeR == 0) {
    if(edgeL == 0) {
      reverse = !(xL[0] == xR[0] && yL[0] == yR[0]);
    } else if(edgeL == 1) {
      reverse = !(xL[1] == xR[0] && yL[1] == yR[0]);
    } else {
      reverse = !(xL[2] == xR[0] && yL[2] == yR[0]);
    }
  } else if(edgeR == 1) {
    if(edgeL == 0) {
      reverse = !(xL[0] == xR[1] && yL[0] == yR[1]);
    } else if(edgeL == 1) {
      reverse = !(xL[1] == xR[1] && yL[1] == yR[1]);
    } else {
      reverse = !(xL[2] == xR[1] && yL[2] == yR[1]);
    }
  } else {
    if(edgeL == 0) {
      reverse = !(xL[0] == xR[2] && yL[0] == yR[2]);
    } else if(edgeL == 1) {
      reverse = !(xL[1] == xR[2] && yL[1] == yR[2]);
    } else {
      reverse = !(xL[2] == xR[2] && yL[2] == yR[2]);
    }
  }

  // Copy data from R to L
  int exIndL = 0;
  int nIndL = 0;
  if(edgeL == 1) {
    exIndL = 4 * 5;
    nIndL = 5;
  } else if(edgeL == 2) {
    exIndL = 2 * 4 * 5;
    nIndL = 2 * 5;
  }

  int *fmask;

  if(edgeR == 0) {
    fmask = FMASK;
  } else if(edgeR == 1) {
    fmask = &FMASK[5];
  } else {
    fmask = &FMASK[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    int rInd;
    if(reverse) {
      rInd = 4 * fmask[5 - i - 1];
    } else {
      rInd = 4 * fmask[i];
    }
    exteriorQL[4 * i]     = qR[rInd];
    exteriorQL[4 * i + 1] = qR[rInd + 1];
    exteriorQL[4 * i + 2] = qR[rInd + 2];
    exteriorQL[4 * i + 3] = qR[rInd + 3];
  }

  // Copy data from L to R
  int exIndR = 0;
  int nIndR = 0;
  if(edgeR == 1) {
    exIndR = 4 * 5;
    nIndR = 5;
  } else if(edgeR == 2) {
    exIndR = 2 * 4 * 5;
    nIndR = 2 * 5;
  }

  if(edgeL == 0) {
    fmask = FMASK;
  } else if(edgeL == 1) {
    fmask = &FMASK[5];
  } else {
    fmask = &FMASK[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    int lInd;
    if(reverse) {
      lInd = 4 * fmask[5 - i - 1];
    } else {
      lInd = 4 * fmask[i];
    }
    exteriorQR[4 * i]     = qL[lInd];
    exteriorQR[4 * i + 1] = qL[lInd + 1];
    exteriorQR[4 * i + 2] = qL[lInd + 2];
    exteriorQR[4 * i + 3] = qL[lInd + 3];
  }

  // Compute numerical fluxes
  // lax_friedrichs(fluxL + exIndL, nxL + nIndL, nyL + nIndL, fscaleL + nIndL, qL, exteriorQL, nIndL);
  roe(fluxL + exIndL, nxL + nIndL, nyL + nIndL, fscaleL + nIndL, qL, exteriorQL, nIndL);

  // lax_friedrichs(fluxR + exIndR, nxR + nIndR, nyR + nIndR, fscaleR + nIndR, qR, exteriorQR, nIndR);
  roe(fluxR + exIndR, nxR + nIndR, nyR + nIndR, fscaleR + nIndR, qR, exteriorQR, nIndR);
}
