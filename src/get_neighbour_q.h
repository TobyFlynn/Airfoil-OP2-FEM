inline void get_neighbour_q(const int *edgeNum, const double *xL,
                            const double *yL, const double *xR,
                            const double *yR, const double *qL0,
                            const double *qL1, const double *qL2,
                            const double *qL3, const double *qR0,
                            const double *qR1, const double *qR2,
                            const double *qR3, double *exteriorQL0,
                            double *exteriorQL1, double *exteriorQL2,
                            double *exteriorQL3, double *exteriorQR0,
                            double *exteriorQR1, double *exteriorQR2,
                            double *exteriorQR3) {
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
  int exInd = 0;
  if(edgeL == 1) exInd = 5;
  else if(edgeL == 2) exInd = 2 * 5;

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
      rInd = fmask[5 - i - 1];
    } else {
      rInd = fmask[i];
    }
    exteriorQL0[exInd + i] += qR0[rInd];
    exteriorQL1[exInd + i] += qR1[rInd];
    exteriorQL2[exInd + i] += qR2[rInd];
    exteriorQL3[exInd + i] += qR3[rInd];
  }

  // Copy data from L to R
  exInd = 0;
  if(edgeR == 1) exInd = 5;
  else if(edgeR == 2) exInd = 2 * 5;

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
      lInd = fmask[5 - i - 1];
    } else {
      lInd = fmask[i];
    }
    exteriorQR0[exInd + i] += qL0[lInd];
    exteriorQR1[exInd + i] += qL1[lInd];
    exteriorQR2[exInd + i] += qL2[lInd];
    exteriorQR3[exInd + i] += qL3[lInd];
  }
}
