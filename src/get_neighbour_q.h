inline void get_neighbour_q(const int *edgeNum, const double *xL,
                            const double *yL, const double *xR, const double *yR,
                            const double *qL, const double *qR,
                            double *exteriorQL, double *exteriorQR) {
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
  if(edgeL == 1) exInd = 4 * NUM_FACE_PTS;
  else if(edgeL == 2) exInd = 2 * 4 * NUM_FACE_PTS;

  int *fmask;

  if(edgeR == 0) {
    fmask = FMASK;
  } else if(edgeR == 1) {
    fmask = &FMASK[NUM_FACE_PTS];
  } else {
    fmask = &FMASK[2 * NUM_FACE_PTS];
  }

  for(int i = 0; i < NUM_FACE_PTS; i++) {
    int rInd;
    if(reverse) {
      rInd = 4 * fmask[NUM_FACE_PTS - i - 1];
    } else {
      rInd = 4 * fmask[i];
    }
    exteriorQL[exInd + 4 * i]     += qR[rInd];
    exteriorQL[exInd + 4 * i + 1] += qR[rInd + 1];
    exteriorQL[exInd + 4 * i + 2] += qR[rInd + 2];
    exteriorQL[exInd + 4 * i + 3] += qR[rInd + 3];
  }

  // Copy data from L to R
  exInd = 0;
  if(edgeR == 1) exInd = 4 * NUM_FACE_PTS;
  else if(edgeR == 2) exInd = 2 * 4 * NUM_FACE_PTS;

  if(edgeL == 0) {
    fmask = FMASK;
  } else if(edgeL == 1) {
    fmask = &FMASK[NUM_FACE_PTS];
  } else {
    fmask = &FMASK[2 * NUM_FACE_PTS];
  }

  for(int i = 0; i < NUM_FACE_PTS; i++) {
    int lInd;
    if(reverse) {
      lInd = 4 * fmask[NUM_FACE_PTS - i - 1];
    } else {
      lInd = 4 * fmask[i];
    }
    exteriorQR[exInd + 4 * i]     += qL[lInd];
    exteriorQR[exInd + 4 * i + 1] += qL[lInd + 1];
    exteriorQR[exInd + 4 * i + 2] += qL[lInd + 2];
    exteriorQR[exInd + 4 * i + 3] += qL[lInd + 3];
  }
}
