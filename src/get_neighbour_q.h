// fmask0 = 0, 1, 2, 3, 4
// fmask1 = 4, 8, 11, 13, 14
// fmask2 = 0, 5, 9, 12, 14

inline void get_neighbour_q(const double *n0, const double *n1, const double *xL,
                            const double *yL, const double *xR, const double *yR,
                            const double *qL, const double *qR,
                            double *exteriorQL, double *exteriorQR) {
  // Work out which edge for each element
  int edgeL;
  if((n0[0] == xL[0] && n0[1] == yL[0]) || (n1[0] == xL[0] && n1[1] == yL[0])) {
    // Check whether edge to node 1 or 2 of cell
    if((n0[0] == xL[1] && n0[1] == yL[1]) || (n1[0] == xL[1] && n1[1] == yL[1])) {
      // Edge from node 0 to 1 of cell, i.e. edge 0
      edgeL = 0;
    } else {
      // Edge from node 2 to 0 of cell, i.e. edge 2
      edgeL = 2;
    }
  } else {
    // If node 0 of cell not one of the points, then edge is between node 1 and 2 (edge 1)
    edgeL = 1;
  }

  int edgeR;
  if((n0[0] == xR[0] && n0[1] == yR[0]) || (n1[0] == xR[0] && n1[1] == yR[0])) {
    // Check whether edge to node 1 or 2 of cell
    if((n0[0] == xR[1] && n0[1] == yR[1]) || (n1[0] == xR[1] && n1[1] == yR[1])) {
      // Edge from node 0 to 1 of cell, i.e. edge 0
      edgeR = 0;
    } else {
      // Edge from node 2 to 0 of cell, i.e. edge 2
      edgeR = 2;
    }
  } else {
    // If node 0 of cell not one of the points, then edge is between node 1 and 2 (edge 1)
    edgeR = 1;
  }

  // Copy data from R to L
  int exInd = 0;
  if(edgeL == 1) exInd = 4 * NUM_FACE_PTS;
  else if(edgeL == 2) exInd = 2 * 4 * NUM_FACE_PTS;

  int *fmask;

  if(edgeR == 0) {
    fmask = fmask0;
  } else if(edgeR == 1) {
    fmask = fmask1;
  } else {
    fmask = fmask2;
  }

  // As all edges go around clockwise, have to do reverse order for one side of the copy
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    int rInd = fmask[NUM_FACE_PTS - i - 1] * 4;
    exteriorQL[exInd + i * 4]     = qR[rInd];
    exteriorQL[exInd + i * 4 + 1] = qR[rInd + 1];
    exteriorQL[exInd + i * 4 + 2] = qR[rInd + 2];
    exteriorQL[exInd + i * 4 + 3] = qR[rInd + 3];
  }

  // Copy data from L to R
  exInd = 0;
  if(edgeR == 1) exInd = 4 * NUM_FACE_PTS;
  else if(edgeR == 2) exInd = 2 * 4 * NUM_FACE_PTS;

  if(edgeL == 0) {
    fmask = fmask0;
  } else if(edgeL == 1) {
    fmask = fmask1;
  } else {
    fmask = fmask2;
  }

  // As all edges go around clockwise, have to do reverse order for one side of the copy
  for(int i = 0; i < NUM_FACE_PTS; i++) {
    int lInd = fmask[NUM_FACE_PTS - i - 1] * 4;
    exteriorQR[exInd + i * 4]     = qL[lInd];
    exteriorQR[exInd + i * 4 + 1] = qL[lInd + 1];
    exteriorQR[exInd + i * 4 + 2] = qL[lInd + 2];
    exteriorQR[exInd + i * 4 + 3] = qL[lInd + 3];
  }
}
