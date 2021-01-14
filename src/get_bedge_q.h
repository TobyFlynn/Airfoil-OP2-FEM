inline void get_bedge_q(const int *bedge_type, const double *n0,
                        const double *n1, const double *x, const double *y,
                        const double *nx, const double *ny, const double *q,
                        double *exteriorQ) {
  // Find which edge of cell
  int edge;
  if((n0[0] == x[0] && n0[1] == y[0]) || (n1[0] == x[0] && n1[1] == y[0])) {
    // Check whether edge to node 1 or 2 of cell
    if((n0[0] == x[1] && n0[1] == y[1]) || (n1[0] == x[1] && n1[1] == y[1])) {
      // Edge from node 0 to 1 of cell, i.e. edge 0
      edge = 0;
    } else {
      // Edge from node 2 to 0 of cell, i.e. edge 2
      edge = 2;
    }
  } else {
    // If node 0 of cell not one of the points, then edge is between node 1 and 2 (edge 1)
    edge = 1;
  }

  // bedge_type: 0 = inflow, 1 = outflow, 2 = wall
  if(*bedge_type == 0) {
    // Inflow
    int exInd = 0;
    if(edge == 1) exInd = 4 * NUM_FACE_PTS;
    else if(edge == 2) exInd = 2 * 4 * NUM_FACE_PTS;
    
    for(int i = 0; i < NUM_FACE_PTS; i++) {
      exteriorQ[exInd + i * 4]     = bc_r;
      exteriorQ[exInd + i * 4 + 1] = bc_r * bc_u;
      exteriorQ[exInd + i * 4 + 2] = 0.0;
      exteriorQ[exInd + i * 4 + 3] = bc_e;
    }
  } else if(*bedge_type == 1) {
    // Outflow
    int exInd = 0;
    if(edge == 1) exInd = 4 * NUM_FACE_PTS;
    else if(edge == 2) exInd = 2 * 4 * NUM_FACE_PTS;

    int *fmask;
    int nInd;

    if(edge == 0) {
      fmask = fmask0;
      nInd = 0;
    } else if(edge == 1) {
      fmask = fmask1;
      nInd = NUM_FACE_PTS;
    } else {
      fmask = fmask2;
      nInd = 2 * NUM_FACE_PTS;
    }

    for(int i = 0; i < NUM_FACE_PTS; i++) {
      int qInd = fmask[i] * 4;
      exteriorQ[exInd + i * 4]     = bc_r;
      exteriorQ[exInd + i * 4 + 1] = bc_r * bc_u;
      exteriorQ[exInd + i * 4 + 2] = 0.0;
      exteriorQ[exInd + i * 4 + 3] = q[qInd + 3];
    }
  } else {
    // Wall
    int exInd = 0;
    if(edge == 1) exInd = 4 * NUM_FACE_PTS;
    else if(edge == 2) exInd = 2 * 4 * NUM_FACE_PTS;

    int *fmask;
    int nInd;

    if(edge == 0) {
      fmask = fmask0;
      nInd = 0;
    } else if(edge == 1) {
      fmask = fmask1;
      nInd = NUM_FACE_PTS;
    } else {
      fmask = fmask2;
      nInd = 2 * NUM_FACE_PTS;
    }

    for(int i = 0; i < NUM_FACE_PTS; i++) {
      int qInd = fmask[i] * 4;
      exteriorQ[exInd + i * 4]     = q[qInd];
      exteriorQ[exInd + i * 4 + 1] = q[qInd + 1] - 2 * (nx[nInd + i] * q[qInd + 1] + ny[nInd + i] * q[qInd + 2]) * nx[nInd + i];
      exteriorQ[exInd + i * 4 + 2] = q[qInd + 2] - 2 * (nx[nInd + i] * q[qInd + 1] + ny[nInd + i] * q[qInd + 2]) * ny[nInd + i];
      exteriorQ[exInd + i * 4 + 3] = q[qInd + 3];
    }
  }
}
