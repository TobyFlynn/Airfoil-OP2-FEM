inline void get_bedge_q(const int *bedge_type, const int *bedgeNum,
                        const double *nx, const double *ny,
                        const double *q, double *exteriorQ) {
  int exInd = 0;
  int nInd = 0;
  if(*bedgeNum == 1) {
    exInd = 4 * 5;
    nInd = 5;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 4 * 5;
    nInd = 2 * 5;
  }

  int *fmask;

  if(*bedgeNum == 0) {
    fmask = FMASK;
  } else if(*bedgeNum == 1) {
    fmask = &FMASK[5];
  } else {
    fmask = &FMASK[2 * 5];
  }

  if(*bedge_type == 0) {
    // Inflow
    for(int i = 0; i < 5; i++) {
      exteriorQ[exInd + i * 4]     += bc_r;
      exteriorQ[exInd + i * 4 + 1] += bc_r * bc_u;
      exteriorQ[exInd + i * 4 + 2] +=  bc_r * bc_v;
      exteriorQ[exInd + i * 4 + 3] += bc_e;
    }
  } else if(*bedge_type == 1) {
    // Outflow
    for(int i = 0; i < 5; i++) {
      int qInd = fmask[i] * 4;
      exteriorQ[exInd + i * 4]     += bc_r;
      exteriorQ[exInd + i * 4 + 1] += bc_r * bc_u;
      exteriorQ[exInd + i * 4 + 2] +=  bc_r * bc_v;
      exteriorQ[exInd + i * 4 + 3] += q[qInd + 3];
    }
  } else {
    // Wall
    for(int i = 0; i < 5; i++) {
      int qInd = fmask[i] * 4;
      exteriorQ[exInd + i * 4]     += q[qInd];
      exteriorQ[exInd + i * 4 + 1] += q[qInd + 1] - 2 * (nx[nInd + i] * q[qInd + 1] + ny[nInd + i] * q[qInd + 2]) * nx[nInd + i];
      exteriorQ[exInd + i * 4 + 2] += q[qInd + 2] - 2 * (nx[nInd + i] * q[qInd + 1] + ny[nInd + i] * q[qInd + 2]) * ny[nInd + i];
      exteriorQ[exInd + i * 4 + 3] += q[qInd + 3];
    }
  }
}
