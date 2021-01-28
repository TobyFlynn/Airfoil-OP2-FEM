inline void get_bedge_q(const int *bedge_type, const int *bedgeNum,
                        const double *nx, const double *ny, const double *q0,
                        const double *q1, const double *q2, const double *q3,
                        double *exteriorQ0, double *exteriorQ1,
                        double *exteriorQ2, double *exteriorQ3) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 5;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 5;
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
      exteriorQ0[exInd + i] += bc_r;
      exteriorQ1[exInd + i] += bc_r * bc_u;
      exteriorQ2[exInd + i] += bc_r * bc_v;
      exteriorQ3[exInd + i] += bc_e;
    }
  } else if(*bedge_type == 1) {
    // Outflow
    for(int i = 0; i < 5; i++) {
      int qInd = fmask[i];
      exteriorQ0[exInd] += bc_r;
      exteriorQ1[exInd] += bc_r * bc_u;
      exteriorQ2[exInd] += bc_r * bc_v;
      exteriorQ3[exInd] += q3[qInd];
    }
  } else {
    // Wall
    for(int i = 0; i < 5; i++) {
      int qInd = fmask[i];
      exteriorQ0[exInd] += q0[qInd];
      exteriorQ1[exInd] += q1[qInd] - 2 * (nx[exInd + i] * q1[qInd] + ny[exInd + i] * q2[qInd]) * nx[exInd + i];
      exteriorQ2[exInd] += q2[qInd] - 2 * (nx[exInd + i] * q1[qInd] + ny[exInd + i] * q2[qInd]) * ny[exInd + i];
      exteriorQ3[exInd] += q3[qInd];
    }
  }
}
