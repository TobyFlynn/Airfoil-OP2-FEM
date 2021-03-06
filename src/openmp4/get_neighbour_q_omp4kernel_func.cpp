//
// auto-generated by op2.py
//

void get_neighbour_q_omp4_kernel(
  int *data0,
  int dat0size,
  int *map1,
  int map1size,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  double *data5,
  int dat5size,
  double *data7,
  int dat7size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size]) \
    map(to: FMASK_ompkernel[:15])\
    map(to:col_reord[0:set_size1],map1[0:map1size],data1[0:dat1size],data2[0:dat2size],data5[0:dat5size],data7[0:dat7size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map1idx;
    int map3idx;
    map1idx = map1[n_op + set_size1 * 0];
    map3idx = map1[n_op + set_size1 * 1];

    //variable mapping
    const int *edgeNum = &data0[2*n_op];
    const double *xL = &data1[3 * map1idx];
    const double *yL = &data2[3 * map1idx];
    const double *xR = &data1[3 * map3idx];
    const double *yR = &data2[3 * map3idx];
    const double *qL = &data5[60 * map1idx];
    const double *qR = &data5[60 * map3idx];
    double *exteriorQL = &data7[60 * map1idx];
    double *exteriorQR = &data7[60 * map3idx];

    //inline function
    

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

    int exInd = 0;
    if(edgeL == 1) exInd = 4 * 5;
    else if(edgeL == 2) exInd = 2 * 4 * 5;

    int *fmask;

    if(edgeR == 0) {
      fmask = FMASK_ompkernel;
    } else if(edgeR == 1) {
      fmask = &FMASK_ompkernel[5];
    } else {
      fmask = &FMASK_ompkernel[2 * 5];
    }

    for(int i = 0; i < 5; i++) {
      int rInd;
      if(reverse) {
        rInd = 4 * fmask[5 - i - 1];
      } else {
        rInd = 4 * fmask[i];
      }
      exteriorQL[exInd + 4 * i]     += qR[rInd];
      exteriorQL[exInd + 4 * i + 1] += qR[rInd + 1];
      exteriorQL[exInd + 4 * i + 2] += qR[rInd + 2];
      exteriorQL[exInd + 4 * i + 3] += qR[rInd + 3];
    }

    exInd = 0;
    if(edgeR == 1) exInd = 4 * 5;
    else if(edgeR == 2) exInd = 2 * 4 * 5;

    if(edgeL == 0) {
      fmask = FMASK_ompkernel;
    } else if(edgeL == 1) {
      fmask = &FMASK_ompkernel[5];
    } else {
      fmask = &FMASK_ompkernel[2 * 5];
    }

    for(int i = 0; i < 5; i++) {
      int lInd;
      if(reverse) {
        lInd = 4 * fmask[5 - i - 1];
      } else {
        lInd = 4 * fmask[i];
      }
      exteriorQR[exInd + 4 * i]     += qL[lInd];
      exteriorQR[exInd + 4 * i + 1] += qL[lInd + 1];
      exteriorQR[exInd + 4 * i + 2] += qL[lInd + 2];
      exteriorQR[exInd + 4 * i + 3] += qL[lInd + 3];
    }
    //end inline func
  }

}
