//
// auto-generated by op2.py
//

//user function
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
  if(edgeL == 1) exInd = 4 * 5;
  else if(edgeL == 2) exInd = 2 * 4 * 5;

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
    exteriorQL[exInd + 4 * i]     += qR[rInd];
    exteriorQL[exInd + 4 * i + 1] += qR[rInd + 1];
    exteriorQL[exInd + 4 * i + 2] += qR[rInd + 2];
    exteriorQL[exInd + 4 * i + 3] += qR[rInd + 3];
  }

  // Copy data from L to R
  exInd = 0;
  if(edgeR == 1) exInd = 4 * 5;
  else if(edgeR == 2) exInd = 2 * 4 * 5;

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
    exteriorQR[exInd + 4 * i]     += qL[lInd];
    exteriorQR[exInd + 4 * i + 1] += qL[lInd + 1];
    exteriorQR[exInd + 4 * i + 2] += qL[lInd + 2];
    exteriorQR[exInd + 4 * i + 3] += qL[lInd + 3];
  }
}
#ifdef VECTORIZE
//user function -- modified for vectorisation
#if defined __clang__ || defined __GNUC__
__attribute__((always_inline))
#endif
inline void get_neighbour_q_vec( const int edgeNum[][SIMD_VEC], const double xL[][SIMD_VEC], const double yL[][SIMD_VEC], const double xR[][SIMD_VEC], const double yR[][SIMD_VEC], const double qL[][SIMD_VEC], const double qR[][SIMD_VEC], double exteriorQL[][SIMD_VEC], double exteriorQR[][SIMD_VEC], int idx ) {

  int edgeL = edgeNum[0][idx];
  int edgeR = edgeNum[1][idx];
  bool reverse;

  if(edgeR == 0) {
    if(edgeL == 0) {
      reverse = !(xL[0][idx] == xR[0][idx] && yL[0][idx] == yR[0][idx]);
    } else if(edgeL == 1) {
      reverse = !(xL[1][idx] == xR[0][idx] && yL[1][idx] == yR[0][idx]);
    } else {
      reverse = !(xL[2][idx] == xR[0][idx] && yL[2][idx] == yR[0][idx]);
    }
  } else if(edgeR == 1) {
    if(edgeL == 0) {
      reverse = !(xL[0][idx] == xR[1][idx] && yL[0][idx] == yR[1][idx]);
    } else if(edgeL == 1) {
      reverse = !(xL[1][idx] == xR[1][idx] && yL[1][idx] == yR[1][idx]);
    } else {
      reverse = !(xL[2][idx] == xR[1][idx] && yL[2][idx] == yR[1][idx]);
    }
  } else {
    if(edgeL == 0) {
      reverse = !(xL[0][idx] == xR[2][idx] && yL[0][idx] == yR[2][idx]);
    } else if(edgeL == 1) {
      reverse = !(xL[1][idx] == xR[2][idx] && yL[1][idx] == yR[2][idx]);
    } else {
      reverse = !(xL[2][idx] == xR[2][idx] && yL[2][idx] == yR[2][idx]);
    }
  }

  int exInd = 0;
  if(edgeL == 1) exInd = 4 * 5;
  else if(edgeL == 2) exInd = 2 * 4 * 5;

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
    exteriorQL[exInd + 4 * i][idx]     = qR[rInd][idx];
    exteriorQL[exInd + 4 * i + 1][idx] = qR[rInd + 1][idx];
    exteriorQL[exInd + 4 * i + 2][idx] = qR[rInd + 2][idx];
    exteriorQL[exInd + 4 * i + 3][idx] = qR[rInd + 3][idx];
  }

  exInd = 0;
  if(edgeR == 1) exInd = 4 * 5;
  else if(edgeR == 2) exInd = 2 * 4 * 5;

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
    exteriorQR[exInd + 4 * i][idx]     = qL[lInd][idx];
    exteriorQR[exInd + 4 * i + 1][idx] = qL[lInd + 1][idx];
    exteriorQR[exInd + 4 * i + 2][idx] = qL[lInd + 2][idx];
    exteriorQR[exInd + 4 * i + 3][idx] = qL[lInd + 3][idx];
  }

}
#endif

// host stub function
void op_par_loop_get_neighbour_q(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  int nargs = 9;
  op_arg args[9];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  //create aligned pointers for dats
  ALIGNED_int const int * __restrict__ ptr0 = (int *) arg0.data;
  DECLARE_PTR_ALIGNED(ptr0,int_ALIGN);
  ALIGNED_double const double * __restrict__ ptr1 = (double *) arg1.data;
  DECLARE_PTR_ALIGNED(ptr1,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr2 = (double *) arg2.data;
  DECLARE_PTR_ALIGNED(ptr2,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr3 = (double *) arg3.data;
  DECLARE_PTR_ALIGNED(ptr3,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr4 = (double *) arg4.data;
  DECLARE_PTR_ALIGNED(ptr4,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr5 = (double *) arg5.data;
  DECLARE_PTR_ALIGNED(ptr5,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr6 = (double *) arg6.data;
  DECLARE_PTR_ALIGNED(ptr6,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr7 = (double *) arg7.data;
  DECLARE_PTR_ALIGNED(ptr7,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr8 = (double *) arg8.data;
  DECLARE_PTR_ALIGNED(ptr8,double_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(4);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: get_neighbour_q\n");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      if ((n+SIMD_VEC >= set->core_size) && (n+SIMD_VEC-set->core_size < SIMD_VEC)) {
        op_mpi_wait_all(nargs, args);
      }
      ALIGNED_int int dat0[2][SIMD_VEC];
      ALIGNED_double double dat1[3][SIMD_VEC];
      ALIGNED_double double dat2[3][SIMD_VEC];
      ALIGNED_double double dat3[3][SIMD_VEC];
      ALIGNED_double double dat4[3][SIMD_VEC];
      ALIGNED_double double dat5[60][SIMD_VEC];
      ALIGNED_double double dat6[60][SIMD_VEC];
      ALIGNED_double double dat7[60][SIMD_VEC];
      ALIGNED_double double dat8[60][SIMD_VEC];
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx0_2 = 2 * (n+i);
        int idx1_3 = 3 * arg1.map_data[(n+i) * arg1.map->dim + 0];
        int idx2_3 = 3 * arg1.map_data[(n+i) * arg1.map->dim + 0];
        int idx3_3 = 3 * arg1.map_data[(n+i) * arg1.map->dim + 1];
        int idx4_3 = 3 * arg1.map_data[(n+i) * arg1.map->dim + 1];
        int idx5_60 = 60 * arg1.map_data[(n+i) * arg1.map->dim + 0];
        int idx6_60 = 60 * arg1.map_data[(n+i) * arg1.map->dim + 1];

        dat0[0][i] = (ptr0)[idx0_2 + 0];
        dat0[1][i] = (ptr0)[idx0_2 + 1];

        dat1[0][i] = (ptr1)[idx1_3 + 0];
        dat1[1][i] = (ptr1)[idx1_3 + 1];
        dat1[2][i] = (ptr1)[idx1_3 + 2];

        dat2[0][i] = (ptr2)[idx2_3 + 0];
        dat2[1][i] = (ptr2)[idx2_3 + 1];
        dat2[2][i] = (ptr2)[idx2_3 + 2];

        dat3[0][i] = (ptr3)[idx3_3 + 0];
        dat3[1][i] = (ptr3)[idx3_3 + 1];
        dat3[2][i] = (ptr3)[idx3_3 + 2];

        dat4[0][i] = (ptr4)[idx4_3 + 0];
        dat4[1][i] = (ptr4)[idx4_3 + 1];
        dat4[2][i] = (ptr4)[idx4_3 + 2];

        dat5[0][i] = (ptr5)[idx5_60 + 0];
        dat5[1][i] = (ptr5)[idx5_60 + 1];
        dat5[2][i] = (ptr5)[idx5_60 + 2];
        dat5[3][i] = (ptr5)[idx5_60 + 3];
        dat5[4][i] = (ptr5)[idx5_60 + 4];
        dat5[5][i] = (ptr5)[idx5_60 + 5];
        dat5[6][i] = (ptr5)[idx5_60 + 6];
        dat5[7][i] = (ptr5)[idx5_60 + 7];
        dat5[8][i] = (ptr5)[idx5_60 + 8];
        dat5[9][i] = (ptr5)[idx5_60 + 9];
        dat5[10][i] = (ptr5)[idx5_60 + 10];
        dat5[11][i] = (ptr5)[idx5_60 + 11];
        dat5[12][i] = (ptr5)[idx5_60 + 12];
        dat5[13][i] = (ptr5)[idx5_60 + 13];
        dat5[14][i] = (ptr5)[idx5_60 + 14];
        dat5[15][i] = (ptr5)[idx5_60 + 15];
        dat5[16][i] = (ptr5)[idx5_60 + 16];
        dat5[17][i] = (ptr5)[idx5_60 + 17];
        dat5[18][i] = (ptr5)[idx5_60 + 18];
        dat5[19][i] = (ptr5)[idx5_60 + 19];
        dat5[20][i] = (ptr5)[idx5_60 + 20];
        dat5[21][i] = (ptr5)[idx5_60 + 21];
        dat5[22][i] = (ptr5)[idx5_60 + 22];
        dat5[23][i] = (ptr5)[idx5_60 + 23];
        dat5[24][i] = (ptr5)[idx5_60 + 24];
        dat5[25][i] = (ptr5)[idx5_60 + 25];
        dat5[26][i] = (ptr5)[idx5_60 + 26];
        dat5[27][i] = (ptr5)[idx5_60 + 27];
        dat5[28][i] = (ptr5)[idx5_60 + 28];
        dat5[29][i] = (ptr5)[idx5_60 + 29];
        dat5[30][i] = (ptr5)[idx5_60 + 30];
        dat5[31][i] = (ptr5)[idx5_60 + 31];
        dat5[32][i] = (ptr5)[idx5_60 + 32];
        dat5[33][i] = (ptr5)[idx5_60 + 33];
        dat5[34][i] = (ptr5)[idx5_60 + 34];
        dat5[35][i] = (ptr5)[idx5_60 + 35];
        dat5[36][i] = (ptr5)[idx5_60 + 36];
        dat5[37][i] = (ptr5)[idx5_60 + 37];
        dat5[38][i] = (ptr5)[idx5_60 + 38];
        dat5[39][i] = (ptr5)[idx5_60 + 39];
        dat5[40][i] = (ptr5)[idx5_60 + 40];
        dat5[41][i] = (ptr5)[idx5_60 + 41];
        dat5[42][i] = (ptr5)[idx5_60 + 42];
        dat5[43][i] = (ptr5)[idx5_60 + 43];
        dat5[44][i] = (ptr5)[idx5_60 + 44];
        dat5[45][i] = (ptr5)[idx5_60 + 45];
        dat5[46][i] = (ptr5)[idx5_60 + 46];
        dat5[47][i] = (ptr5)[idx5_60 + 47];
        dat5[48][i] = (ptr5)[idx5_60 + 48];
        dat5[49][i] = (ptr5)[idx5_60 + 49];
        dat5[50][i] = (ptr5)[idx5_60 + 50];
        dat5[51][i] = (ptr5)[idx5_60 + 51];
        dat5[52][i] = (ptr5)[idx5_60 + 52];
        dat5[53][i] = (ptr5)[idx5_60 + 53];
        dat5[54][i] = (ptr5)[idx5_60 + 54];
        dat5[55][i] = (ptr5)[idx5_60 + 55];
        dat5[56][i] = (ptr5)[idx5_60 + 56];
        dat5[57][i] = (ptr5)[idx5_60 + 57];
        dat5[58][i] = (ptr5)[idx5_60 + 58];
        dat5[59][i] = (ptr5)[idx5_60 + 59];

        dat6[0][i] = (ptr6)[idx6_60 + 0];
        dat6[1][i] = (ptr6)[idx6_60 + 1];
        dat6[2][i] = (ptr6)[idx6_60 + 2];
        dat6[3][i] = (ptr6)[idx6_60 + 3];
        dat6[4][i] = (ptr6)[idx6_60 + 4];
        dat6[5][i] = (ptr6)[idx6_60 + 5];
        dat6[6][i] = (ptr6)[idx6_60 + 6];
        dat6[7][i] = (ptr6)[idx6_60 + 7];
        dat6[8][i] = (ptr6)[idx6_60 + 8];
        dat6[9][i] = (ptr6)[idx6_60 + 9];
        dat6[10][i] = (ptr6)[idx6_60 + 10];
        dat6[11][i] = (ptr6)[idx6_60 + 11];
        dat6[12][i] = (ptr6)[idx6_60 + 12];
        dat6[13][i] = (ptr6)[idx6_60 + 13];
        dat6[14][i] = (ptr6)[idx6_60 + 14];
        dat6[15][i] = (ptr6)[idx6_60 + 15];
        dat6[16][i] = (ptr6)[idx6_60 + 16];
        dat6[17][i] = (ptr6)[idx6_60 + 17];
        dat6[18][i] = (ptr6)[idx6_60 + 18];
        dat6[19][i] = (ptr6)[idx6_60 + 19];
        dat6[20][i] = (ptr6)[idx6_60 + 20];
        dat6[21][i] = (ptr6)[idx6_60 + 21];
        dat6[22][i] = (ptr6)[idx6_60 + 22];
        dat6[23][i] = (ptr6)[idx6_60 + 23];
        dat6[24][i] = (ptr6)[idx6_60 + 24];
        dat6[25][i] = (ptr6)[idx6_60 + 25];
        dat6[26][i] = (ptr6)[idx6_60 + 26];
        dat6[27][i] = (ptr6)[idx6_60 + 27];
        dat6[28][i] = (ptr6)[idx6_60 + 28];
        dat6[29][i] = (ptr6)[idx6_60 + 29];
        dat6[30][i] = (ptr6)[idx6_60 + 30];
        dat6[31][i] = (ptr6)[idx6_60 + 31];
        dat6[32][i] = (ptr6)[idx6_60 + 32];
        dat6[33][i] = (ptr6)[idx6_60 + 33];
        dat6[34][i] = (ptr6)[idx6_60 + 34];
        dat6[35][i] = (ptr6)[idx6_60 + 35];
        dat6[36][i] = (ptr6)[idx6_60 + 36];
        dat6[37][i] = (ptr6)[idx6_60 + 37];
        dat6[38][i] = (ptr6)[idx6_60 + 38];
        dat6[39][i] = (ptr6)[idx6_60 + 39];
        dat6[40][i] = (ptr6)[idx6_60 + 40];
        dat6[41][i] = (ptr6)[idx6_60 + 41];
        dat6[42][i] = (ptr6)[idx6_60 + 42];
        dat6[43][i] = (ptr6)[idx6_60 + 43];
        dat6[44][i] = (ptr6)[idx6_60 + 44];
        dat6[45][i] = (ptr6)[idx6_60 + 45];
        dat6[46][i] = (ptr6)[idx6_60 + 46];
        dat6[47][i] = (ptr6)[idx6_60 + 47];
        dat6[48][i] = (ptr6)[idx6_60 + 48];
        dat6[49][i] = (ptr6)[idx6_60 + 49];
        dat6[50][i] = (ptr6)[idx6_60 + 50];
        dat6[51][i] = (ptr6)[idx6_60 + 51];
        dat6[52][i] = (ptr6)[idx6_60 + 52];
        dat6[53][i] = (ptr6)[idx6_60 + 53];
        dat6[54][i] = (ptr6)[idx6_60 + 54];
        dat6[55][i] = (ptr6)[idx6_60 + 55];
        dat6[56][i] = (ptr6)[idx6_60 + 56];
        dat6[57][i] = (ptr6)[idx6_60 + 57];
        dat6[58][i] = (ptr6)[idx6_60 + 58];
        dat6[59][i] = (ptr6)[idx6_60 + 59];

        dat7[0][i] = 0.0;
        dat7[1][i] = 0.0;
        dat7[2][i] = 0.0;
        dat7[3][i] = 0.0;
        dat7[4][i] = 0.0;
        dat7[5][i] = 0.0;
        dat7[6][i] = 0.0;
        dat7[7][i] = 0.0;
        dat7[8][i] = 0.0;
        dat7[9][i] = 0.0;
        dat7[10][i] = 0.0;
        dat7[11][i] = 0.0;
        dat7[12][i] = 0.0;
        dat7[13][i] = 0.0;
        dat7[14][i] = 0.0;
        dat7[15][i] = 0.0;
        dat7[16][i] = 0.0;
        dat7[17][i] = 0.0;
        dat7[18][i] = 0.0;
        dat7[19][i] = 0.0;
        dat7[20][i] = 0.0;
        dat7[21][i] = 0.0;
        dat7[22][i] = 0.0;
        dat7[23][i] = 0.0;
        dat7[24][i] = 0.0;
        dat7[25][i] = 0.0;
        dat7[26][i] = 0.0;
        dat7[27][i] = 0.0;
        dat7[28][i] = 0.0;
        dat7[29][i] = 0.0;
        dat7[30][i] = 0.0;
        dat7[31][i] = 0.0;
        dat7[32][i] = 0.0;
        dat7[33][i] = 0.0;
        dat7[34][i] = 0.0;
        dat7[35][i] = 0.0;
        dat7[36][i] = 0.0;
        dat7[37][i] = 0.0;
        dat7[38][i] = 0.0;
        dat7[39][i] = 0.0;
        dat7[40][i] = 0.0;
        dat7[41][i] = 0.0;
        dat7[42][i] = 0.0;
        dat7[43][i] = 0.0;
        dat7[44][i] = 0.0;
        dat7[45][i] = 0.0;
        dat7[46][i] = 0.0;
        dat7[47][i] = 0.0;
        dat7[48][i] = 0.0;
        dat7[49][i] = 0.0;
        dat7[50][i] = 0.0;
        dat7[51][i] = 0.0;
        dat7[52][i] = 0.0;
        dat7[53][i] = 0.0;
        dat7[54][i] = 0.0;
        dat7[55][i] = 0.0;
        dat7[56][i] = 0.0;
        dat7[57][i] = 0.0;
        dat7[58][i] = 0.0;
        dat7[59][i] = 0.0;

        dat8[0][i] = 0.0;
        dat8[1][i] = 0.0;
        dat8[2][i] = 0.0;
        dat8[3][i] = 0.0;
        dat8[4][i] = 0.0;
        dat8[5][i] = 0.0;
        dat8[6][i] = 0.0;
        dat8[7][i] = 0.0;
        dat8[8][i] = 0.0;
        dat8[9][i] = 0.0;
        dat8[10][i] = 0.0;
        dat8[11][i] = 0.0;
        dat8[12][i] = 0.0;
        dat8[13][i] = 0.0;
        dat8[14][i] = 0.0;
        dat8[15][i] = 0.0;
        dat8[16][i] = 0.0;
        dat8[17][i] = 0.0;
        dat8[18][i] = 0.0;
        dat8[19][i] = 0.0;
        dat8[20][i] = 0.0;
        dat8[21][i] = 0.0;
        dat8[22][i] = 0.0;
        dat8[23][i] = 0.0;
        dat8[24][i] = 0.0;
        dat8[25][i] = 0.0;
        dat8[26][i] = 0.0;
        dat8[27][i] = 0.0;
        dat8[28][i] = 0.0;
        dat8[29][i] = 0.0;
        dat8[30][i] = 0.0;
        dat8[31][i] = 0.0;
        dat8[32][i] = 0.0;
        dat8[33][i] = 0.0;
        dat8[34][i] = 0.0;
        dat8[35][i] = 0.0;
        dat8[36][i] = 0.0;
        dat8[37][i] = 0.0;
        dat8[38][i] = 0.0;
        dat8[39][i] = 0.0;
        dat8[40][i] = 0.0;
        dat8[41][i] = 0.0;
        dat8[42][i] = 0.0;
        dat8[43][i] = 0.0;
        dat8[44][i] = 0.0;
        dat8[45][i] = 0.0;
        dat8[46][i] = 0.0;
        dat8[47][i] = 0.0;
        dat8[48][i] = 0.0;
        dat8[49][i] = 0.0;
        dat8[50][i] = 0.0;
        dat8[51][i] = 0.0;
        dat8[52][i] = 0.0;
        dat8[53][i] = 0.0;
        dat8[54][i] = 0.0;
        dat8[55][i] = 0.0;
        dat8[56][i] = 0.0;
        dat8[57][i] = 0.0;
        dat8[58][i] = 0.0;
        dat8[59][i] = 0.0;

      }
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        get_neighbour_q_vec(
          dat0,
          dat1,
          dat2,
          dat3,
          dat4,
          dat5,
          dat6,
          dat7,
          dat8,
          i);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx7_60 = 60 * arg1.map_data[(n+i) * arg1.map->dim + 0];
        int idx8_60 = 60 * arg1.map_data[(n+i) * arg1.map->dim + 1];

        (ptr7)[idx7_60 + 0] += dat7[0][i];
        (ptr7)[idx7_60 + 1] += dat7[1][i];
        (ptr7)[idx7_60 + 2] += dat7[2][i];
        (ptr7)[idx7_60 + 3] += dat7[3][i];
        (ptr7)[idx7_60 + 4] += dat7[4][i];
        (ptr7)[idx7_60 + 5] += dat7[5][i];
        (ptr7)[idx7_60 + 6] += dat7[6][i];
        (ptr7)[idx7_60 + 7] += dat7[7][i];
        (ptr7)[idx7_60 + 8] += dat7[8][i];
        (ptr7)[idx7_60 + 9] += dat7[9][i];
        (ptr7)[idx7_60 + 10] += dat7[10][i];
        (ptr7)[idx7_60 + 11] += dat7[11][i];
        (ptr7)[idx7_60 + 12] += dat7[12][i];
        (ptr7)[idx7_60 + 13] += dat7[13][i];
        (ptr7)[idx7_60 + 14] += dat7[14][i];
        (ptr7)[idx7_60 + 15] += dat7[15][i];
        (ptr7)[idx7_60 + 16] += dat7[16][i];
        (ptr7)[idx7_60 + 17] += dat7[17][i];
        (ptr7)[idx7_60 + 18] += dat7[18][i];
        (ptr7)[idx7_60 + 19] += dat7[19][i];
        (ptr7)[idx7_60 + 20] += dat7[20][i];
        (ptr7)[idx7_60 + 21] += dat7[21][i];
        (ptr7)[idx7_60 + 22] += dat7[22][i];
        (ptr7)[idx7_60 + 23] += dat7[23][i];
        (ptr7)[idx7_60 + 24] += dat7[24][i];
        (ptr7)[idx7_60 + 25] += dat7[25][i];
        (ptr7)[idx7_60 + 26] += dat7[26][i];
        (ptr7)[idx7_60 + 27] += dat7[27][i];
        (ptr7)[idx7_60 + 28] += dat7[28][i];
        (ptr7)[idx7_60 + 29] += dat7[29][i];
        (ptr7)[idx7_60 + 30] += dat7[30][i];
        (ptr7)[idx7_60 + 31] += dat7[31][i];
        (ptr7)[idx7_60 + 32] += dat7[32][i];
        (ptr7)[idx7_60 + 33] += dat7[33][i];
        (ptr7)[idx7_60 + 34] += dat7[34][i];
        (ptr7)[idx7_60 + 35] += dat7[35][i];
        (ptr7)[idx7_60 + 36] += dat7[36][i];
        (ptr7)[idx7_60 + 37] += dat7[37][i];
        (ptr7)[idx7_60 + 38] += dat7[38][i];
        (ptr7)[idx7_60 + 39] += dat7[39][i];
        (ptr7)[idx7_60 + 40] += dat7[40][i];
        (ptr7)[idx7_60 + 41] += dat7[41][i];
        (ptr7)[idx7_60 + 42] += dat7[42][i];
        (ptr7)[idx7_60 + 43] += dat7[43][i];
        (ptr7)[idx7_60 + 44] += dat7[44][i];
        (ptr7)[idx7_60 + 45] += dat7[45][i];
        (ptr7)[idx7_60 + 46] += dat7[46][i];
        (ptr7)[idx7_60 + 47] += dat7[47][i];
        (ptr7)[idx7_60 + 48] += dat7[48][i];
        (ptr7)[idx7_60 + 49] += dat7[49][i];
        (ptr7)[idx7_60 + 50] += dat7[50][i];
        (ptr7)[idx7_60 + 51] += dat7[51][i];
        (ptr7)[idx7_60 + 52] += dat7[52][i];
        (ptr7)[idx7_60 + 53] += dat7[53][i];
        (ptr7)[idx7_60 + 54] += dat7[54][i];
        (ptr7)[idx7_60 + 55] += dat7[55][i];
        (ptr7)[idx7_60 + 56] += dat7[56][i];
        (ptr7)[idx7_60 + 57] += dat7[57][i];
        (ptr7)[idx7_60 + 58] += dat7[58][i];
        (ptr7)[idx7_60 + 59] += dat7[59][i];

        (ptr8)[idx8_60 + 0] += dat8[0][i];
        (ptr8)[idx8_60 + 1] += dat8[1][i];
        (ptr8)[idx8_60 + 2] += dat8[2][i];
        (ptr8)[idx8_60 + 3] += dat8[3][i];
        (ptr8)[idx8_60 + 4] += dat8[4][i];
        (ptr8)[idx8_60 + 5] += dat8[5][i];
        (ptr8)[idx8_60 + 6] += dat8[6][i];
        (ptr8)[idx8_60 + 7] += dat8[7][i];
        (ptr8)[idx8_60 + 8] += dat8[8][i];
        (ptr8)[idx8_60 + 9] += dat8[9][i];
        (ptr8)[idx8_60 + 10] += dat8[10][i];
        (ptr8)[idx8_60 + 11] += dat8[11][i];
        (ptr8)[idx8_60 + 12] += dat8[12][i];
        (ptr8)[idx8_60 + 13] += dat8[13][i];
        (ptr8)[idx8_60 + 14] += dat8[14][i];
        (ptr8)[idx8_60 + 15] += dat8[15][i];
        (ptr8)[idx8_60 + 16] += dat8[16][i];
        (ptr8)[idx8_60 + 17] += dat8[17][i];
        (ptr8)[idx8_60 + 18] += dat8[18][i];
        (ptr8)[idx8_60 + 19] += dat8[19][i];
        (ptr8)[idx8_60 + 20] += dat8[20][i];
        (ptr8)[idx8_60 + 21] += dat8[21][i];
        (ptr8)[idx8_60 + 22] += dat8[22][i];
        (ptr8)[idx8_60 + 23] += dat8[23][i];
        (ptr8)[idx8_60 + 24] += dat8[24][i];
        (ptr8)[idx8_60 + 25] += dat8[25][i];
        (ptr8)[idx8_60 + 26] += dat8[26][i];
        (ptr8)[idx8_60 + 27] += dat8[27][i];
        (ptr8)[idx8_60 + 28] += dat8[28][i];
        (ptr8)[idx8_60 + 29] += dat8[29][i];
        (ptr8)[idx8_60 + 30] += dat8[30][i];
        (ptr8)[idx8_60 + 31] += dat8[31][i];
        (ptr8)[idx8_60 + 32] += dat8[32][i];
        (ptr8)[idx8_60 + 33] += dat8[33][i];
        (ptr8)[idx8_60 + 34] += dat8[34][i];
        (ptr8)[idx8_60 + 35] += dat8[35][i];
        (ptr8)[idx8_60 + 36] += dat8[36][i];
        (ptr8)[idx8_60 + 37] += dat8[37][i];
        (ptr8)[idx8_60 + 38] += dat8[38][i];
        (ptr8)[idx8_60 + 39] += dat8[39][i];
        (ptr8)[idx8_60 + 40] += dat8[40][i];
        (ptr8)[idx8_60 + 41] += dat8[41][i];
        (ptr8)[idx8_60 + 42] += dat8[42][i];
        (ptr8)[idx8_60 + 43] += dat8[43][i];
        (ptr8)[idx8_60 + 44] += dat8[44][i];
        (ptr8)[idx8_60 + 45] += dat8[45][i];
        (ptr8)[idx8_60 + 46] += dat8[46][i];
        (ptr8)[idx8_60 + 47] += dat8[47][i];
        (ptr8)[idx8_60 + 48] += dat8[48][i];
        (ptr8)[idx8_60 + 49] += dat8[49][i];
        (ptr8)[idx8_60 + 50] += dat8[50][i];
        (ptr8)[idx8_60 + 51] += dat8[51][i];
        (ptr8)[idx8_60 + 52] += dat8[52][i];
        (ptr8)[idx8_60 + 53] += dat8[53][i];
        (ptr8)[idx8_60 + 54] += dat8[54][i];
        (ptr8)[idx8_60 + 55] += dat8[55][i];
        (ptr8)[idx8_60 + 56] += dat8[56][i];
        (ptr8)[idx8_60 + 57] += dat8[57][i];
        (ptr8)[idx8_60 + 58] += dat8[58][i];
        (ptr8)[idx8_60 + 59] += dat8[59][i];

      }
    }

    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      if (n==set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      int map1idx;
      int map3idx;
      map1idx = arg1.map_data[n * arg1.map->dim + 0];
      map3idx = arg1.map_data[n * arg1.map->dim + 1];

      get_neighbour_q(
        &(ptr0)[2 * n],
        &(ptr1)[3 * map1idx],
        &(ptr2)[3 * map1idx],
        &(ptr3)[3 * map3idx],
        &(ptr4)[3 * map3idx],
        &(ptr5)[60 * map1idx],
        &(ptr6)[60 * map3idx],
        &(ptr7)[60 * map1idx],
        &(ptr8)[60 * map3idx]);
    }
  }

  if (exec_size == 0 || exec_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[4].name      = name;
  OP_kernels[4].count    += 1;
  OP_kernels[4].time     += wall_t2 - wall_t1;
  OP_kernels[4].transfer += (float)set->size * arg1.size;
  OP_kernels[4].transfer += (float)set->size * arg2.size;
  OP_kernels[4].transfer += (float)set->size * arg5.size;
  OP_kernels[4].transfer += (float)set->size * arg0.size;
  OP_kernels[4].transfer += (float)set->size * arg1.map->dim * 4.0f;
}