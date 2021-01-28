//
// auto-generated by op2.py
//

//user function
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
      exteriorQ0[exInd + i] += bc_r;
      exteriorQ1[exInd + i] += bc_r * bc_u;
      exteriorQ2[exInd + i] += bc_r * bc_v;
      exteriorQ3[exInd + i] += q3[qInd];
    }
  } else {
    // Wall
    for(int i = 0; i < 5; i++) {
      int qInd = fmask[i];
      exteriorQ0[exInd + i] += q0[qInd];
      exteriorQ1[exInd + i] += q1[qInd] - 2 * (nx[exInd + i] * q1[qInd] + ny[exInd + i] * q2[qInd]) * nx[exInd + i];
      exteriorQ2[exInd + i] += q2[qInd] - 2 * (nx[exInd + i] * q1[qInd] + ny[exInd + i] * q2[qInd]) * ny[exInd + i];
      exteriorQ3[exInd + i] += q3[qInd];
    }
  }
}
#ifdef VECTORIZE
//user function -- modified for vectorisation
#if defined __clang__ || defined __GNUC__
__attribute__((always_inline))
#endif
inline void get_bedge_q_vec( const int bedge_type[][SIMD_VEC], const int bedgeNum[][SIMD_VEC], const double nx[][SIMD_VEC], const double ny[][SIMD_VEC], const double q0[][SIMD_VEC], const double q1[][SIMD_VEC], const double q2[][SIMD_VEC], const double q3[][SIMD_VEC], double exteriorQ0[][SIMD_VEC], double exteriorQ1[][SIMD_VEC], double exteriorQ2[][SIMD_VEC], double exteriorQ3[][SIMD_VEC], int idx ) {
  int exInd = 0;
  if(bedgeNum[0][idx]== 1) {
    exInd = 5;
  } else if(bedgeNum[0][idx]== 2) {
    exInd = 2 * 5;
  }

  int *fmask;

  if(bedgeNum[0][idx]== 0) {
    fmask = FMASK;
  } else if(bedgeNum[0][idx]== 1) {
    fmask = &FMASK[5];
  } else {
    fmask = &FMASK[2 * 5];
  }

  if(bedge_type[0][idx]== 0) {

    for(int i = 0; i < 5; i++) {
      exteriorQ0[exInd + i][idx] = bc_r;
      exteriorQ1[exInd + i][idx] = bc_r * bc_u;
      exteriorQ2[exInd + i][idx] = bc_r * bc_v;
      exteriorQ3[exInd + i][idx] = bc_e;
    }
  } else if(bedge_type[0][idx]== 1) {

    for(int i = 0; i < 5; i++) {
      int qInd = fmask[i];
      exteriorQ0[exInd + i][idx] = bc_r;
      exteriorQ1[exInd + i][idx] = bc_r * bc_u;
      exteriorQ2[exInd + i][idx] = bc_r * bc_v;
      exteriorQ3[exInd + i][idx] = q3[qInd][idx];
    }
  } else {

    for(int i = 0; i < 5; i++) {
      int qInd = fmask[i];
      exteriorQ0[exInd + i][idx] = q0[qInd][idx];
      exteriorQ1[exInd + i][idx] = q1[qInd][idx] - 2 * (nx[exInd + i][idx] * q1[qInd][idx] + ny[exInd + i][idx] * q2[qInd][idx]) * nx[exInd + i][idx];
      exteriorQ2[exInd + i][idx] = q2[qInd][idx] - 2 * (nx[exInd + i][idx] * q1[qInd][idx] + ny[exInd + i][idx] * q2[qInd][idx]) * ny[exInd + i][idx];
      exteriorQ3[exInd + i][idx] = q3[qInd][idx];
    }
  }

}
#endif

// host stub function
void op_par_loop_get_bedge_q(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11){

  int nargs = 12;
  op_arg args[12];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  args[10] = arg10;
  args[11] = arg11;
  //create aligned pointers for dats
  ALIGNED_int const int * __restrict__ ptr0 = (int *) arg0.data;
  DECLARE_PTR_ALIGNED(ptr0,int_ALIGN);
  ALIGNED_int const int * __restrict__ ptr1 = (int *) arg1.data;
  DECLARE_PTR_ALIGNED(ptr1,int_ALIGN);
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
  ALIGNED_double const double * __restrict__ ptr7 = (double *) arg7.data;
  DECLARE_PTR_ALIGNED(ptr7,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr8 = (double *) arg8.data;
  DECLARE_PTR_ALIGNED(ptr8,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr9 = (double *) arg9.data;
  DECLARE_PTR_ALIGNED(ptr9,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr10 = (double *) arg10.data;
  DECLARE_PTR_ALIGNED(ptr10,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr11 = (double *) arg11.data;
  DECLARE_PTR_ALIGNED(ptr11,double_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(4);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: get_bedge_q\n");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      if ((n+SIMD_VEC >= set->core_size) && (n+SIMD_VEC-set->core_size < SIMD_VEC)) {
        op_mpi_wait_all(nargs, args);
      }
      ALIGNED_int int dat0[1][SIMD_VEC];
      ALIGNED_int int dat1[1][SIMD_VEC];
      ALIGNED_double double dat2[15][SIMD_VEC];
      ALIGNED_double double dat3[15][SIMD_VEC];
      ALIGNED_double double dat4[15][SIMD_VEC];
      ALIGNED_double double dat5[15][SIMD_VEC];
      ALIGNED_double double dat6[15][SIMD_VEC];
      ALIGNED_double double dat7[15][SIMD_VEC];
      ALIGNED_double double dat8[15][SIMD_VEC];
      ALIGNED_double double dat9[15][SIMD_VEC];
      ALIGNED_double double dat10[15][SIMD_VEC];
      ALIGNED_double double dat11[15][SIMD_VEC];
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx0_1 = 1 * (n+i);
        int idx1_1 = 1 * (n+i);
        int idx2_15 = 15 * arg2.map_data[(n+i) * arg2.map->dim + 0];
        int idx3_15 = 15 * arg2.map_data[(n+i) * arg2.map->dim + 0];
        int idx4_15 = 15 * arg2.map_data[(n+i) * arg2.map->dim + 0];
        int idx5_15 = 15 * arg2.map_data[(n+i) * arg2.map->dim + 0];
        int idx6_15 = 15 * arg2.map_data[(n+i) * arg2.map->dim + 0];
        int idx7_15 = 15 * arg2.map_data[(n+i) * arg2.map->dim + 0];

        dat0[0][i] = (ptr0)[idx0_1 + 0];

        dat1[0][i] = (ptr1)[idx1_1 + 0];

        dat2[0][i] = (ptr2)[idx2_15 + 0];
        dat2[1][i] = (ptr2)[idx2_15 + 1];
        dat2[2][i] = (ptr2)[idx2_15 + 2];
        dat2[3][i] = (ptr2)[idx2_15 + 3];
        dat2[4][i] = (ptr2)[idx2_15 + 4];
        dat2[5][i] = (ptr2)[idx2_15 + 5];
        dat2[6][i] = (ptr2)[idx2_15 + 6];
        dat2[7][i] = (ptr2)[idx2_15 + 7];
        dat2[8][i] = (ptr2)[idx2_15 + 8];
        dat2[9][i] = (ptr2)[idx2_15 + 9];
        dat2[10][i] = (ptr2)[idx2_15 + 10];
        dat2[11][i] = (ptr2)[idx2_15 + 11];
        dat2[12][i] = (ptr2)[idx2_15 + 12];
        dat2[13][i] = (ptr2)[idx2_15 + 13];
        dat2[14][i] = (ptr2)[idx2_15 + 14];

        dat3[0][i] = (ptr3)[idx3_15 + 0];
        dat3[1][i] = (ptr3)[idx3_15 + 1];
        dat3[2][i] = (ptr3)[idx3_15 + 2];
        dat3[3][i] = (ptr3)[idx3_15 + 3];
        dat3[4][i] = (ptr3)[idx3_15 + 4];
        dat3[5][i] = (ptr3)[idx3_15 + 5];
        dat3[6][i] = (ptr3)[idx3_15 + 6];
        dat3[7][i] = (ptr3)[idx3_15 + 7];
        dat3[8][i] = (ptr3)[idx3_15 + 8];
        dat3[9][i] = (ptr3)[idx3_15 + 9];
        dat3[10][i] = (ptr3)[idx3_15 + 10];
        dat3[11][i] = (ptr3)[idx3_15 + 11];
        dat3[12][i] = (ptr3)[idx3_15 + 12];
        dat3[13][i] = (ptr3)[idx3_15 + 13];
        dat3[14][i] = (ptr3)[idx3_15 + 14];

        dat4[0][i] = (ptr4)[idx4_15 + 0];
        dat4[1][i] = (ptr4)[idx4_15 + 1];
        dat4[2][i] = (ptr4)[idx4_15 + 2];
        dat4[3][i] = (ptr4)[idx4_15 + 3];
        dat4[4][i] = (ptr4)[idx4_15 + 4];
        dat4[5][i] = (ptr4)[idx4_15 + 5];
        dat4[6][i] = (ptr4)[idx4_15 + 6];
        dat4[7][i] = (ptr4)[idx4_15 + 7];
        dat4[8][i] = (ptr4)[idx4_15 + 8];
        dat4[9][i] = (ptr4)[idx4_15 + 9];
        dat4[10][i] = (ptr4)[idx4_15 + 10];
        dat4[11][i] = (ptr4)[idx4_15 + 11];
        dat4[12][i] = (ptr4)[idx4_15 + 12];
        dat4[13][i] = (ptr4)[idx4_15 + 13];
        dat4[14][i] = (ptr4)[idx4_15 + 14];

        dat5[0][i] = (ptr5)[idx5_15 + 0];
        dat5[1][i] = (ptr5)[idx5_15 + 1];
        dat5[2][i] = (ptr5)[idx5_15 + 2];
        dat5[3][i] = (ptr5)[idx5_15 + 3];
        dat5[4][i] = (ptr5)[idx5_15 + 4];
        dat5[5][i] = (ptr5)[idx5_15 + 5];
        dat5[6][i] = (ptr5)[idx5_15 + 6];
        dat5[7][i] = (ptr5)[idx5_15 + 7];
        dat5[8][i] = (ptr5)[idx5_15 + 8];
        dat5[9][i] = (ptr5)[idx5_15 + 9];
        dat5[10][i] = (ptr5)[idx5_15 + 10];
        dat5[11][i] = (ptr5)[idx5_15 + 11];
        dat5[12][i] = (ptr5)[idx5_15 + 12];
        dat5[13][i] = (ptr5)[idx5_15 + 13];
        dat5[14][i] = (ptr5)[idx5_15 + 14];

        dat6[0][i] = (ptr6)[idx6_15 + 0];
        dat6[1][i] = (ptr6)[idx6_15 + 1];
        dat6[2][i] = (ptr6)[idx6_15 + 2];
        dat6[3][i] = (ptr6)[idx6_15 + 3];
        dat6[4][i] = (ptr6)[idx6_15 + 4];
        dat6[5][i] = (ptr6)[idx6_15 + 5];
        dat6[6][i] = (ptr6)[idx6_15 + 6];
        dat6[7][i] = (ptr6)[idx6_15 + 7];
        dat6[8][i] = (ptr6)[idx6_15 + 8];
        dat6[9][i] = (ptr6)[idx6_15 + 9];
        dat6[10][i] = (ptr6)[idx6_15 + 10];
        dat6[11][i] = (ptr6)[idx6_15 + 11];
        dat6[12][i] = (ptr6)[idx6_15 + 12];
        dat6[13][i] = (ptr6)[idx6_15 + 13];
        dat6[14][i] = (ptr6)[idx6_15 + 14];

        dat7[0][i] = (ptr7)[idx7_15 + 0];
        dat7[1][i] = (ptr7)[idx7_15 + 1];
        dat7[2][i] = (ptr7)[idx7_15 + 2];
        dat7[3][i] = (ptr7)[idx7_15 + 3];
        dat7[4][i] = (ptr7)[idx7_15 + 4];
        dat7[5][i] = (ptr7)[idx7_15 + 5];
        dat7[6][i] = (ptr7)[idx7_15 + 6];
        dat7[7][i] = (ptr7)[idx7_15 + 7];
        dat7[8][i] = (ptr7)[idx7_15 + 8];
        dat7[9][i] = (ptr7)[idx7_15 + 9];
        dat7[10][i] = (ptr7)[idx7_15 + 10];
        dat7[11][i] = (ptr7)[idx7_15 + 11];
        dat7[12][i] = (ptr7)[idx7_15 + 12];
        dat7[13][i] = (ptr7)[idx7_15 + 13];
        dat7[14][i] = (ptr7)[idx7_15 + 14];

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

        dat9[0][i] = 0.0;
        dat9[1][i] = 0.0;
        dat9[2][i] = 0.0;
        dat9[3][i] = 0.0;
        dat9[4][i] = 0.0;
        dat9[5][i] = 0.0;
        dat9[6][i] = 0.0;
        dat9[7][i] = 0.0;
        dat9[8][i] = 0.0;
        dat9[9][i] = 0.0;
        dat9[10][i] = 0.0;
        dat9[11][i] = 0.0;
        dat9[12][i] = 0.0;
        dat9[13][i] = 0.0;
        dat9[14][i] = 0.0;

        dat10[0][i] = 0.0;
        dat10[1][i] = 0.0;
        dat10[2][i] = 0.0;
        dat10[3][i] = 0.0;
        dat10[4][i] = 0.0;
        dat10[5][i] = 0.0;
        dat10[6][i] = 0.0;
        dat10[7][i] = 0.0;
        dat10[8][i] = 0.0;
        dat10[9][i] = 0.0;
        dat10[10][i] = 0.0;
        dat10[11][i] = 0.0;
        dat10[12][i] = 0.0;
        dat10[13][i] = 0.0;
        dat10[14][i] = 0.0;

        dat11[0][i] = 0.0;
        dat11[1][i] = 0.0;
        dat11[2][i] = 0.0;
        dat11[3][i] = 0.0;
        dat11[4][i] = 0.0;
        dat11[5][i] = 0.0;
        dat11[6][i] = 0.0;
        dat11[7][i] = 0.0;
        dat11[8][i] = 0.0;
        dat11[9][i] = 0.0;
        dat11[10][i] = 0.0;
        dat11[11][i] = 0.0;
        dat11[12][i] = 0.0;
        dat11[13][i] = 0.0;
        dat11[14][i] = 0.0;

      }
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        get_bedge_q_vec(
          dat0,
          dat1,
          dat2,
          dat3,
          dat4,
          dat5,
          dat6,
          dat7,
          dat8,
          dat9,
          dat10,
          dat11,
          i);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx8_15 = 15 * arg2.map_data[(n+i) * arg2.map->dim + 0];
        int idx9_15 = 15 * arg2.map_data[(n+i) * arg2.map->dim + 0];
        int idx10_15 = 15 * arg2.map_data[(n+i) * arg2.map->dim + 0];
        int idx11_15 = 15 * arg2.map_data[(n+i) * arg2.map->dim + 0];

        (ptr8)[idx8_15 + 0] += dat8[0][i];
        (ptr8)[idx8_15 + 1] += dat8[1][i];
        (ptr8)[idx8_15 + 2] += dat8[2][i];
        (ptr8)[idx8_15 + 3] += dat8[3][i];
        (ptr8)[idx8_15 + 4] += dat8[4][i];
        (ptr8)[idx8_15 + 5] += dat8[5][i];
        (ptr8)[idx8_15 + 6] += dat8[6][i];
        (ptr8)[idx8_15 + 7] += dat8[7][i];
        (ptr8)[idx8_15 + 8] += dat8[8][i];
        (ptr8)[idx8_15 + 9] += dat8[9][i];
        (ptr8)[idx8_15 + 10] += dat8[10][i];
        (ptr8)[idx8_15 + 11] += dat8[11][i];
        (ptr8)[idx8_15 + 12] += dat8[12][i];
        (ptr8)[idx8_15 + 13] += dat8[13][i];
        (ptr8)[idx8_15 + 14] += dat8[14][i];

        (ptr9)[idx9_15 + 0] += dat9[0][i];
        (ptr9)[idx9_15 + 1] += dat9[1][i];
        (ptr9)[idx9_15 + 2] += dat9[2][i];
        (ptr9)[idx9_15 + 3] += dat9[3][i];
        (ptr9)[idx9_15 + 4] += dat9[4][i];
        (ptr9)[idx9_15 + 5] += dat9[5][i];
        (ptr9)[idx9_15 + 6] += dat9[6][i];
        (ptr9)[idx9_15 + 7] += dat9[7][i];
        (ptr9)[idx9_15 + 8] += dat9[8][i];
        (ptr9)[idx9_15 + 9] += dat9[9][i];
        (ptr9)[idx9_15 + 10] += dat9[10][i];
        (ptr9)[idx9_15 + 11] += dat9[11][i];
        (ptr9)[idx9_15 + 12] += dat9[12][i];
        (ptr9)[idx9_15 + 13] += dat9[13][i];
        (ptr9)[idx9_15 + 14] += dat9[14][i];

        (ptr10)[idx10_15 + 0] += dat10[0][i];
        (ptr10)[idx10_15 + 1] += dat10[1][i];
        (ptr10)[idx10_15 + 2] += dat10[2][i];
        (ptr10)[idx10_15 + 3] += dat10[3][i];
        (ptr10)[idx10_15 + 4] += dat10[4][i];
        (ptr10)[idx10_15 + 5] += dat10[5][i];
        (ptr10)[idx10_15 + 6] += dat10[6][i];
        (ptr10)[idx10_15 + 7] += dat10[7][i];
        (ptr10)[idx10_15 + 8] += dat10[8][i];
        (ptr10)[idx10_15 + 9] += dat10[9][i];
        (ptr10)[idx10_15 + 10] += dat10[10][i];
        (ptr10)[idx10_15 + 11] += dat10[11][i];
        (ptr10)[idx10_15 + 12] += dat10[12][i];
        (ptr10)[idx10_15 + 13] += dat10[13][i];
        (ptr10)[idx10_15 + 14] += dat10[14][i];

        (ptr11)[idx11_15 + 0] += dat11[0][i];
        (ptr11)[idx11_15 + 1] += dat11[1][i];
        (ptr11)[idx11_15 + 2] += dat11[2][i];
        (ptr11)[idx11_15 + 3] += dat11[3][i];
        (ptr11)[idx11_15 + 4] += dat11[4][i];
        (ptr11)[idx11_15 + 5] += dat11[5][i];
        (ptr11)[idx11_15 + 6] += dat11[6][i];
        (ptr11)[idx11_15 + 7] += dat11[7][i];
        (ptr11)[idx11_15 + 8] += dat11[8][i];
        (ptr11)[idx11_15 + 9] += dat11[9][i];
        (ptr11)[idx11_15 + 10] += dat11[10][i];
        (ptr11)[idx11_15 + 11] += dat11[11][i];
        (ptr11)[idx11_15 + 12] += dat11[12][i];
        (ptr11)[idx11_15 + 13] += dat11[13][i];
        (ptr11)[idx11_15 + 14] += dat11[14][i];

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
      int map2idx;
      map2idx = arg2.map_data[n * arg2.map->dim + 0];

      get_bedge_q(
        &(ptr0)[1 * n],
        &(ptr1)[1 * n],
        &(ptr2)[15 * map2idx],
        &(ptr3)[15 * map2idx],
        &(ptr4)[15 * map2idx],
        &(ptr5)[15 * map2idx],
        &(ptr6)[15 * map2idx],
        &(ptr7)[15 * map2idx],
        &(ptr8)[15 * map2idx],
        &(ptr9)[15 * map2idx],
        &(ptr10)[15 * map2idx],
        &(ptr11)[15 * map2idx]);
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
  OP_kernels[4].transfer += (float)set->size * arg2.size;
  OP_kernels[4].transfer += (float)set->size * arg3.size;
  OP_kernels[4].transfer += (float)set->size * arg4.size;
  OP_kernels[4].transfer += (float)set->size * arg5.size;
  OP_kernels[4].transfer += (float)set->size * arg6.size;
  OP_kernels[4].transfer += (float)set->size * arg8.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg9.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg11.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg0.size;
  OP_kernels[4].transfer += (float)set->size * arg1.size;
  OP_kernels[4].transfer += (float)set->size * arg2.map->dim * 4.0f;
}
