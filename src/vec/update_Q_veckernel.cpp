//
// auto-generated by op2.py
//

//user function
inline void update_Q(const double *dt, double *q0, double *q1, double *q2,
                     double *q3, const double *rk10, const double *rk11,
                     const double *rk12, const double *rk13, const double *rk20,
                     const double *rk21, const double *rk22, const double *rk23,
                     const double *rk30, const double *rk31, const double *rk32,
                     const double *rk33, double *workingQ0, double *workingQ1,
                     double *workingQ2, double *workingQ3) {
  for(int i = 0; i < 15; i++) {
    q0[i] = q0[i] + (*dt) * (rk10[i]/ 6.0 + rk20[i] / 6.0 + 2.0 * rk30[i] / 3.0);
    workingQ0[i] = q0[i];
    q1[i] = q1[i] + (*dt) * (rk11[i]/ 6.0 + rk21[i] / 6.0 + 2.0 * rk31[i] / 3.0);
    workingQ1[i] = q1[i];
    q2[i] = q2[i] + (*dt) * (rk12[i]/ 6.0 + rk22[i] / 6.0 + 2.0 * rk32[i] / 3.0);
    workingQ2[i] = q2[i];
    q3[i] = q3[i] + (*dt) * (rk13[i]/ 6.0 + rk23[i] / 6.0 + 2.0 * rk33[i] / 3.0);
    workingQ3[i] = q3[i];
  }
}

// host stub function
void op_par_loop_update_Q(char const *name, op_set set,
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
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  op_arg arg15,
  op_arg arg16,
  op_arg arg17,
  op_arg arg18,
  op_arg arg19,
  op_arg arg20){

  int nargs = 21;
  op_arg args[21];

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
  args[12] = arg12;
  args[13] = arg13;
  args[14] = arg14;
  args[15] = arg15;
  args[16] = arg16;
  args[17] = arg17;
  args[18] = arg18;
  args[19] = arg19;
  args[20] = arg20;
  //create aligned pointers for dats
  ALIGNED_double       double * __restrict__ ptr1 = (double *) arg1.data;
  DECLARE_PTR_ALIGNED(ptr1,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr2 = (double *) arg2.data;
  DECLARE_PTR_ALIGNED(ptr2,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr3 = (double *) arg3.data;
  DECLARE_PTR_ALIGNED(ptr3,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr4 = (double *) arg4.data;
  DECLARE_PTR_ALIGNED(ptr4,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr5 = (double *) arg5.data;
  DECLARE_PTR_ALIGNED(ptr5,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr6 = (double *) arg6.data;
  DECLARE_PTR_ALIGNED(ptr6,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr7 = (double *) arg7.data;
  DECLARE_PTR_ALIGNED(ptr7,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr8 = (double *) arg8.data;
  DECLARE_PTR_ALIGNED(ptr8,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr9 = (double *) arg9.data;
  DECLARE_PTR_ALIGNED(ptr9,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr10 = (double *) arg10.data;
  DECLARE_PTR_ALIGNED(ptr10,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr11 = (double *) arg11.data;
  DECLARE_PTR_ALIGNED(ptr11,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr12 = (double *) arg12.data;
  DECLARE_PTR_ALIGNED(ptr12,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr13 = (double *) arg13.data;
  DECLARE_PTR_ALIGNED(ptr13,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr14 = (double *) arg14.data;
  DECLARE_PTR_ALIGNED(ptr14,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr15 = (double *) arg15.data;
  DECLARE_PTR_ALIGNED(ptr15,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr16 = (double *) arg16.data;
  DECLARE_PTR_ALIGNED(ptr16,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr17 = (double *) arg17.data;
  DECLARE_PTR_ALIGNED(ptr17,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr18 = (double *) arg18.data;
  DECLARE_PTR_ALIGNED(ptr18,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr19 = (double *) arg19.data;
  DECLARE_PTR_ALIGNED(ptr19,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr20 = (double *) arg20.data;
  DECLARE_PTR_ALIGNED(ptr20,double_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(4);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  update_Q");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      double dat0[SIMD_VEC];
      for ( int i=0; i<SIMD_VEC; i++ ){
        dat0[i] = *((double*)arg0.data);
      }
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        update_Q(
          &dat0[i],
          &(ptr1)[15 * (n+i)],
          &(ptr2)[15 * (n+i)],
          &(ptr3)[15 * (n+i)],
          &(ptr4)[15 * (n+i)],
          &(ptr5)[15 * (n+i)],
          &(ptr6)[15 * (n+i)],
          &(ptr7)[15 * (n+i)],
          &(ptr8)[15 * (n+i)],
          &(ptr9)[15 * (n+i)],
          &(ptr10)[15 * (n+i)],
          &(ptr11)[15 * (n+i)],
          &(ptr12)[15 * (n+i)],
          &(ptr13)[15 * (n+i)],
          &(ptr14)[15 * (n+i)],
          &(ptr15)[15 * (n+i)],
          &(ptr16)[15 * (n+i)],
          &(ptr17)[15 * (n+i)],
          &(ptr18)[15 * (n+i)],
          &(ptr19)[15 * (n+i)],
          &(ptr20)[15 * (n+i)]);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
    }
    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      update_Q(
        (double*)arg0.data,
        &(ptr1)[15*n],
        &(ptr2)[15*n],
        &(ptr3)[15*n],
        &(ptr4)[15*n],
        &(ptr5)[15*n],
        &(ptr6)[15*n],
        &(ptr7)[15*n],
        &(ptr8)[15*n],
        &(ptr9)[15*n],
        &(ptr10)[15*n],
        &(ptr11)[15*n],
        &(ptr12)[15*n],
        &(ptr13)[15*n],
        &(ptr14)[15*n],
        &(ptr15)[15*n],
        &(ptr16)[15*n],
        &(ptr17)[15*n],
        &(ptr18)[15*n],
        &(ptr19)[15*n],
        &(ptr20)[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[4].name      = name;
  OP_kernels[4].count    += 1;
  OP_kernels[4].time     += wall_t2 - wall_t1;
  OP_kernels[4].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg5.size;
  OP_kernels[4].transfer += (float)set->size * arg6.size;
  OP_kernels[4].transfer += (float)set->size * arg7.size;
  OP_kernels[4].transfer += (float)set->size * arg8.size;
  OP_kernels[4].transfer += (float)set->size * arg9.size;
  OP_kernels[4].transfer += (float)set->size * arg10.size;
  OP_kernels[4].transfer += (float)set->size * arg11.size;
  OP_kernels[4].transfer += (float)set->size * arg12.size;
  OP_kernels[4].transfer += (float)set->size * arg13.size;
  OP_kernels[4].transfer += (float)set->size * arg14.size;
  OP_kernels[4].transfer += (float)set->size * arg15.size;
  OP_kernels[4].transfer += (float)set->size * arg16.size;
  OP_kernels[4].transfer += (float)set->size * arg17.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg18.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg19.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg20.size * 2.0f;
}
