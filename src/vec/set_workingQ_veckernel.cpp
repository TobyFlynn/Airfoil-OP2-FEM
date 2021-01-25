//
// auto-generated by op2.py
//

//user function
inline void set_workingQ(const double *dt, const int *stage, const double *q,
                         const double *k1, const double *k2, double *workingQ) {
  if(*stage == 0) {
    for(int i = 0; i < 4 * 15; i++) {
      workingQ[i] = q[i] + (*dt) * k1[i];
    }
  } else {
    for(int i = 0; i < 4 * 15; i++) {
      workingQ[i] = q[i] + (*dt) * (k1[i] / 4.0 + k2[i] / 4.0);
    }
  }
}

// host stub function
void op_par_loop_set_workingQ(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  //create aligned pointers for dats
  ALIGNED_double const double * __restrict__ ptr2 = (double *) arg2.data;
  DECLARE_PTR_ALIGNED(ptr2,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr3 = (double *) arg3.data;
  DECLARE_PTR_ALIGNED(ptr3,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr4 = (double *) arg4.data;
  DECLARE_PTR_ALIGNED(ptr4,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr5 = (double *) arg5.data;
  DECLARE_PTR_ALIGNED(ptr5,double_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(7);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  set_workingQ");
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
      int dat1[SIMD_VEC];
      for ( int i=0; i<SIMD_VEC; i++ ){
        dat1[i] = *((int*)arg1.data);
      }
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        set_workingQ(
          &dat0[i],
          &dat1[i],
          &(ptr2)[60 * (n+i)],
          &(ptr3)[60 * (n+i)],
          &(ptr4)[60 * (n+i)],
          &(ptr5)[60 * (n+i)]);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
    }
    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      set_workingQ(
        (double*)arg0.data,
        (int*)arg1.data,
        &(ptr2)[60*n],
        &(ptr3)[60*n],
        &(ptr4)[60*n],
        &(ptr5)[60*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[7].name      = name;
  OP_kernels[7].count    += 1;
  OP_kernels[7].time     += wall_t2 - wall_t1;
  OP_kernels[7].transfer += (float)set->size * arg2.size;
  OP_kernels[7].transfer += (float)set->size * arg3.size;
  OP_kernels[7].transfer += (float)set->size * arg4.size;
  OP_kernels[7].transfer += (float)set->size * arg5.size * 2.0f;
}
