//
// auto-generated by op2.py
//

//user function
inline void backwards_euler_update_Q(const double *dt, const double *q0,
                     const double *q1, const double *q2, const double *q3,
                     double *rhs0, double *rhs1, double *rhs2, double *rhs3) {
  for(int i = 0; i < 15; i++) {
    rhs0[i] = q0[i] - (*dt) * rhs0[i];
    rhs1[i] = q1[i] - (*dt) * rhs1[i];
    rhs2[i] = q2[i] - (*dt) * rhs2[i];
    rhs3[i] = q3[i] - (*dt) * rhs3[i];
  }
}

// host stub function
void op_par_loop_backwards_euler_update_Q(char const *name, op_set set,
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
  ALIGNED_double const double * __restrict__ ptr1 = (double *) arg1.data;
  DECLARE_PTR_ALIGNED(ptr1,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr2 = (double *) arg2.data;
  DECLARE_PTR_ALIGNED(ptr2,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr3 = (double *) arg3.data;
  DECLARE_PTR_ALIGNED(ptr3,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr4 = (double *) arg4.data;
  DECLARE_PTR_ALIGNED(ptr4,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr5 = (double *) arg5.data;
  DECLARE_PTR_ALIGNED(ptr5,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr6 = (double *) arg6.data;
  DECLARE_PTR_ALIGNED(ptr6,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr7 = (double *) arg7.data;
  DECLARE_PTR_ALIGNED(ptr7,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr8 = (double *) arg8.data;
  DECLARE_PTR_ALIGNED(ptr8,double_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(9);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  backwards_euler_update_Q");
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
        backwards_euler_update_Q(
          &dat0[i],
          &(ptr1)[15 * (n+i)],
          &(ptr2)[15 * (n+i)],
          &(ptr3)[15 * (n+i)],
          &(ptr4)[15 * (n+i)],
          &(ptr5)[15 * (n+i)],
          &(ptr6)[15 * (n+i)],
          &(ptr7)[15 * (n+i)],
          &(ptr8)[15 * (n+i)]);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
    }
    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      backwards_euler_update_Q(
        (double*)arg0.data,
        &(ptr1)[15*n],
        &(ptr2)[15*n],
        &(ptr3)[15*n],
        &(ptr4)[15*n],
        &(ptr5)[15*n],
        &(ptr6)[15*n],
        &(ptr7)[15*n],
        &(ptr8)[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[9].name      = name;
  OP_kernels[9].count    += 1;
  OP_kernels[9].time     += wall_t2 - wall_t1;
  OP_kernels[9].transfer += (float)set->size * arg1.size;
  OP_kernels[9].transfer += (float)set->size * arg2.size;
  OP_kernels[9].transfer += (float)set->size * arg3.size;
  OP_kernels[9].transfer += (float)set->size * arg4.size;
  OP_kernels[9].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[9].transfer += (float)set->size * arg6.size * 2.0f;
  OP_kernels[9].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[9].transfer += (float)set->size * arg8.size * 2.0f;
}
