//
// auto-generated by op2.py
//

//user function
inline void set_ic(double *q, double *workingQ) {
  for(int i = 0; i < 15; i++) {
    q[i * 4]     = bc_r;
    q[i * 4 + 1] = bc_r * bc_u;
    q[i * 4 + 2] = bc_r * bc_v;
    q[i * 4 + 3] = bc_e;
    workingQ[i * 4]     = q[i * 4];
    workingQ[i * 4 + 1] = q[i * 4 + 1];
    workingQ[i * 4 + 2] = q[i * 4 + 2];
    workingQ[i * 4 + 3] = q[i * 4 + 3];
  }
}

// host stub function
void op_par_loop_set_ic(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;
  //create aligned pointers for dats
  ALIGNED_double       double * __restrict__ ptr0 = (double *) arg0.data;
  DECLARE_PTR_ALIGNED(ptr0,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr1 = (double *) arg1.data;
  DECLARE_PTR_ALIGNED(ptr1,double_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(1);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  set_ic");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        set_ic(
          &(ptr0)[60 * (n+i)],
          &(ptr1)[60 * (n+i)]);
      }
    }
    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      set_ic(
        &(ptr0)[60*n],
        &(ptr1)[60*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[1].name      = name;
  OP_kernels[1].count    += 1;
  OP_kernels[1].time     += wall_t2 - wall_t1;
  OP_kernels[1].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[1].transfer += (float)set->size * arg1.size * 2.0f;
}
