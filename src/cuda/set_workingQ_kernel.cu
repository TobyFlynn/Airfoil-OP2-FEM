//
// auto-generated by op2.py
//

//user function
__device__ void set_workingQ_gpu( const double *dt, const int *stage, const double *q0,
                         const double *q1, const double *q2, const double *q3,
                         const double *k10, const double *k11,
                         const double *k12, const double *k13,
                         const double *k20, const double *k21,
                         const double *k22, const double *k23,
                         double *workingQ0, double *workingQ1,
                         double *workingQ2, double *workingQ3) {
  if(*stage == 0) {
    for(int i = 0; i < 15; i++) {
      workingQ0[i] = q0[i] + (*dt) * k10[i];
      workingQ1[i] = q1[i] + (*dt) * k11[i];
      workingQ2[i] = q2[i] + (*dt) * k12[i];
      workingQ3[i] = q3[i] + (*dt) * k13[i];
    }
  } else {
    for(int i = 0; i < 15; i++) {
      workingQ0[i] = q0[i] + (*dt) * (k10[i] / 4.0 + k20[i] / 4.0);
      workingQ1[i] = q1[i] + (*dt) * (k11[i] / 4.0 + k21[i] / 4.0);
      workingQ2[i] = q2[i] + (*dt) * (k12[i] / 4.0 + k22[i] / 4.0);
      workingQ3[i] = q3[i] + (*dt) * (k13[i] / 4.0 + k23[i] / 4.0);
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_set_workingQ(
  const double *arg0,
  const int *arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  const double *__restrict arg5,
  const double *__restrict arg6,
  const double *__restrict arg7,
  const double *__restrict arg8,
  const double *__restrict arg9,
  const double *__restrict arg10,
  const double *__restrict arg11,
  const double *__restrict arg12,
  const double *__restrict arg13,
  double *arg14,
  double *arg15,
  double *arg16,
  double *arg17,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    set_workingQ_gpu(arg0,
                 arg1,
                 arg2+n*15,
                 arg3+n*15,
                 arg4+n*15,
                 arg5+n*15,
                 arg6+n*15,
                 arg7+n*15,
                 arg8+n*15,
                 arg9+n*15,
                 arg10+n*15,
                 arg11+n*15,
                 arg12+n*15,
                 arg13+n*15,
                 arg14+n*15,
                 arg15+n*15,
                 arg16+n*15,
                 arg17+n*15);
  }
}


//host stub function
void op_par_loop_set_workingQ(char const *name, op_set set,
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
  op_arg arg17){

  double*arg0h = (double *)arg0.data;
  int*arg1h = (int *)arg1.data;
  int nargs = 18;
  op_arg args[18];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(7);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[7].name      = name;
  OP_kernels[7].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  set_workingQ");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(double));
    consts_bytes += ROUND_UP(1*sizeof(int));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg0.data   = OP_consts_h + consts_bytes;
    arg0.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg0.data)[d] = arg0h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    arg1.data   = OP_consts_h + consts_bytes;
    arg1.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((int *)arg1.data)[d] = arg1h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_7
      int nthread = OP_BLOCK_SIZE_7;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_set_workingQ<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (int *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      (double *) arg5.data_d,
      (double *) arg6.data_d,
      (double *) arg7.data_d,
      (double *) arg8.data_d,
      (double *) arg9.data_d,
      (double *) arg10.data_d,
      (double *) arg11.data_d,
      (double *) arg12.data_d,
      (double *) arg13.data_d,
      (double *) arg14.data_d,
      (double *) arg15.data_d,
      (double *) arg16.data_d,
      (double *) arg17.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[7].time     += wall_t2 - wall_t1;
  OP_kernels[7].transfer += (float)set->size * arg2.size;
  OP_kernels[7].transfer += (float)set->size * arg3.size;
  OP_kernels[7].transfer += (float)set->size * arg4.size;
  OP_kernels[7].transfer += (float)set->size * arg5.size;
  OP_kernels[7].transfer += (float)set->size * arg6.size;
  OP_kernels[7].transfer += (float)set->size * arg7.size;
  OP_kernels[7].transfer += (float)set->size * arg8.size;
  OP_kernels[7].transfer += (float)set->size * arg9.size;
  OP_kernels[7].transfer += (float)set->size * arg10.size;
  OP_kernels[7].transfer += (float)set->size * arg11.size;
  OP_kernels[7].transfer += (float)set->size * arg12.size;
  OP_kernels[7].transfer += (float)set->size * arg13.size;
  OP_kernels[7].transfer += (float)set->size * arg14.size * 2.0f;
  OP_kernels[7].transfer += (float)set->size * arg15.size * 2.0f;
  OP_kernels[7].transfer += (float)set->size * arg16.size * 2.0f;
  OP_kernels[7].transfer += (float)set->size * arg17.size * 2.0f;
}
