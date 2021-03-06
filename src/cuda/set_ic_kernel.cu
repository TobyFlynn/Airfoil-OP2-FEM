//
// auto-generated by op2.py
//

//user function
__device__ void set_ic_gpu( double *q, double *workingQ) {
  for(int i = 0; i < 15; i++) {
    q[i * 4]     = bc_r_cuda;
    q[i * 4 + 1] = bc_r_cuda * bc_u_cuda;
    q[i * 4 + 2] = bc_r_cuda * bc_v_cuda;
    q[i * 4 + 3] = bc_e_cuda;
    workingQ[i * 4]     = q[i * 4];
    workingQ[i * 4 + 1] = q[i * 4 + 1];
    workingQ[i * 4 + 2] = q[i * 4 + 2];
    workingQ[i * 4 + 3] = q[i * 4 + 3];
  }

}

// CUDA kernel function
__global__ void op_cuda_set_ic(
  double *arg0,
  double *arg1,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    set_ic_gpu(arg0+n*60,
           arg1+n*60);
  }
}


//host stub function
void op_par_loop_set_ic(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(1);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[1].name      = name;
  OP_kernels[1].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  set_ic");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_1
      int nthread = OP_BLOCK_SIZE_1;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_set_ic<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[1].time     += wall_t2 - wall_t1;
  OP_kernels[1].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[1].transfer += (float)set->size * arg1.size * 2.0f;
}
