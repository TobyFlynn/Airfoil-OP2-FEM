//
// auto-generated by op2.py
//

#include <cmath>

//user function
__device__ void calc_dt_gpu( const double *q, const double *fscale, double *dt1) {
  double dt1_arr[3 * 5];
  for(int i = 0; i < 3 * 5; i++) {
    double rho = q[FMASK_cuda[i] * 4];
    double u = q[FMASK_cuda[i] * 4 + 1] / rho;
    double v = q[FMASK_cuda[i] * 4 + 2] / rho;
    double p = (gam_cuda - 1.0) * (q[FMASK_cuda[i] * 4 + 3] - rho * (u * u + v * v) * 0.5);
    double c = sqrt(abs(gam_cuda * p / rho));
    dt1_arr[i] = ((4 + 1) * (4 + 1)) * 0.5 * fscale[FMASK_cuda[i]] *(sqrt(u * u + v * v) + c);
  }

  double max = *dt1;
  for(int i = 0; i < 3 * 5; i++) {
    if(dt1_arr[i] > max) {
      max = dt1_arr[i];
    }
  }
  *dt1 = max;

}

// CUDA kernel function
__global__ void op_cuda_calc_dt(
  const double *__restrict arg0,
  const double *__restrict arg1,
  double *arg2,
  int   set_size ) {

  double arg2_l[1];
  for ( int d=0; d<1; d++ ){
    arg2_l[d]=arg2[d+blockIdx.x*1];
  }

  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    calc_dt_gpu(arg0+n*60,
            arg1+n*15,
            arg2_l);
  }

  //global reductions

  for ( int d=0; d<1; d++ ){
    op_reduction<OP_MAX>(&arg2[d+blockIdx.x*1],arg2_l[d]);
  }
}


//host stub function
void op_par_loop_calc_dt(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  double*arg2h = (double *)arg2.data;
  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(3);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[3].name      = name;
  OP_kernels[3].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  calc_dt");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_3
      int nthread = OP_BLOCK_SIZE_3;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    //transfer global reduction data to GPU
    int maxblocks = nblocks;
    int reduct_bytes = 0;
    int reduct_size  = 0;
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    reduct_size   = MAX(reduct_size,sizeof(double));
    reallocReductArrays(reduct_bytes);
    reduct_bytes = 0;
    arg2.data   = OP_reduct_h + reduct_bytes;
    arg2.data_d = OP_reduct_d + reduct_bytes;
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        ((double *)arg2.data)[d+b*1] = arg2h[d];
      }
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    mvReductArraysToDevice(reduct_bytes);

    int nshared = reduct_size*nthread;
    op_cuda_calc_dt<<<nblocks,nthread,nshared>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      set->size );
    //transfer global reduction data back to CPU
    mvReductArraysToHost(reduct_bytes);
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        arg2h[d] = MAX(arg2h[d],((double *)arg2.data)[d+b*1]);
      }
    }
    arg2.data = (char *)arg2h;
    op_mpi_reduce(&arg2,arg2h);
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[3].time     += wall_t2 - wall_t1;
  OP_kernels[3].transfer += (float)set->size * arg0.size;
  OP_kernels[3].transfer += (float)set->size * arg1.size;
}
