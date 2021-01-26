//
// auto-generated by op2.py
//

// #include <cblas.h>
#include <algorithm>
#include <cmath>

//user function
__device__ void euler_rhs_gpu( const double *q, double *exteriorQ,
                      const double *rx, const double *ry, const double *sx,
                      const double *sy, const double *fscale, const double *nx,
                      const double *ny, double *qRHS) {
  double F[4 * 15];
  double G[4 * 15];
  for(int i = 0; i < 15; i++) {
    double rho, u, v, p;
    euler_flux(&q[i * 4], &F[i * 4], &G[i * 4], &rho, &u, &v, &p);
  }

  for(int i = 0; i < 4; i++) {
    double dFdr[15];
    double dFds[15];
    double dGdr[15];
    double dGds[15];





    for(int j = 0; j < 15; j++) {

    }
  }

  double mQ[4 * 3 * 5];
  double mF[4 * 3 * 5];
  double mG[4 * 3 * 5];
  double mRho[3 * 5];
  double mU[3 * 5];
  double mV[3 * 5];
  double mP[3 * 5];

  for(int i = 0; i < 3 * 5; i++) {
    int ind = FMASK_cuda[i] * 4;
    mQ[i * 4]     = q[ind];
    mQ[i * 4 + 1] = q[ind + 1];
    mQ[i * 4 + 2] = q[ind + 2];
    mQ[i * 4 + 3] = q[ind + 3];

    euler_flux(&mQ[i * 4], &mF[i * 4], &mG[i * 4], &mRho[i], &mU[i], &mV[i], &mP[i]);
  }

  double pF[4 * 3 * 5];
  double pG[4 * 3 * 5];
  double pRho[3 * 5];
  double pU[3 * 5];
  double pV[3 * 5];
  double pP[3 * 5];
  for(int i = 0; i < 3 * 5; i++) {
    euler_flux(&exteriorQ[i * 4], &pF[i * 4], &pG[i * 4], &pRho[i], &pU[i], &pV[i], &pP[i]);
  }

  double flux[4 * 3 * 5];

  roe(flux, nx, ny, fscale, q, exteriorQ);

  for(int i = 0; i < 4; i++) {

  }

  for(int i = 0; i < 4 * 3 * 5; i++) {
    exteriorQ[i] = 0.0;
  }

}

// CUDA kernel function
__global__ void op_cuda_euler_rhs(
  const double *__restrict arg0,
  double *arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  const double *__restrict arg5,
  const double *__restrict arg6,
  const double *__restrict arg7,
  const double *__restrict arg8,
  double *arg9,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    euler_rhs_gpu(arg0+n*60,
              arg1+n*60,
              arg2+n*15,
              arg3+n*15,
              arg4+n*15,
              arg5+n*15,
              arg6+n*15,
              arg7+n*15,
              arg8+n*15,
              arg9+n*60);
  }
}


//host stub function
void op_par_loop_euler_rhs(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  int nargs = 10;
  op_arg args[10];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(6);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[6].name      = name;
  OP_kernels[6].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  euler_rhs");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_6
      int nthread = OP_BLOCK_SIZE_6;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_euler_rhs<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      (double *) arg5.data_d,
      (double *) arg6.data_d,
      (double *) arg7.data_d,
      (double *) arg8.data_d,
      (double *) arg9.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[6].time     += wall_t2 - wall_t1;
  OP_kernels[6].transfer += (float)set->size * arg0.size;
  OP_kernels[6].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg2.size;
  OP_kernels[6].transfer += (float)set->size * arg3.size;
  OP_kernels[6].transfer += (float)set->size * arg4.size;
  OP_kernels[6].transfer += (float)set->size * arg5.size;
  OP_kernels[6].transfer += (float)set->size * arg6.size;
  OP_kernels[6].transfer += (float)set->size * arg7.size;
  OP_kernels[6].transfer += (float)set->size * arg8.size;
  OP_kernels[6].transfer += (float)set->size * arg9.size * 2.0f;
}
