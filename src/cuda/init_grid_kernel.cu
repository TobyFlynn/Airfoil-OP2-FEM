//
// auto-generated by op2.py
//

#include <cblas.h>

//user function
__device__ void init_grid_gpu( const double *n0, const double *n1, const double *n2,
                      double *nodeX, double *nodeY, double *x, double *y,
                      double *xr, double *yr, double *xs, double *ys,
                      double *rx, double *ry, double *sx, double *sy,
                      double *nx, double *ny, double *fscale) {


  nodeX[0] = n0[0];
  nodeX[1] = n1[0];
  nodeX[2] = n2[0];
  nodeY[0] = n0[1];
  nodeY[1] = n1[1];
  nodeY[2] = n2[1];



  cblas_dcopy(15, ones_cuda, 1, x, 1);
  cblas_daxpy(15, 1.0, r_cuda, 1, x, 1);
  cblas_dscal(15, 0.5 * n1[0], x, 1);

  double temp[15];
  cblas_dcopy(15, ones_cuda, 1, temp, 1);
  cblas_daxpy(15, 1.0, s_cuda, 1, temp, 1);

  cblas_daxpy(15, 0.5 * n2[0], temp, 1, x, 1);

  cblas_dcopy(15, s_cuda, 1, temp, 1);
  cblas_daxpy(15, 1.0, r_cuda, 1, temp, 1);

  cblas_daxpy(15, -0.5 * n0[0], temp, 1, x, 1);



  cblas_dcopy(15, ones_cuda, 1, y, 1);
  cblas_daxpy(15, 1.0, r_cuda, 1, y, 1);
  cblas_dscal(15, 0.5 * n1[1], y, 1);

  cblas_dcopy(15, ones_cuda, 1, temp, 1);
  cblas_daxpy(15, 1.0, s_cuda, 1, temp, 1);

  cblas_daxpy(15, 0.5 * n2[1], temp, 1, y, 1);

  cblas_dcopy(15, s_cuda, 1, temp, 1);
  cblas_daxpy(15, 1.0, r_cuda, 1, temp, 1);

  cblas_daxpy(15, -0.5 * n0[1], temp, 1, y, 1);


  cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dr_cuda, 15, x, 1, 0.0, xr, 1);

  cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Ds_cuda, 15, x, 1, 0.0, xs, 1);

  cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dr_cuda, 15, y, 1, 0.0, yr, 1);

  cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Ds_cuda, 15, y, 1, 0.0, ys, 1);

  double J[15];
  for(int i = 0; i < 15; i++) {
    J[i] = -xs[i] * yr[i] + xr[i] * ys[i];
  }

  for(int i = 0; i < 15; i++) {
    rx[i] = ys[i] / J[i];
    sx[i] = -yr[i] / J[i];
    ry[i] = -xs[i] / J[i];
    sy[i] = xr[i] / J[i];
  }


  for(int i = 0; i < 5; i++) {
    nx[i] = yr[FMASK_cuda[i]];
    ny[i] = -xr[FMASK_cuda[i]];
  }

  for(int i = 0; i < 5; i++) {
    nx[5 + i] = ys[FMASK_cuda[5 + i]] - yr[FMASK_cuda[5 + i]];
    ny[5 + i] = xr[FMASK_cuda[5 + i]] - xs[FMASK_cuda[5 + i]];
  }

  for(int i = 0; i < 5; i++) {
    nx[2 * 5 + i] = -ys[FMASK_cuda[2 * 5 + i]];
    ny[2 * 5 + i] = xs[FMASK_cuda[2 * 5 + i]];
  }

  for(int i = 0; i < 3 * 5; i++) {
    double sJ = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    nx[i] = nx[i] / sJ;
    ny[i] = ny[i] / sJ;
    fscale[i] = sJ / J[FMASK_cuda[i]];
  }

}

// CUDA kernel function
__global__ void op_cuda_init_grid(
  const double *__restrict ind_arg0,
  const int *__restrict opDat0Map,
  double *arg3,
  double *arg4,
  double *arg5,
  double *arg6,
  double *arg7,
  double *arg8,
  double *arg9,
  double *arg10,
  double *arg11,
  double *arg12,
  double *arg13,
  double *arg14,
  double *arg15,
  double *arg16,
  double *arg17,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    int map0idx;
    int map1idx;
    int map2idx;
    map0idx = opDat0Map[n + set_size * 0];
    map1idx = opDat0Map[n + set_size * 1];
    map2idx = opDat0Map[n + set_size * 2];

    //user-supplied kernel call
    init_grid_gpu(ind_arg0+map0idx*2,
              ind_arg0+map1idx*2,
              ind_arg0+map2idx*2,
              arg3+n*3,
              arg4+n*3,
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
void op_par_loop_init_grid(char const *name, op_set set,
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
  op_timing_realloc(0);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[0].name      = name;
  OP_kernels[0].count    += 1;


  int    ninds   = 1;
  int    inds[18] = {0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: init_grid\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_0
      int nthread = OP_BLOCK_SIZE_0;
    #else
      int nthread = OP_block_size;
    #endif

    for ( int round=0; round<2; round++ ){
      if (round==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = round==0 ? 0 : set->core_size;
      int end = round==0 ? set->core_size : set->size + set->exec_size;
      if (end-start>0) {
        int nblocks = (end-start-1)/nthread+1;
        op_cuda_init_grid<<<nblocks,nthread>>>(
        (double *)arg0.data_d,
        arg0.map_data_d,
        (double*)arg3.data_d,
        (double*)arg4.data_d,
        (double*)arg5.data_d,
        (double*)arg6.data_d,
        (double*)arg7.data_d,
        (double*)arg8.data_d,
        (double*)arg9.data_d,
        (double*)arg10.data_d,
        (double*)arg11.data_d,
        (double*)arg12.data_d,
        (double*)arg13.data_d,
        (double*)arg14.data_d,
        (double*)arg15.data_d,
        (double*)arg16.data_d,
        (double*)arg17.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[0].time     += wall_t2 - wall_t1;
}
