//
// auto-generated by op2.py
//

//user function
#include <cblas.h>
#include <cmath>

#include "fluxes.h"

inline void euler_rhs(const double *q, double *flux,
                      const double *rx, const double *ry, const double *sx,
                      const double *sy, double *qRHS) {
  double F[4 * 15];
  double G[4 * 15];
  for(int i = 0; i < 15; i++) {
    double rho, u, v, p;
    euler_flux(&q[i * 4], &F[i * 4], &G[i * 4], &rho, &u, &v, &p);
  }

  // Compute weak derivatives
  for(int i = 0; i < 4; i++) {
    double dFdr[15];
    double dFds[15];
    double dGdr[15];
    double dGds[15];

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Drw, 15, &F[i], 4, 0.0, dFdr, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dsw, 15, &F[i], 4, 0.0, dFds, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Drw, 15, &G[i], 4, 0.0, dGdr, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dsw, 15, &G[i], 4, 0.0, dGds, 1);

    for(int j = 0; j < 15; j++) {
      qRHS[i + j * 4] = (rx[j] * dFdr[j] + sx[j] * dFds[j]) + (ry[j] * dGdr[j] + sy[j] * dGds[j]);
    }
  }

  for(int i = 0; i < 4; i++) {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, -1.0, LIFT, 15, &flux[i], 4, 1.0, qRHS + i, 4);
  }

  for(int i = 0; i < 4 * 3 * 5; i++) {
    flux[i] = 0.0;
  }
}

// host stub function
void op_par_loop_euler_rhs(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6){

  int nargs = 7;
  op_arg args[7];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  //create aligned pointers for dats
  ALIGNED_double const double * __restrict__ ptr0 = (double *) arg0.data;
  DECLARE_PTR_ALIGNED(ptr0,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr1 = (double *) arg1.data;
  DECLARE_PTR_ALIGNED(ptr1,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr2 = (double *) arg2.data;
  DECLARE_PTR_ALIGNED(ptr2,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr3 = (double *) arg3.data;
  DECLARE_PTR_ALIGNED(ptr3,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr4 = (double *) arg4.data;
  DECLARE_PTR_ALIGNED(ptr4,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr5 = (double *) arg5.data;
  DECLARE_PTR_ALIGNED(ptr5,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr6 = (double *) arg6.data;
  DECLARE_PTR_ALIGNED(ptr6,double_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(6);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  euler_rhs");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        euler_rhs(
          &(ptr0)[60 * (n+i)],
          &(ptr1)[60 * (n+i)],
          &(ptr2)[15 * (n+i)],
          &(ptr3)[15 * (n+i)],
          &(ptr4)[15 * (n+i)],
          &(ptr5)[15 * (n+i)],
          &(ptr6)[60 * (n+i)]);
      }
    }
    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      euler_rhs(
        &(ptr0)[60*n],
        &(ptr1)[60*n],
        &(ptr2)[15*n],
        &(ptr3)[15*n],
        &(ptr4)[15*n],
        &(ptr5)[15*n],
        &(ptr6)[60*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[6].name      = name;
  OP_kernels[6].count    += 1;
  OP_kernels[6].time     += wall_t2 - wall_t1;
  OP_kernels[6].transfer += (float)set->size * arg0.size;
  OP_kernels[6].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg2.size;
  OP_kernels[6].transfer += (float)set->size * arg3.size;
  OP_kernels[6].transfer += (float)set->size * arg4.size;
  OP_kernels[6].transfer += (float)set->size * arg5.size;
  OP_kernels[6].transfer += (float)set->size * arg6.size * 2.0f;
}
