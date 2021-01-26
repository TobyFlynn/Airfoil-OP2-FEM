//
// auto-generated by op2.py
//

//user function
#include "fluxes.h"

//user function
//#pragma acc routine
inline void euler_rhs_openacc( const double *q, double *exteriorQ,
                      const double *rx, const double *ry, const double *sx,
                      const double *sy, const double *fscale, const double *nx,
                      const double *ny, const double *dFdr, const double *dFds,
                      const double *dGdr, const double *dGds, double *flux,
                      double *qRHS) {

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 15; j++) {
      qRHS[i + j * 4] = (rx[j] * dFdr[i + j * 4] + sx[j] * dFds[i + j * 4]) + (ry[j] * dGdr[i + j * 4] + sy[j] * dGds[i + j * 4]);
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
    int ind = FMASK[i] * 4;
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

  roe(flux, nx, ny, fscale, q, exteriorQ);

  for(int i = 0; i < 4 * 3 * 5; i++) {
    exteriorQ[i] = 0.0;
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
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14){

  int nargs = 15;
  op_arg args[15];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(7);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[7].name      = name;
  OP_kernels[7].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  euler_rhs");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    double* data5 = (double*)arg5.data_d;
    double* data6 = (double*)arg6.data_d;
    double* data7 = (double*)arg7.data_d;
    double* data8 = (double*)arg8.data_d;
    double* data9 = (double*)arg9.data_d;
    double* data10 = (double*)arg10.data_d;
    double* data11 = (double*)arg11.data_d;
    double* data12 = (double*)arg12.data_d;
    double* data13 = (double*)arg13.data_d;
    double* data14 = (double*)arg14.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14)
    for ( int n=0; n<set->size; n++ ){
      euler_rhs_openacc(
        &data0[60*n],
        &data1[60*n],
        &data2[15*n],
        &data3[15*n],
        &data4[15*n],
        &data5[15*n],
        &data6[15*n],
        &data7[15*n],
        &data8[15*n],
        &data9[60*n],
        &data10[60*n],
        &data11[60*n],
        &data12[60*n],
        &data13[60*n],
        &data14[60*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[7].time     += wall_t2 - wall_t1;
  OP_kernels[7].transfer += (float)set->size * arg0.size;
  OP_kernels[7].transfer += (float)set->size * arg1.size * 2.0f;
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
  OP_kernels[7].transfer += (float)set->size * arg13.size * 2.0f;
  OP_kernels[7].transfer += (float)set->size * arg14.size * 2.0f;
}
