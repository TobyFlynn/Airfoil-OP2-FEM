//
// auto-generated by op2.py
//

//user function
#include <cmath>
#include "fluxes.h"

//user function
//#pragma acc routine
inline void get_bedge_q_openacc( const int *bedge_type, const int *bedgeNum,
                        const double *nx, const double *ny,
                        const double *fscale, const double *q, double *flux) {
  double exteriorQ[4 * 5];
  int exInd = 0;
  int nInd = 0;
  if(*bedgeNum == 1) {
    exInd = 4 * 5;
    nInd = 5;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 4 * 5;
    nInd = 2 * 5;
  }

  int *fmask;

  if(*bedgeNum == 0) {
    fmask = FMASK;
  } else if(*bedgeNum == 1) {
    fmask = &FMASK[5];
  } else {
    fmask = &FMASK[2 * 5];
  }

  if(*bedge_type == 0) {

    for(int i = 0; i < 5; i++) {
      exteriorQ[i * 4]     = bc_r;
      exteriorQ[i * 4 + 1] = bc_r * bc_u;
      exteriorQ[i * 4 + 2] =  bc_r * bc_v;
      exteriorQ[i * 4 + 3] = bc_e;
    }
  } else if(*bedge_type == 1) {

    for(int i = 0; i < 5; i++) {
      int qInd = fmask[i] * 4;
      exteriorQ[i * 4]     = bc_r;
      exteriorQ[i * 4 + 1] = bc_r * bc_u;
      exteriorQ[i * 4 + 2] =  bc_r * bc_v;
      exteriorQ[i * 4 + 3] = q[qInd + 3];
    }
  } else {

    for(int i = 0; i < 5; i++) {
      int qInd = fmask[i] * 4;
      exteriorQ[i * 4]     = q[qInd];
      exteriorQ[i * 4 + 1] = q[qInd + 1] - 2 * (nx[nInd + i] * q[qInd + 1] + ny[nInd + i] * q[qInd + 2]) * nx[nInd + i];
      exteriorQ[i * 4 + 2] = q[qInd + 2] - 2 * (nx[nInd + i] * q[qInd + 1] + ny[nInd + i] * q[qInd + 2]) * ny[nInd + i];
      exteriorQ[i * 4 + 3] = q[qInd + 3];
    }
  }


  roe(flux + exInd, nx + nInd, ny + nInd, fscale + nInd, q, exteriorQ, nInd);
}

// host stub function
void op_par_loop_get_bedge_q(char const *name, op_set set,
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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(5);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[5].name      = name;
  OP_kernels[5].count    += 1;

  int  ninds   = 5;
  int  inds[7] = {-1,-1,0,1,2,3,4};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: get_bedge_q\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_5
    int part_size = OP_PART_SIZE_5;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map2 = arg2.map_data_d;

    int* data0 = (int*)arg0.data_d;
    int* data1 = (int*)arg1.data_d;
    double *data2 = (double *)arg2.data_d;
    double *data3 = (double *)arg3.data_d;
    double *data4 = (double *)arg4.data_d;
    double *data5 = (double *)arg5.data_d;
    double *data6 = (double *)arg6.data_d;

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);
    ncolors = Plan->ncolors;
    int *col_reord = Plan->col_reord;
    int set_size1 = set->size + set->exec_size;

    // execute plan
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = Plan->col_offsets[0][col];
      int end = Plan->col_offsets[0][col+1];

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data0,data1,data2,data3,data4,data5,data6)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        map2idx = map2[n + set_size1 * 0];


        get_bedge_q_openacc(
          &data0[1 * n],
          &data1[1 * n],
          &data2[15 * map2idx],
          &data3[15 * map2idx],
          &data4[15 * map2idx],
          &data5[60 * map2idx],
          &data6[60 * map2idx]);
      }

    }
    OP_kernels[5].transfer  += Plan->transfer;
    OP_kernels[5].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[5].time     += wall_t2 - wall_t1;
}
