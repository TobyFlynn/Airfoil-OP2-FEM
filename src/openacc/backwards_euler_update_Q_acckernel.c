//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void backwards_euler_update_Q_openacc( const double *dt, const double *q0,
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

  double*arg0h = (double *)arg0.data;
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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(9);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[9].name      = name;
  OP_kernels[9].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  backwards_euler_update_Q");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg0_l = arg0h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    double* data5 = (double*)arg5.data_d;
    double* data6 = (double*)arg6.data_d;
    double* data7 = (double*)arg7.data_d;
    double* data8 = (double*)arg8.data_d;
    #pragma acc parallel loop independent deviceptr(data1,data2,data3,data4,data5,data6,data7,data8)
    for ( int n=0; n<set->size; n++ ){
      backwards_euler_update_Q_openacc(
        &arg0_l,
        &data1[15*n],
        &data2[15*n],
        &data3[15*n],
        &data4[15*n],
        &data5[15*n],
        &data6[15*n],
        &data7[15*n],
        &data8[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
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
