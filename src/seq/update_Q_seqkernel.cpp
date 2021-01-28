//
// auto-generated by op2.py
//

//user function
#include "../update_Q.h"

// host stub function
void op_par_loop_update_Q(char const *name, op_set set,
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
  op_arg arg17,
  op_arg arg18,
  op_arg arg19,
  op_arg arg20){

  int nargs = 21;
  op_arg args[21];

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
  args[18] = arg18;
  args[19] = arg19;
  args[20] = arg20;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(8);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  update_Q");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    for ( int n=0; n<set_size; n++ ){
      update_Q(
        (double*)arg0.data,
        &((double*)arg1.data)[15*n],
        &((double*)arg2.data)[15*n],
        &((double*)arg3.data)[15*n],
        &((double*)arg4.data)[15*n],
        &((double*)arg5.data)[15*n],
        &((double*)arg6.data)[15*n],
        &((double*)arg7.data)[15*n],
        &((double*)arg8.data)[15*n],
        &((double*)arg9.data)[15*n],
        &((double*)arg10.data)[15*n],
        &((double*)arg11.data)[15*n],
        &((double*)arg12.data)[15*n],
        &((double*)arg13.data)[15*n],
        &((double*)arg14.data)[15*n],
        &((double*)arg15.data)[15*n],
        &((double*)arg16.data)[15*n],
        &((double*)arg17.data)[15*n],
        &((double*)arg18.data)[15*n],
        &((double*)arg19.data)[15*n],
        &((double*)arg20.data)[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[8].name      = name;
  OP_kernels[8].count    += 1;
  OP_kernels[8].time     += wall_t2 - wall_t1;
  OP_kernels[8].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[8].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[8].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[8].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[8].transfer += (float)set->size * arg5.size;
  OP_kernels[8].transfer += (float)set->size * arg6.size;
  OP_kernels[8].transfer += (float)set->size * arg7.size;
  OP_kernels[8].transfer += (float)set->size * arg8.size;
  OP_kernels[8].transfer += (float)set->size * arg9.size;
  OP_kernels[8].transfer += (float)set->size * arg10.size;
  OP_kernels[8].transfer += (float)set->size * arg11.size;
  OP_kernels[8].transfer += (float)set->size * arg12.size;
  OP_kernels[8].transfer += (float)set->size * arg13.size;
  OP_kernels[8].transfer += (float)set->size * arg14.size;
  OP_kernels[8].transfer += (float)set->size * arg15.size;
  OP_kernels[8].transfer += (float)set->size * arg16.size;
  OP_kernels[8].transfer += (float)set->size * arg17.size * 2.0f;
  OP_kernels[8].transfer += (float)set->size * arg18.size * 2.0f;
  OP_kernels[8].transfer += (float)set->size * arg19.size * 2.0f;
  OP_kernels[8].transfer += (float)set->size * arg20.size * 2.0f;
}
