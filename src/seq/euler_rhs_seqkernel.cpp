//
// auto-generated by op2.py
//

//user function
#include "../euler_rhs.h"

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


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  euler_rhs");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    for ( int n=0; n<set_size; n++ ){
      euler_rhs(
        &((double*)arg0.data)[60*n],
        &((double*)arg1.data)[60*n],
        &((double*)arg2.data)[15*n],
        &((double*)arg3.data)[15*n],
        &((double*)arg4.data)[15*n],
        &((double*)arg5.data)[15*n],
        &((double*)arg6.data)[15*n],
        &((double*)arg7.data)[15*n],
        &((double*)arg8.data)[15*n],
        &((double*)arg9.data)[60*n],
        &((double*)arg10.data)[60*n],
        &((double*)arg11.data)[60*n],
        &((double*)arg12.data)[60*n],
        &((double*)arg13.data)[60*n],
        &((double*)arg14.data)[60*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[7].name      = name;
  OP_kernels[7].count    += 1;
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
