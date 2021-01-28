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
  op_arg arg14,
  op_arg arg15,
  op_arg arg16,
  op_arg arg17,
  op_arg arg18,
  op_arg arg19,
  op_arg arg20,
  op_arg arg21,
  op_arg arg22,
  op_arg arg23,
  op_arg arg24,
  op_arg arg25,
  op_arg arg26,
  op_arg arg27,
  op_arg arg28,
  op_arg arg29,
  op_arg arg30,
  op_arg arg31,
  op_arg arg32,
  op_arg arg33,
  op_arg arg34,
  op_arg arg35,
  op_arg arg36,
  op_arg arg37,
  op_arg arg38){

  int nargs = 39;
  op_arg args[39];

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
  args[21] = arg21;
  args[22] = arg22;
  args[23] = arg23;
  args[24] = arg24;
  args[25] = arg25;
  args[26] = arg26;
  args[27] = arg27;
  args[28] = arg28;
  args[29] = arg29;
  args[30] = arg30;
  args[31] = arg31;
  args[32] = arg32;
  args[33] = arg33;
  args[34] = arg34;
  args[35] = arg35;
  args[36] = arg36;
  args[37] = arg37;
  args[38] = arg38;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(6);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  euler_rhs");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    for ( int n=0; n<set_size; n++ ){
      euler_rhs(
        &((double*)arg0.data)[15*n],
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
        &((double*)arg20.data)[15*n],
        &((double*)arg21.data)[15*n],
        &((double*)arg22.data)[15*n],
        &((double*)arg23.data)[15*n],
        &((double*)arg24.data)[15*n],
        &((double*)arg25.data)[15*n],
        &((double*)arg26.data)[15*n],
        &((double*)arg27.data)[15*n],
        &((double*)arg28.data)[15*n],
        &((double*)arg29.data)[15*n],
        &((double*)arg30.data)[15*n],
        &((double*)arg31.data)[15*n],
        &((double*)arg32.data)[15*n],
        &((double*)arg33.data)[15*n],
        &((double*)arg34.data)[15*n],
        &((double*)arg35.data)[15*n],
        &((double*)arg36.data)[15*n],
        &((double*)arg37.data)[15*n],
        &((double*)arg38.data)[15*n]);
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
  OP_kernels[6].transfer += (float)set->size * arg1.size;
  OP_kernels[6].transfer += (float)set->size * arg2.size;
  OP_kernels[6].transfer += (float)set->size * arg3.size;
  OP_kernels[6].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg6.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg8.size;
  OP_kernels[6].transfer += (float)set->size * arg9.size;
  OP_kernels[6].transfer += (float)set->size * arg10.size;
  OP_kernels[6].transfer += (float)set->size * arg11.size;
  OP_kernels[6].transfer += (float)set->size * arg12.size;
  OP_kernels[6].transfer += (float)set->size * arg13.size;
  OP_kernels[6].transfer += (float)set->size * arg14.size;
  OP_kernels[6].transfer += (float)set->size * arg15.size;
  OP_kernels[6].transfer += (float)set->size * arg16.size;
  OP_kernels[6].transfer += (float)set->size * arg17.size;
  OP_kernels[6].transfer += (float)set->size * arg18.size;
  OP_kernels[6].transfer += (float)set->size * arg19.size;
  OP_kernels[6].transfer += (float)set->size * arg20.size;
  OP_kernels[6].transfer += (float)set->size * arg21.size;
  OP_kernels[6].transfer += (float)set->size * arg22.size;
  OP_kernels[6].transfer += (float)set->size * arg23.size;
  OP_kernels[6].transfer += (float)set->size * arg24.size;
  OP_kernels[6].transfer += (float)set->size * arg25.size;
  OP_kernels[6].transfer += (float)set->size * arg26.size;
  OP_kernels[6].transfer += (float)set->size * arg27.size;
  OP_kernels[6].transfer += (float)set->size * arg28.size;
  OP_kernels[6].transfer += (float)set->size * arg29.size;
  OP_kernels[6].transfer += (float)set->size * arg30.size;
  OP_kernels[6].transfer += (float)set->size * arg31.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg32.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg33.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg34.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg35.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg36.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg37.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg38.size * 2.0f;
}
