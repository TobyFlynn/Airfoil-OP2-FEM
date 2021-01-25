//
// auto-generated by op2.py
//

//user function
#include "../init_grid.h"

// host stub function
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
  OP_kernels[0].name      = name;
  OP_kernels[0].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);

  int  ninds   = 1;
  int  inds[18] = {0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: init_grid\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_0
    int part_size = OP_PART_SIZE_0;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {

    op_plan *Plan = op_plan_get_stage_upload(name,set,part_size,nargs,args,ninds,inds,OP_STAGE_ALL,0);

    // execute plan
    int block_offset = 0;
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all(nargs, args);
      }
      int nblocks = Plan->ncolblk[col];

      #pragma omp parallel for
      for ( int blockIdx=0; blockIdx<nblocks; blockIdx++ ){
        int blockId  = Plan->blkmap[blockIdx + block_offset];
        int nelem    = Plan->nelems[blockId];
        int offset_b = Plan->offset[blockId];
        for ( int n=offset_b; n<offset_b+nelem; n++ ){
          int map0idx;
          int map1idx;
          int map2idx;
          map0idx = arg0.map_data[n * arg0.map->dim + 0];
          map1idx = arg0.map_data[n * arg0.map->dim + 1];
          map2idx = arg0.map_data[n * arg0.map->dim + 2];


          init_grid(
            &((double*)arg0.data)[2 * map0idx],
            &((double*)arg0.data)[2 * map1idx],
            &((double*)arg0.data)[2 * map2idx],
            &((double*)arg3.data)[3 * n],
            &((double*)arg4.data)[3 * n],
            &((double*)arg5.data)[15 * n],
            &((double*)arg6.data)[15 * n],
            &((double*)arg7.data)[15 * n],
            &((double*)arg8.data)[15 * n],
            &((double*)arg9.data)[15 * n],
            &((double*)arg10.data)[15 * n],
            &((double*)arg11.data)[15 * n],
            &((double*)arg12.data)[15 * n],
            &((double*)arg13.data)[15 * n],
            &((double*)arg14.data)[15 * n],
            &((double*)arg15.data)[15 * n],
            &((double*)arg16.data)[15 * n],
            &((double*)arg17.data)[15 * n]);
        }
      }

      block_offset += nblocks;
    }
    OP_kernels[0].transfer  += Plan->transfer;
    OP_kernels[0].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[0].time     += wall_t2 - wall_t1;
}
