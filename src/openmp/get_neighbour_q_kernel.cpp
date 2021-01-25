//
// auto-generated by op2.py
//

//user function
#include "../get_neighbour_q.h"

// host stub function
void op_par_loop_get_neighbour_q(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

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
  op_timing_realloc(4);
  OP_kernels[4].name      = name;
  OP_kernels[4].count    += 1;
  op_timers_core(&cpu_t1, &wall_t1);

  int  ninds   = 4;
  int  inds[9] = {-1,0,1,0,1,2,2,3,3};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: get_neighbour_q\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_4
    int part_size = OP_PART_SIZE_4;
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
          int map1idx;
          int map3idx;
          map1idx = arg1.map_data[n * arg1.map->dim + 0];
          map3idx = arg1.map_data[n * arg1.map->dim + 1];


          get_neighbour_q(
            &((int*)arg0.data)[2 * n],
            &((double*)arg1.data)[3 * map1idx],
            &((double*)arg2.data)[3 * map1idx],
            &((double*)arg1.data)[3 * map3idx],
            &((double*)arg2.data)[3 * map3idx],
            &((double*)arg5.data)[60 * map1idx],
            &((double*)arg5.data)[60 * map3idx],
            &((double*)arg7.data)[60 * map1idx],
            &((double*)arg7.data)[60 * map3idx]);
        }
      }

      block_offset += nblocks;
    }
    OP_kernels[4].transfer  += Plan->transfer;
    OP_kernels[4].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[4].time     += wall_t2 - wall_t1;
}
