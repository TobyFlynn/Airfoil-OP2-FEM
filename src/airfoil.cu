// Include OP2 stuff
#include "op_seq.h"
// Include C++ stuff
#include <string>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <getopt.h>

#include "petscksp.h"
#include "petscvec.h"
#include "cublas_v2.h"

#include "constants/constant_r.h"
#include "constants/constant_s.h"
#include "constants/constant_Dr.h"
#include "constants/constant_Ds.h"
#include "constants/constant_Drw.h"
#include "constants/constant_Dsw.h"
#include "constants/constant_LIFT.h"
#include "constants/constant_util.h"
#include "load_mesh.h"
#include "save_solution.h"
#include "airfoil_cuda_matrices.h"
#include "airfoil_data.h"

// Include kernels
// #include "init_grid.h"
// #include "euler_rhs.h"
// #include "get_neighbour_q.h"
// #include "get_bedge_q.h"
// #include "set_ic.h"
// #include "set_workingQ.h"
// #include "update_Q.h"
// #include "calc_dt.h"
// #include "internal_fluxes.h"

using namespace std;

// Stuff for parsing command line arguments
extern char *optarg;
extern int  optind, opterr, optopt;
static struct option options[] = {
  {"iter", required_argument, 0, 0},
  {"alpha", required_argument, 0, 0},
  {0,    0,                  0,  0}
};

PetscErrorCode matAMult(Mat A, Vec x, Vec y);

inline void get_RHS(op_dat *Q, op_dat *rhs);

// Global OP2 sets
op_set nodes, cells, edges, bedges;

// Global OP2 maps
op_map cell2nodes, edge2nodes, edge2cells, bedge2nodes, bedge2cells;

// Global OP2 datasets
op_dat node_coords, nodeX, nodeY, x, y, xr, yr, xs, ys, rx, ry, sx, sy, nx, ny, fscale, bedge_type, edgeNum, bedgeNum;
op_dat rhs[4], Q[4], F[4], G[4], dFdr[4], dFds[4], dGdr[4], dGds[4], workingQ[4], exteriorQ[4], flux[4], rk[3][4];

cublasHandle_t cublas_handle;

int main(int argc, char **argv) {
  // Initialise cuBLAS
  cublasCreate(&cublas_handle);
  cublasSetPointerMode(cublas_handle, CUBLAS_POINTER_MODE_HOST);

  // Initialise PETSc
  char help[] = "TODO";
  int ierr = PetscInitialize(&argc, &argv, (char *)0, help);
  if(ierr) {
    cout << "Error initialising PETSc" << endl;
    return ierr;
  }

  double cpu_1, wall_1, cpu_2, wall_2;

  op_timers(&cpu_1, &wall_1);

  double *coords;
  int *cgnsCells;
  int *edge2node_data;
  int *edge2cell_data;
  int *bedge2node_data;
  int *bedge2cell_data;
  int *bedge_type_data;
  int *edgeNum_data;
  int *bedgeNum_data;
  int numNodes, numCells, numEdges, numBoundaryEdges;

  load_mesh("./naca0012.cgns", &numNodes, &numCells, &numEdges,
            &numBoundaryEdges, &coords, &cgnsCells, &edge2node_data,
            &edge2cell_data, &edgeNum_data, &bedge2node_data, &bedge2cell_data,
            &bedgeNum_data, &bedge_type_data);

  // Initialise OP2
  op_init(argc, argv, 2);

  // Get input from args
  int iter = 1;
  bc_alpha = 0.0;

  int opt_index = 0;
  while(getopt_long_only(argc, argv, "", options, &opt_index) != -1) {
    if(strcmp((char*)options[opt_index].name,"iter") == 0) iter = atoi(optarg);
    if(strcmp((char*)options[opt_index].name,"alpha") == 0) bc_alpha = stod(optarg);
  }

  if(bc_alpha != 0.0) {
    gam = 1.4;
    bc_mach = 0.4;
    // bc_mach = 0.1;
    bc_p = 1.0;
    bc_r = 1.0;
    bc_u = sin((M_PI/2.0) - (bc_alpha * M_PI / 180.0)) * sqrt(gam * bc_p / bc_r) * bc_mach;
    bc_v = cos((M_PI/2.0) - (bc_alpha * M_PI / 180.0)) * sqrt(gam * bc_p / bc_r) * bc_mach;
    //bc_e = bc_p / (bc_r * (gam - 1.0)) + 0.5f * bc_u * bc_u;
    bc_e = bc_p / (gam - 1.0) + 0.5 * bc_r * (bc_u * bc_u + bc_v * bc_v);
  } else {
    gam = 1.4;
    bc_mach = 0.4;
    // bc_mach = 0.1;
    bc_p = 1.0;
    bc_r = 1.0;
    bc_u = sqrt(gam * bc_p / bc_r) * bc_mach;
    bc_v = 0.0;
    //bc_e = bc_p / (bc_r * (gam - 1.0)) + 0.5f * bc_u * bc_u;
    bc_e = bc_p / (gam - 1.0) + 0.5 * bc_r * (bc_u * bc_u);
  }

  cout << "gam: " << gam << endl;
  cout << "bc_mach: " << bc_mach << endl;
  cout << "bc_alpha: " << bc_alpha << endl;
  cout << "bc_p: " << bc_p << endl;
  cout << "bc_r: " << bc_r << endl;
  cout << "bc_u: " << bc_u << endl;
  cout << "bc_v: " << bc_v << endl;
  cout << "bc_e: " << bc_e << endl;

  // Declare memory for data that will be calculated in initialisation kernel
  AirfoilData *data = new AirfoilData(numCells);

  // Declare OP2 sets
  nodes  = op_decl_set(numNodes, "nodes");
  cells  = op_decl_set(numCells, "cells");
  edges  = op_decl_set(numEdges, "edges");
  bedges = op_decl_set(numBoundaryEdges, "bedges");

  // Declare OP2 maps
  cell2nodes  = op_decl_map(cells, nodes, 3, cgnsCells, "cell2nodes");
  edge2nodes  = op_decl_map(edges, nodes, 2, edge2node_data, "edge2nodes");
  edge2cells  = op_decl_map(edges, cells, 2, edge2cell_data, "edge2cells");
  bedge2nodes = op_decl_map(bedges, nodes, 2, bedge2node_data, "bedge2nodes");
  bedge2cells = op_decl_map(bedges, cells, 1, bedge2cell_data, "bedge2cells");

  // Declare OP2 datasets
    // Structure: {x, y}
  node_coords = op_decl_dat(nodes, 2, "double", coords, "node_coords");
    // Coords of nodes per cell
  nodeX = op_decl_dat(cells, 3, "double", data->nodeX_data, "nodeX");
  nodeY = op_decl_dat(cells, 3, "double", data->nodeY_data, "nodeY");
    // The x and y coordinates of all the solution points in a cell
  x = op_decl_dat(cells, 15, "double", data->x_data, "x");
  y = op_decl_dat(cells, 15, "double", data->y_data, "y");
    // Geometric factors that relate to mapping between global and local (cell) coordinates
  xr = op_decl_dat(cells, 15, "double", data->xr_data, "xr");
  yr = op_decl_dat(cells, 15, "double", data->yr_data, "yr");
  xs = op_decl_dat(cells, 15, "double", data->xs_data, "xs");
  ys = op_decl_dat(cells, 15, "double", data->ys_data, "ys");
  rx = op_decl_dat(cells, 15, "double", data->rx_data, "rx");
  ry = op_decl_dat(cells, 15, "double", data->ry_data, "ry");
  sx = op_decl_dat(cells, 15, "double", data->sx_data, "sx");
  sy = op_decl_dat(cells, 15, "double", data->sy_data, "sy");
    // Normals for each cell (calculated for each node on each edge, nodes can appear on multiple edges)
  nx = op_decl_dat(cells, 3 * 5, "double", data->nx_data, "nx");
  ny = op_decl_dat(cells, 3 * 5, "double", data->ny_data, "ny");
    // surface Jacobian / Jacobian (used when lifting the boundary fluxes)
  fscale = op_decl_dat(cells, 3 * 5, "double", data->fscale_data, "fscale");
    // Values for compressible Euler equations in vectors
    // Structure: {q0_0, q1_0, q2_0, q3_0, q0_1, q1_1, ..., q3_15}
  for(int i = 0; i < 4; i++) {
    // Q[i] = op_decl_dat(cells, 15, "double", data->Q_data[i], "Q" + i);
    F[i] = op_decl_dat(cells, 15, "double", data->F_data[i], "F" + i);
    G[i] = op_decl_dat(cells, 15, "double", data->G_data[i], "G" + i);
    dFdr[i] = op_decl_dat(cells, 15, "double", data->dFdr_data[i], "dFdr" + i);
    dFds[i] = op_decl_dat(cells, 15, "double", data->dFds_data[i], "dFds" + i);
    dGdr[i] = op_decl_dat(cells, 15, "double", data->dGdr_data[i], "dGdr" + i);
    dGds[i] = op_decl_dat(cells, 15, "double", data->dGds_data[i], "dGds" + i);
    workingQ[i] = op_decl_dat(cells, 15, "double", data->workingQ_data[i], "workingQ" + i);
    exteriorQ[i] = op_decl_dat(cells, 3 * 5, "double", data->exteriorQ_data[i], "exteriorQ" + i);
    flux[i] = op_decl_dat(cells, 3 * 5, "double", data->flux_data[i], "flux" + i);
    rk[0][i] = op_decl_dat(cells, 15, "double", data->rk1_data[i], "rk1" + i);
    rk[1][i] = op_decl_dat(cells, 15, "double", data->rk2_data[i], "rk2" + i);
    rk[2][i] = op_decl_dat(cells, 15, "double", data->rk3_data[i], "rk3" + i);
    Q[i] = op_decl_dat(cells, 15, "double", data->Q_data[i], "Q" + i);
    rhs[i] = op_decl_dat(cells, 15, "double", data->rhs_data[i], "rhs" + i);
  }
  bedge_type = op_decl_dat(bedges, 1, "int", bedge_type_data, "bedge_type");
  edgeNum = op_decl_dat(edges, 2, "int", edgeNum_data, "edgeNum");
  bedgeNum  = op_decl_dat(bedges, 1, "int", bedgeNum_data, "bedgeNum");

  // Declare OP2 constants
  op_decl_const(1, "double", &gam);
  op_decl_const(1, "double", &bc_mach);
  op_decl_const(1, "double", &bc_alpha);
  op_decl_const(1, "double", &bc_p);
  op_decl_const(1, "double", &bc_r);
  op_decl_const(1, "double", &bc_u);
  op_decl_const(1, "double", &bc_v);
  op_decl_const(1, "double", &bc_e);
  op_decl_const(15, "double", ones);
  op_decl_const(15, "double", r);
  op_decl_const(15, "double", s);
  op_decl_const(15 * 15, "double", Dr);
  op_decl_const(15 * 15, "double", Ds);
  op_decl_const(15 * 15, "double", Drw);
  op_decl_const(15 * 15, "double", Dsw);
  op_decl_const(3 * 5, "int", FMASK);
  op_decl_const(15 * 15, "double", LIFT);

  // Matrix multiplications using cuBLAS
  op_arg init_grid_args[] = {
    op_arg_dat(x, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(y, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(xr, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(xs, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(yr, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(ys, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(cells, 6, init_grid_args);
  init_grid_matrices(cublas_handle, numCells, (double *)node_coords->data,
                     (int *)cell2nodes->map, (double *)x->data_d,
                     (double *)y->data_d, (double *)xr->data_d,
                     (double *)xs->data_d, (double *)yr->data_d,
                     (double *)ys->data_d);
  // Check this
  op_mpi_set_dirtybit_cuda(6, init_grid_args);
  // x->dirty_hd = 2;
  // y->dirty_hd = 2;
  // xr->dirty_hd = 2;
  // xs->dirty_hd = 2;
  // yr->dirty_hd = 2;
  // ys->dirty_hd = 2;

  // Initialisation kernels
  op_par_loop(init_grid, "init_grid", cells,
              op_arg_dat(node_coords, 0, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(node_coords, 1, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(node_coords, 2, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(nodeX, -1, OP_ID, 3, "double", OP_WRITE),
              op_arg_dat(nodeY, -1, OP_ID, 3, "double", OP_WRITE),
              op_arg_dat(xr, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(yr, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(xs, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(ys, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(rx, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(ry, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(sx, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(sy, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(nx, -1, OP_ID, 3 * 5, "double", OP_WRITE),
              op_arg_dat(ny, -1, OP_ID, 3 * 5, "double", OP_WRITE),
              op_arg_dat(fscale, -1, OP_ID, 3 * 5, "double", OP_WRITE));


  op_par_loop(set_ic, "set_ic", cells,
              op_arg_dat(Q[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(Q[1], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(Q[2], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(Q[3], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(workingQ[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(workingQ[1], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(workingQ[2], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(workingQ[3], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(exteriorQ[0], -1, OP_ID, 3 * 5, "double", OP_WRITE),
              op_arg_dat(exteriorQ[1], -1, OP_ID, 3 * 5, "double", OP_WRITE),
              op_arg_dat(exteriorQ[2], -1, OP_ID, 3 * 5, "double", OP_WRITE),
              op_arg_dat(exteriorQ[3], -1, OP_ID, 3 * 5, "double", OP_WRITE));

  double dt1 = 0.0;
  op_par_loop(calc_dt, "calc_dt", cells,
              op_arg_dat(Q[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[2], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[3], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(fscale, -1, OP_ID, 3 * 5, "double", OP_READ),
              op_arg_gbl(&dt1, 1, "double", OP_MAX));

  dt = 1.0 / dt1;
  cout << "dt: " << dt << endl;
  double t = 0.0;

  double rk_t_frac[3] = {0.0, 1.0, 0.5};
  // Timers
  double cpu_loop_start, wall_loop_start, cpu_loop_end, wall_loop_end, cpu_loop_1, wall_loop_1, cpu_loop_2, wall_loop_2;
  double get_neighbour_q_t = 0.0;
  double get_bedge_q_t = 0.0;
  double euler_rhs_t = 0.0;
  double set_workingQ_t = 0.0;
  double update_Q_t = 0.0;
  double calc_dt_t = 0.0;
  double internal_fluxes_t = 0.0;
  double internal_fluxes_mat_t = 0.0;
  double face_fluxes_mat_t = 0.0;
  op_timers(&cpu_loop_start, &wall_loop_start);

  Vec q;
  VecCreateSeqCUDA(PETSC_COMM_SELF, 4 * 15 * numCells, &q);
  double *q_d;
  VecCUDAGetArray(q, &q_d);
  op_arg q_vec_petsc_args[] = {
    op_arg_dat(Q[0], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(Q[1], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(Q[2], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(Q[3], -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(cells, 4, q_vec_petsc_args);
  cudaMemcpy(q_d, (double *)Q[0]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(q_d + 15 *numCells, (double *)Q[1]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(q_d + 2 * 15 *numCells, (double *)Q[2]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(q_d + 3 * 15 *numCells, (double *)Q[3]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(4, q_vec_petsc_args);
  VecCUDARestoreArray(q, &q_d);
  Vec q1;
  VecCreateSeqCUDA(PETSC_COMM_SELF, 4 * 15 * numCells, &q1);
  double *q1_d;
  VecCUDAGetArray(q1, &q1_d);
  op_mpi_halo_exchanges_cuda(cells, 4, q_vec_petsc_args);
  cudaMemcpy(q1_d, (double *)Q[0]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(q1_d + 15 *numCells, (double *)Q[1]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(q1_d + 2 * 15 *numCells, (double *)Q[2]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(q1_d + 3 * 15 *numCells, (double *)Q[3]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(4, q_vec_petsc_args);
  VecCUDARestoreArray(q1, &q1_d);
  // Create A matrix
  Mat Amat;
  MatCreateShell(PETSC_COMM_SELF, 4 * 15 * numCells, 4 * 15 * numCells, PETSC_DETERMINE, PETSC_DETERMINE, NULL, &Amat);
  MatShellSetOperation(Amat, MATOP_MULT, (void(*)(void))matAMult);
  // Create KSP
  KSP ksp;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetOperators(ksp, Amat, Amat);
  // Backwards Euler
  for(int i = 0; i < iter; i++) {
    // Solve
    KSPSolve(ksp, q, q1);

    // Update solution
    VecCUDAGetArray(q, &q_d);
    VecCUDAGetArray(q1, &q1_d);

    cudaMemcpy(q_d, q1_d, 4 * 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);

    VecCUDARestoreArray(q, &q_d);
    VecCUDARestoreArray(q1, &q1_d);
  }
  VecCUDAGetArray(q, &q_d);
  op_arg q_write_vec_petsc_args[] = {
    op_arg_dat(Q[0], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(Q[1], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(Q[2], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(Q[3], -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(cells, 4, q_write_vec_petsc_args);
  cudaMemcpy((double *)Q[0]->data_d, q_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy((double *)Q[1]->data_d, q_d + 15 *numCells, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy((double *)Q[2]->data_d, q_d + 2 * 15 *numCells, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy((double *)Q[3]->data_d, q_d + 3 * 15 *numCells, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(4, q_write_vec_petsc_args);
  VecCUDARestoreArray(q, &q_d);

  // VecDestory(q);
  // VecDestory(q1);
  // KSPDestory(ksp);

  /*
  // SSP RK-3
  for(int i = 0; i < iter; i++) {
    for(int j = 0; j < 3; j++) {
      get_RHS(workingQ, rk[j]);

      if(j != 2) {
        // op_timers(&cpu_loop_1, &wall_loop_1);
        op_par_loop(set_workingQ, "set_workingQ", cells,
                    op_arg_gbl(&dt, 1, "double", OP_READ),
                    op_arg_gbl(&j, 1, "int", OP_READ),
                    op_arg_dat(Q[0], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(Q[1], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(Q[2], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(Q[3], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(rk[0][0], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(rk[0][1], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(rk[0][2], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(rk[0][3], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(rk[1][0], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(rk[1][1], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(rk[1][2], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(rk[1][3], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(workingQ[0], -1, OP_ID, 15, "double", OP_WRITE),
                    op_arg_dat(workingQ[1], -1, OP_ID, 15, "double", OP_WRITE),
                    op_arg_dat(workingQ[2], -1, OP_ID, 15, "double", OP_WRITE),
                    op_arg_dat(workingQ[3], -1, OP_ID, 15, "double", OP_WRITE));
        // op_timers(&cpu_loop_2, &wall_loop_2);
        // set_workingQ_t += wall_loop_2 - wall_loop_1;
      }
    }

    op_timers(&cpu_loop_1, &wall_loop_1);
    op_par_loop(update_Q, "update_Q", cells,
                op_arg_gbl(&dt, 1, "double", OP_READ),
                op_arg_dat(Q[0], -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(Q[1], -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(Q[2], -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(Q[3], -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(rk[0][0], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[0][1], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[0][2], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[0][3], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[1][0], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[1][1], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[1][2], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[1][3], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[2][0], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[2][1], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[2][2], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[2][3], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(workingQ[0], -1, OP_ID, 15, "double", OP_WRITE),
                op_arg_dat(workingQ[1], -1, OP_ID, 15, "double", OP_WRITE),
                op_arg_dat(workingQ[2], -1, OP_ID, 15, "double", OP_WRITE),
                op_arg_dat(workingQ[3], -1, OP_ID, 15, "double", OP_WRITE));
    op_timers(&cpu_loop_2, &wall_loop_2);
    update_Q_t += wall_loop_2 - wall_loop_1;

    t += dt;
    dt1 = 0.0;
    op_timers(&cpu_loop_1, &wall_loop_1);
    op_par_loop(calc_dt, "calc_dt", cells,
                op_arg_dat(Q[0], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(Q[1], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(Q[2], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(Q[3], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(fscale, -1, OP_ID, 3 * 5, "double", OP_READ),
                op_arg_gbl(&dt1, 1, "double", OP_MAX));
    op_timers(&cpu_loop_2, &wall_loop_2);
    calc_dt_t += wall_loop_2 - wall_loop_1;

    dt = 1.0 / dt1;
    if(i % 1000 == 0)
      cout << "iter: " << i << " time: " << t <<  " dt: " << dt << endl;
  }
  op_timers(&cpu_loop_end, &wall_loop_end);

  cout << "Time: " << t << endl;
  */

  // Save info for python test script
  // op_fetch_data_hdf5_file(node_coords, "points.h5");
  // op_fetch_data_hdf5_file(x, "points.h5");
  // op_fetch_data_hdf5_file(y, "points.h5");
  // op_fetch_data_hdf5_file(nx, "points.h5");
  // op_fetch_data_hdf5_file(ny, "points.h5");

  // Save solution to CGNS file
  double *sol_q0 = (double *)malloc(15 * op_get_size(cells) * sizeof(double));
  double *sol_q1 = (double *)malloc(15 * op_get_size(cells) * sizeof(double));
  double *sol_q2 = (double *)malloc(15 * op_get_size(cells) * sizeof(double));
  double *sol_q3 = (double *)malloc(15 * op_get_size(cells) * sizeof(double));
  op_fetch_data(Q[0], sol_q0);
  op_fetch_data(Q[1], sol_q1);
  op_fetch_data(Q[2], sol_q2);
  op_fetch_data(Q[3], sol_q3);
  save_solution("./naca0012.cgns", op_get_size(nodes), op_get_size(cells),
                sol_q0, sol_q1, sol_q2, sol_q3, cgnsCells, gam);

  free(sol_q0);
  free(sol_q1);
  free(sol_q2);
  free(sol_q3);

  op_timers(&cpu_2, &wall_2);

  op_timing_output();

  cout << endl << "Total execution time: " << wall_2 - wall_1 << endl;
  cout << "Total time in main loop: " << wall_loop_end - wall_loop_start << endl;
  cout << "Average time per iteration: " << (wall_loop_end - wall_loop_start) / iter << endl;
  cout << "*** Kernels ***" << endl;
  cout << "get_neighbour_q" << endl;
  cout << "  Total: " << get_neighbour_q_t << endl;
  // cout << "  Per iter: " << get_neighbour_q_t / iter << endl;
  cout << "get_bedge_q" << endl;
  cout << "  Total: " << get_bedge_q_t << endl;
  // cout << "  Per iter: " << get_bedge_q_t / iter << endl;
  cout << "euler_rhs" << endl;
  cout << "  Total: " << euler_rhs_t << endl;
  // cout << "  Per iter: " << euler_rhs_t / iter << endl;
  cout << "set_workingQ" << endl;
  cout << "  Total: " << set_workingQ_t << endl;
  // cout << "  Per iter: " << set_workingQ_t / iter << endl;
  cout << "update_Q" << endl;
  cout << "  Total: " << update_Q_t << endl;
  // cout << "  Per iter: " << update_Q_t / iter << endl;
  cout << "calc_dt" << endl;
  cout << "  Total: " << calc_dt_t << endl;
  // cout << "  Per iter: " << calc_dt_t / iter << endl;
  cout << "internal_fluxes" << endl;
  cout << "  Total: " << internal_fluxes_t << endl;
  // cout << "  Per iter: " << set_workingQ_t / iter << endl;
  cout << "internal_fluxes_mat" << endl;
  cout << "  Total: " << internal_fluxes_mat_t << endl;
  // cout << "  Per iter: " << update_Q_t / iter << endl;
  cout << "face_fluxes_mat" << endl;
  cout << "  Total: " << face_fluxes_mat_t << endl;
  // cout << "  Per iter: " << calc_dt_t / iter << endl;

  cout << endl << "Estimate wall time to simulate 1 second: " << (wall_loop_end - wall_loop_start) / t << endl;


  // Clean up OP2
  op_exit();

  free(coords);
  free(cgnsCells);
  free(edge2node_data);
  free(edge2cell_data);
  free(bedge2node_data);
  free(bedge2cell_data);
  free(edgeNum_data);
  free(bedgeNum_data);
  free(bedge_type_data);

  delete data;

  cublasDestroy(cublas_handle);

  ierr = PetscFinalize();
  return ierr;
}

inline void get_RHS(op_dat *Q, op_dat *rhs) {
  int numCells = op_get_size(cells);
  // Get neighbouring values of q on internal edges
  // op_timers(&cpu_loop_1, &wall_loop_1);
  op_par_loop(get_neighbour_q, "get_neighbour_q", edges,
              op_arg_dat(edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(nodeX, 0, edge2cells, 3, "double", OP_READ),
              op_arg_dat(nodeY, 0, edge2cells, 3, "double", OP_READ),
              op_arg_dat(nodeX, 1, edge2cells, 3, "double", OP_READ),
              op_arg_dat(nodeY, 1, edge2cells, 3, "double", OP_READ),
              op_arg_dat(Q[0], 0, edge2cells, 15, "double", OP_READ),
              op_arg_dat(Q[1], 0, edge2cells, 15, "double", OP_READ),
              op_arg_dat(Q[2], 0, edge2cells, 15, "double", OP_READ),
              op_arg_dat(Q[3], 0, edge2cells, 15, "double", OP_READ),
              op_arg_dat(Q[0], 1, edge2cells, 15, "double", OP_READ),
              op_arg_dat(Q[1], 1, edge2cells, 15, "double", OP_READ),
              op_arg_dat(Q[2], 1, edge2cells, 15, "double", OP_READ),
              op_arg_dat(Q[3], 1, edge2cells, 15, "double", OP_READ),
              op_arg_dat(exteriorQ[0], 0, edge2cells, 3 * 5, "double", OP_INC),
              op_arg_dat(exteriorQ[1], 0, edge2cells, 3 * 5, "double", OP_INC),
              op_arg_dat(exteriorQ[2], 0, edge2cells, 3 * 5, "double", OP_INC),
              op_arg_dat(exteriorQ[3], 0, edge2cells, 3 * 5, "double", OP_INC),
              op_arg_dat(exteriorQ[0], 1, edge2cells, 3 * 5, "double", OP_INC),
              op_arg_dat(exteriorQ[1], 1, edge2cells, 3 * 5, "double", OP_INC),
              op_arg_dat(exteriorQ[2], 1, edge2cells, 3 * 5, "double", OP_INC),
              op_arg_dat(exteriorQ[3], 1, edge2cells, 3 * 5, "double", OP_INC));
  // op_timers(&cpu_loop_2, &wall_loop_2);
  // get_neighbour_q_t += wall_loop_2 - wall_loop_1;

  // Enforce boundary conditions
  // op_timers(&cpu_loop_1, &wall_loop_1);
  op_par_loop(get_bedge_q, "get_bedge_q", bedges,
              op_arg_dat(bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(nx, 0, bedge2cells, 3 * 5, "double", OP_READ),
              op_arg_dat(ny, 0, bedge2cells, 3 * 5, "double", OP_READ),
              op_arg_dat(Q[0], 0, bedge2cells, 15, "double", OP_READ),
              op_arg_dat(Q[1], 0, bedge2cells, 15, "double", OP_READ),
              op_arg_dat(Q[2], 0, bedge2cells, 15, "double", OP_READ),
              op_arg_dat(Q[3], 0, bedge2cells, 15, "double", OP_READ),
              op_arg_dat(exteriorQ[0], 0, bedge2cells, 3 * 5, "double", OP_INC),
              op_arg_dat(exteriorQ[1], 0, bedge2cells, 3 * 5, "double", OP_INC),
              op_arg_dat(exteriorQ[2], 0, bedge2cells, 3 * 5, "double", OP_INC),
              op_arg_dat(exteriorQ[3], 0, bedge2cells, 3 * 5, "double", OP_INC));
  // op_timers(&cpu_loop_2, &wall_loop_2);
  // get_bedge_q_t += wall_loop_2 - wall_loop_1;

  // op_timers(&cpu_loop_1, &wall_loop_1);
  op_par_loop(internal_fluxes, "internal_fluxes", cells,
              op_arg_dat(Q[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[2], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[3], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(F[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(F[1], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(F[2], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(F[3], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(G[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(G[1], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(G[2], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(G[3], -1, OP_ID, 15, "double", OP_WRITE));
  // op_timers(&cpu_loop_2, &wall_loop_2);
  // internal_fluxes_t += wall_loop_2 - wall_loop_1;

  // TODO matrix mult
  // op_timers(&cpu_loop_1, &wall_loop_1);
  op_arg internal_fluxes_args[] = {
    op_arg_dat(F[0], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(F[1], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(F[2], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(F[3], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(G[0], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(G[1], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(G[2], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(G[3], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(dFdr[0], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dFdr[1], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dFdr[2], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dFdr[3], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dFds[0], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dFds[1], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dFds[2], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dFds[3], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dGdr[0], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dGdr[1], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dGdr[2], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dGdr[3], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dGds[0], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dGds[1], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dGds[2], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(dGds[3], -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(cells, 24, internal_fluxes_args);
  internal_fluxes_matrices(cublas_handle, numCells, (double *)F[0]->data_d,
                           (double *)F[1]->data_d, (double *)F[2]->data_d,
                           (double *)F[3]->data_d, (double *)G[0]->data_d,
                           (double *)G[1]->data_d, (double *)G[2]->data_d,
                           (double *)G[3]->data_d, (double *)dFdr[0]->data_d,
                           (double *)dFdr[1]->data_d, (double *)dFdr[2]->data_d,
                           (double *)dFdr[3]->data_d, (double *)dFds[0]->data_d,
                           (double *)dFds[1]->data_d, (double *)dFds[2]->data_d,
                           (double *)dFds[3]->data_d, (double *)dGdr[0]->data_d,
                           (double *)dGdr[1]->data_d, (double *)dGdr[2]->data_d,
                           (double *)dGdr[3]->data_d, (double *)dGds[0]->data_d,
                           (double *)dGds[1]->data_d, (double *)dGds[2]->data_d,
                           (double *)dGds[3]->data_d);

  // Check this
  op_mpi_set_dirtybit_cuda(24, internal_fluxes_args);
  // dFdr->dirty_hd = 2;
  // dFds->dirty_hd = 2;
  // dGdr->dirty_hd = 2;
  // dGds->dirty_hd = 2;
  // op_timers(&cpu_loop_2, &wall_loop_2);
  // internal_fluxes_mat_t += wall_loop_2 - wall_loop_1;

  // Calculate vectors F an G from q for each cell
  // op_timers(&cpu_loop_1, &wall_loop_1);

  op_par_loop(euler_rhs, "euler_rhs", cells,
              op_arg_dat(Q[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[2], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[3], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(exteriorQ[0], -1, OP_ID, 3 * 5, "double", OP_RW),
              op_arg_dat(exteriorQ[1], -1, OP_ID, 3 * 5, "double", OP_RW),
              op_arg_dat(exteriorQ[2], -1, OP_ID, 3 * 5, "double", OP_RW),
              op_arg_dat(exteriorQ[3], -1, OP_ID, 3 * 5, "double", OP_RW),
              op_arg_dat(rx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(ry, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(sx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(sy, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(fscale, -1, OP_ID, 3 * 5, "double", OP_READ),
              op_arg_dat(nx, -1, OP_ID, 3 * 5, "double", OP_READ),
              op_arg_dat(ny, -1, OP_ID, 3 * 5, "double", OP_READ),
              op_arg_dat(dFdr[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dFdr[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dFdr[2], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dFdr[3], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dFds[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dFds[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dFds[2], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dFds[3], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dGdr[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dGdr[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dGdr[2], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dGdr[3], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dGds[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dGds[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dGds[2], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dGds[3], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(flux[0], -1, OP_ID, 3 * 5, "double", OP_WRITE),
              op_arg_dat(flux[1], -1, OP_ID, 3 * 5, "double", OP_WRITE),
              op_arg_dat(flux[2], -1, OP_ID, 3 * 5, "double", OP_WRITE),
              op_arg_dat(flux[3], -1, OP_ID, 3 * 5, "double", OP_WRITE),
              op_arg_dat(rhs[0], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(rhs[1], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(rhs[2], -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(rhs[3], -1, OP_ID, 15, "double", OP_WRITE));
  // op_timers(&cpu_loop_2, &wall_loop_2);
  // euler_rhs_t += wall_loop_2 - wall_loop_1;

  // TODO matrix mult
  // op_timers(&cpu_loop_1, &wall_loop_1);
  op_arg face_fluxes_args[] = {
    op_arg_dat(flux[0], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(flux[1], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(flux[2], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(flux[3], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(rhs[0], -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(rhs[1], -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(rhs[2], -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(rhs[3], -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges_cuda(cells, 8, face_fluxes_args);
  face_fluxes_matrices(cublas_handle, numCells, (double *)flux[0]->data_d,
                       (double *)flux[1]->data_d, (double *)flux[2]->data_d,
                       (double *)flux[3]->data_d, (double *)rhs[0]->data_d,
                       (double *)rhs[1]->data_d, (double *)rhs[2]->data_d,
                       (double *)rhs[3]->data_d);

  // Check this
  op_mpi_set_dirtybit_cuda(8, face_fluxes_args);
  // rk[j]->dirty_hd = 2;
  // op_timers(&cpu_loop_2, &wall_loop_2);
  // face_fluxes_mat_t += wall_loop_2 - wall_loop_1;
}

PetscErrorCode matAMult(Mat A, Vec x, Vec y) {
  int numCells = op_get_size(cells);
  // Get PETSC data on GPU
  double *q;
  double *q1;
  VecCUDAGetArray(x, &q);
  VecCUDAGetArray(y, &q1);

  // Construct OP2 temp dats

  // Set data for OP2 dats
  op_arg set_temp_args[] = {
    op_arg_dat(Q[0], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(Q[1], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(Q[2], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(Q[3], -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(cells, 4, set_temp_args);

  cudaMemcpy((double *)Q[0]->data_d, q, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy((double *)Q[1]->data_d, q + 15 * numCells, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy((double *)Q[2]->data_d, q + 2 * 15 * numCells, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy((double *)Q[3]->data_d, q + 3 * 15 * numCells, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);

  op_mpi_set_dirtybit_cuda(4, set_temp_args);

  // Calc Euler RHS
  get_RHS(Q, rhs);

  // Set y vec
  double dt = 1e-6;
  op_par_loop(backwards_euler_update_Q, "backwards_euler_update_Q", cells,
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(Q[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[2], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(Q[3], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(rhs[0], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(rhs[1], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(rhs[2], -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(rhs[3], -1, OP_ID, 15, "double", OP_RW));

  op_arg set_y_args[] = {
    op_arg_dat(rhs[0], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(rhs[1], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(rhs[2], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(rhs[3], -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(cells, 4, set_y_args);

  cudaMemcpy(q1, (double *)rhs[0]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(q1 + 15 * numCells, (double *)rhs[1]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(q1 + 2 * 15 * numCells, (double *)rhs[2]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(q1 + 3 * 15 * numCells, (double *)rhs[3]->data_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);

  op_mpi_set_dirtybit_cuda(4, set_y_args);

  // Release PETSC data
  VecCUDARestoreArray(x, &q);
  VecCUDARestoreArray(y, &q1);

  return 0;
}
