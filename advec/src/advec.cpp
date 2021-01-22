#define NUM_SOLUTION_PTS 15
#define NUM_FACE_PTS 5
#define ORDER 4
#define _USE_MATH_DEFINES

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

// Include kernels
#include "advec_rhs.h"
#include "get_neighbour_q.h"
#include "get_bedge_q.h"
#include "set_ic.h"
#include "set_workingQ.h"
#include "update_Q.h"
#include "init_grid.h"
#include "neighbour_zero.h"

using namespace std;

// Stuff for parsing command line arguments
extern char *optarg;
extern int  optind, opterr, optopt;
static struct option options[] = {
  {"iter", required_argument, 0, 0},
  {"vx", required_argument,  0,  0},
  {"vy", required_argument,  0,  0},
  {0,    0,                  0,  0}
};

int main(int argc, char **argv) {
  // Load mesh
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

  load_mesh("./grid.cgns", &numNodes, &numCells, &numEdges,
            &numBoundaryEdges, &coords, &cgnsCells, &edge2node_data,
            &edge2cell_data, &edgeNum_data, &bedge2node_data, &bedge2cell_data,
            &bedgeNum_data);

  // Initialise OP2
  op_init(argc, argv, 2);

  // Get input from args
  int iter = 100;

  int opt_index = 0;
  while(getopt_long_only(argc, argv, "", options, &opt_index) != -1) {
    if(strcmp((char*)options[opt_index].name,"vx"  ) == 0) velX = stod(optarg);
    if(strcmp((char*)options[opt_index].name,"vy"  ) == 0) velY = stod(optarg);
    if(strcmp((char*)options[opt_index].name,"iter") == 0) iter = atoi(optarg);
  }

  // Declare memory for data that will be calculated in initialisation kernel
  double *nodeX_data = (double*)malloc(3 * numCells * sizeof(double));
  double *nodeY_data = (double*)malloc(3 * numCells * sizeof(double));
  double *x_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *y_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *xr_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *yr_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *xs_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *ys_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *rx_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *ry_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *sx_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *sy_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *nx_data = (double *)malloc(3 * NUM_FACE_PTS * numCells * sizeof(double));
  double *ny_data = (double *)malloc(3 * NUM_FACE_PTS * numCells * sizeof(double));
  double *fscale_data = (double *)malloc(3 * NUM_FACE_PTS * numCells * sizeof(double));
  double *q_data  = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *workingQ_data  = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *exteriorQ_data = (double *)malloc(3 * NUM_FACE_PTS * numCells * sizeof(double));
  // RK4 data
  double *rk1_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *rk2_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *rk3_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *rk4_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));

  // Declare OP2 sets
  op_set nodes = op_decl_set(numNodes, "nodes");
  op_set cells = op_decl_set(numCells, "cells");
  op_set edges = op_decl_set(numEdges, "edges");
  op_set bedges = op_decl_set(numBoundaryEdges, "bedges");

  // Declare OP2 maps
  op_map cell2nodes = op_decl_map(cells, nodes, 3, cgnsCells, "cell2nodes");
  op_map edge2nodes = op_decl_map(edges, nodes, 2, edge2node_data, "edge2nodes");
  op_map edge2cells = op_decl_map(edges, cells, 2, edge2cell_data, "edge2cells");
  op_map bedge2nodes = op_decl_map(bedges, nodes, 2, bedge2node_data, "bedge2nodes");
  op_map bedge2cells = op_decl_map(bedges, cells, 1, bedge2cell_data, "bedge2cells");

  // Declare OP2 datasets
    // Structure: {x, y}
  op_dat node_coords = op_decl_dat(nodes, 2, "double", coords, "node_coords");
    // Coords of nodes per cell
  op_dat nodeX = op_decl_dat(cells, 3, "double", nodeX_data, "nodeX");
  op_dat nodeY = op_decl_dat(cells, 3, "double", nodeY_data, "nodeY");
    // The x and y coordinates of all the solution points in a cell
  op_dat x = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", x_data, "x");
  op_dat y = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", y_data, "y");
    // Geometric factors that relate to mapping between global and local (cell) coordinates
  op_dat xr = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", xr_data, "xr");
  op_dat yr = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", yr_data, "yr");
  op_dat xs = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", xs_data, "xs");
  op_dat ys = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", ys_data, "ys");
  op_dat rx = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", rx_data, "rx");
  op_dat ry = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", ry_data, "ry");
  op_dat sx = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", sx_data, "sx");
  op_dat sy = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", sy_data, "sy");
    // Normals for each cell (calculated for each node on each edge, nodes can appear on multiple edges)
  op_dat nx = op_decl_dat(cells, 3 * NUM_FACE_PTS, "double", nx_data, "nx");
  op_dat ny = op_decl_dat(cells, 3 * NUM_FACE_PTS, "double", ny_data, "ny");
    // surface Jacobian / Jacobian (used when lifting the boundary fluxes)
  op_dat fscale = op_decl_dat(cells, 3 * NUM_FACE_PTS, "double", fscale_data, "fscale");
    // Q data sets
  op_dat q        = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", q_data, "q");
  op_dat workingQ = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", workingQ_data, "workingQ");
  op_dat rk[4];
  rk[0] = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", rk1_data, "rk1");
  rk[1] = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", rk2_data, "rk2");
  rk[2] = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", rk3_data, "rk3");
  rk[3] = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", rk4_data, "rk4");
    // Holds neighbouring values of Q for nodes on faces
  op_dat exteriorQ = op_decl_dat(cells, 3 * NUM_FACE_PTS, "double", exteriorQ_data, "exteriorQ");
  op_dat edgeNum   = op_decl_dat(edges, 2, "int", edgeNum_data, "edgeNum");
  op_dat bedgeNum  = op_decl_dat(bedges, 1, "int", bedgeNum_data, "bedgeNum");

  // Declare OP2 constants
  op_decl_const(1, "double", &velX);
  op_decl_const(1, "double", &velY);
  op_decl_const(NUM_SOLUTION_PTS, "double", ones);
  op_decl_const(NUM_SOLUTION_PTS, "double", r);
  op_decl_const(NUM_SOLUTION_PTS, "double", s);
  op_decl_const(NUM_SOLUTION_PTS * NUM_SOLUTION_PTS, "double", Dr);
  op_decl_const(NUM_SOLUTION_PTS * NUM_SOLUTION_PTS, "double", Ds);
  op_decl_const(NUM_SOLUTION_PTS * NUM_SOLUTION_PTS, "double", Drw);
  op_decl_const(NUM_SOLUTION_PTS * NUM_SOLUTION_PTS, "double", Dsw);
  op_decl_const(3 * NUM_FACE_PTS, "int", FMASK);
  op_decl_const(NUM_SOLUTION_PTS * NUM_SOLUTION_PTS, "double", LIFT);

  double t = 0.0;

  // Initialisation kernels
  op_par_loop(init_grid, "init_grid", cells,
              op_arg_dat(node_coords, 0, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(node_coords, 1, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(node_coords, 2, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(nodeX, -1, OP_ID, 3, "double", OP_WRITE),
              op_arg_dat(nodeY, -1, OP_ID, 3, "double", OP_WRITE),
              op_arg_dat(x, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(y, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(xr, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(yr, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(xs, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(ys, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(rx, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(ry, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(sx, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(sy, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(nx, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_WRITE),
              op_arg_dat(ny, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_WRITE),
              op_arg_dat(fscale, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_WRITE));

  op_par_loop(set_ic, "set_ic", cells,
              op_arg_dat(x, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
              op_arg_dat(y, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
              op_arg_dat(q, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(workingQ, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE));

  dt = 0.001;
  cout << "dt: " << dt << endl;

  double rk_frac[3] = {0.5, 0.5, 1.0};
  // Run the simulation
  for(int i = 0; i < iter; i++) {
    for(int j = 0; j < 4; j++) {
      op_par_loop(neighbour_zero, "neighbour_zero", cells,
                  op_arg_dat(exteriorQ, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_WRITE));
      // Get neighbouring values of q on internal edges
      op_par_loop(get_neighbour_q, "get_neighbour_q", edges,
                  op_arg_dat(edgeNum, -1, OP_ID, 2, "int", OP_READ),
                  op_arg_dat(node_coords, 0, edge2nodes, 2, "double", OP_READ),
                  op_arg_dat(node_coords, 1, edge2nodes, 2, "double", OP_READ),
                  op_arg_dat(nodeX, 0, edge2cells, 3, "double", OP_READ),
                  op_arg_dat(nodeY, 0, edge2cells, 3, "double", OP_READ),
                  op_arg_dat(nodeX, 1, edge2cells, 3, "double", OP_READ),
                  op_arg_dat(nodeY, 1, edge2cells, 3, "double", OP_READ),
                  op_arg_dat(workingQ, 0, edge2cells, NUM_SOLUTION_PTS, "double", OP_READ),
                  op_arg_dat(workingQ, 1, edge2cells, NUM_SOLUTION_PTS, "double", OP_READ),
                  op_arg_dat(exteriorQ, 0, edge2cells, 3 * NUM_FACE_PTS, "double", OP_INC),
                  op_arg_dat(exteriorQ, 1, edge2cells, 3 * NUM_FACE_PTS, "double", OP_INC));

      // Enforce boundary conditions
      op_par_loop(get_bedge_q, "get_bedge_q", bedges,
                  op_arg_dat(bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                  op_arg_dat(node_coords, 0, bedge2nodes, 2, "double", OP_READ),
                  op_arg_dat(node_coords, 1, bedge2nodes, 2, "double", OP_READ),
                  op_arg_dat(nodeX, 0, bedge2cells, 3, "double", OP_READ),
                  op_arg_dat(nodeY, 0, bedge2cells, 3, "double", OP_READ),
                  op_arg_dat(x, 0, bedge2cells, NUM_SOLUTION_PTS, "double", OP_READ),
                  op_arg_dat(y, 0, bedge2cells, NUM_SOLUTION_PTS, "double", OP_READ),
                  op_arg_dat(exteriorQ, 0, bedge2cells, 3 * NUM_FACE_PTS, "double", OP_INC));

      // Calculate vectors F an G from q for each cell
      op_par_loop(advec_rhs, "advec_rhs", cells,
                  op_arg_dat(workingQ, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                  op_arg_dat(exteriorQ, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_READ),
                  op_arg_dat(rx, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                  op_arg_dat(ry, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                  op_arg_dat(sx, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                  op_arg_dat(sy, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                  op_arg_dat(fscale, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_READ),
                  op_arg_dat(nx, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_READ),
                  op_arg_dat(ny, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_READ),
                  op_arg_dat(rk[j], -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE));

      if(j != 3) {
        op_par_loop(set_workingQ, "set_workingQ", cells,
                    op_arg_gbl(&dt, 1, "double", OP_READ),
                    op_arg_dat(q, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                    op_arg_dat(rk[j], -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                    op_arg_gbl(&rk_frac[j], 1, "double", OP_READ),
                    op_arg_dat(workingQ, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE));
      }
    }
    op_par_loop(update_Q, "update_Q", cells,
                op_arg_gbl(&dt, 1, "double", OP_READ),
                op_arg_dat(q, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_RW),
                op_arg_dat(rk[0], -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(rk[1], -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(rk[2], -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(rk[3], -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(workingQ, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE));

    t += dt;
  }

  cout << "t: " << t << endl;

  // Save solution to CGNS file
  double *sol_q = (double *)malloc(NUM_SOLUTION_PTS * op_get_size(cells) * sizeof(double));
  op_fetch_data(q, sol_q);
  save_solution("./grid.cgns", op_get_size(nodes), op_get_size(cells), sol_q, cgnsCells, gam);

  free(sol_q);

  // Clean up OP2
  op_exit();

  free(coords);
  free(cgnsCells);
  free(nodeX_data);
  free(nodeY_data);
  free(x_data);
  free(y_data);
  free(xr_data);
  free(yr_data);
  free(xs_data);
  free(ys_data);
  free(rx_data);
  free(ry_data);
  free(sx_data);
  free(sy_data);
  free(nx_data);
  free(ny_data);
  free(fscale_data);
  free(q_data);
  free(workingQ_data);
  free(edge2node_data);
  free(edge2cell_data);
  free(bedge2node_data);
  free(bedge2cell_data);
  free(exteriorQ_data);
  free(rk1_data);
  free(rk2_data);
  free(rk3_data);
  free(rk4_data);
  free(edgeNum_data);
  free(bedgeNum_data);
}
