#define NUM_SOLUTION_PTS 15
#define NUM_FACE_PTS 5

// Include OP2 stuff
#include "op_seq.h"
// Include C++ stuff
#include <string>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>

#include "constantsDG.h"
#include "load_mesh.h"

// Include kernels
#include "calc_solution_pt_coords.h"
#include "calc_geometric_factors.h"
#include "calc_normals.h"
#include "euler_rhs.h"
#include "get_neighbour_q.h"

using namespace std;

int main(int argc, char **argv) {
  double *coords;
  int *cgnsCells;
  int *edge2node_data;
  int *edge2cell_data;
  int *bedge2node_data;
  int *bedge2cell_data;
  int numNodes, numCells, numEdges, numBoundaryEdges;

  load_mesh("./naca0012.cgns", &numNodes, &numCells, &numEdges,
            &numBoundaryEdges, &coords, &cgnsCells, &edge2node_data,
            &edge2cell_data, &bedge2node_data, &bedge2cell_data);

  // Initialise OP2
  op_init(argc, argv, 2);

  // TODO get input from args
  int iter = 1;

  // Declare memory for data that will be calculated in initialisation kernel
  double *nodeX_data = (double*)malloc(3 * numCells * sizeof(double));
  double *nodeY_data = (double*)malloc(3 * numCells * sizeof(double));
  double *x_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *y_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *xr_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *yr_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *xs_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *ys_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *J_data  = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *rx_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *ry_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *sx_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *sy_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *nx_data = (double *)malloc(3 * NUM_FACE_PTS * numCells * sizeof(double));
  double *ny_data = (double *)malloc(3 * NUM_FACE_PTS * numCells * sizeof(double));
  double *sJ_data = (double *)malloc(3 * NUM_FACE_PTS * numCells * sizeof(double));
  double *q_data  = (double *)malloc(4 * NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *qRHS_data  = (double *)malloc(4 * NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *exteriorQ_data = (double *)malloc(4 * 3 * NUM_FACE_PTS * numCells * sizeof(double));

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
    // Jacobian
  op_dat J  = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", J_data, "J");
  op_dat rx = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", rx_data, "rx");
  op_dat ry = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", ry_data, "ry");
  op_dat sx = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", sx_data, "sx");
  op_dat sy = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", sy_data, "sy");
    // Normals for each cell (calculated for each node on each edge, nodes can appear on multiple edges)
  op_dat nx = op_decl_dat(cells, 3 * NUM_FACE_PTS, "double", nx_data, "nx");
  op_dat ny = op_decl_dat(cells, 3 * NUM_FACE_PTS, "double", ny_data, "ny");
    // Surface Jacobian
  op_dat sJ = op_decl_dat(cells, 3 * NUM_FACE_PTS, "double", sJ_data, "sJ");
    // Values for compressible Euler equations in vectors
    // Structure: {q0_0, q1_0, q2_0, q3_0, q0_1, q1_1, ..., q3_NUM_SOLUTION_PTS}
  op_dat q    = op_decl_dat(cells, 4 * NUM_SOLUTION_PTS, "double", q_data, "q");
  op_dat qRHS = op_decl_dat(cells, 4 * NUM_SOLUTION_PTS, "double", qRHS_data, "qRHS");
    // Holds neighbouring values of Q for nodes on faces
  op_dat exteriorQ = op_decl_dat(cells, 4 * 3 * NUM_FACE_PTS, "double", exteriorQ_data, "exteriorQ");

  // Declare OP2 constants
  op_decl_const(NUM_SOLUTION_PTS, "double", ones);
  op_decl_const(NUM_SOLUTION_PTS, "double", solution_pts_r);
  op_decl_const(NUM_SOLUTION_PTS, "double", solution_pts_s);
  op_decl_const(NUM_SOLUTION_PTS * NUM_SOLUTION_PTS, "double", Dr);
  op_decl_const(NUM_SOLUTION_PTS * NUM_SOLUTION_PTS, "double", Ds);
  op_decl_const(1, "double", &gam);
  op_decl_const(NUM_SOLUTION_PTS * NUM_SOLUTION_PTS, "double", Drw);
  op_decl_const(NUM_SOLUTION_PTS * NUM_SOLUTION_PTS, "double", Dsw);
  op_decl_const(NUM_FACE_PTS, "int", fmask0);
  op_decl_const(NUM_FACE_PTS, "int", fmask1);
  op_decl_const(NUM_FACE_PTS, "int", fmask2);
  op_decl_const(NUM_SOLUTION_PTS * NUM_SOLUTION_PTS, "double", lift);

  // Initialisation kernels
  op_par_loop(calc_solution_pt_coords, "calc_solution_pt_coords", cells,
              op_arg_dat(node_coords, 0, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(node_coords, 1, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(node_coords, 2, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(nodeX, -1, OP_ID, 3, "double", OP_WRITE),
              op_arg_dat(nodeY, -1, OP_ID, 3, "double", OP_WRITE),
              op_arg_dat(x, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(y, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE));

  op_par_loop(calc_geometric_factors, "calc_geometric_factors", cells,
              op_arg_dat(x, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
              op_arg_dat(y, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
              op_arg_dat(xr, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(yr, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(xs, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(ys, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(J, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(rx, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(ry, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(sx, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(sy, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE));

  op_par_loop(calc_normals, "calc_normals", cells,
              op_arg_dat(xr, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
              op_arg_dat(yr, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
              op_arg_dat(xs, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
              op_arg_dat(ys, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
              op_arg_dat(nx, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_WRITE),
              op_arg_dat(ny, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_WRITE),
              op_arg_dat(sJ, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_WRITE));

  // Run the simulation
  for(int i = 0; i < iter; i++) {
    // Get neighbouring values of q on internal edges
    op_par_loop(get_neighbour_q, "get_neighbour_q", edges,
                op_arg_dat(node_coords, 0, edge2nodes, 2, "double", OP_READ),
                op_arg_dat(node_coords, 1, edge2nodes, 2, "double", OP_READ),
                op_arg_dat(nodeX, 0, edge2cells, 3, "double", OP_READ),
                op_arg_dat(nodeY, 0, edge2cells, 3, "double", OP_READ),
                op_arg_dat(nodeX, 1, edge2cells, 3, "double", OP_READ),
                op_arg_dat(nodeY, 1, edge2cells, 3, "double", OP_READ),
                op_arg_dat(q, 0, edge2cells, 4 * NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(q, 1, edge2cells, 4 * NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(exteriorQ, 0, edge2cells, 4 * 3 * NUM_FACE_PTS, "double", OP_WRITE),
                op_arg_dat(exteriorQ, 1, edge2cells, 4 * 3 * NUM_FACE_PTS, "double", OP_WRITE));
    // TODO handle boundary conditions

    // Calculate vectors F an G from q for each cell
    op_par_loop(euler_rhs, "euler_rhs", cells,
                op_arg_dat(q, -1, OP_ID, 4 * NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(exteriorQ, -1, OP_ID, 4 * 3 * NUM_FACE_PTS, "double", OP_READ),
                op_arg_dat(rx, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(ry, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(sx, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(sy, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(J, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_READ),
                op_arg_dat(sJ, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_READ),
                op_arg_dat(nx, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_READ),
                op_arg_dat(ny, -1, OP_ID, 3 * NUM_FACE_PTS, "double", OP_READ),
                op_arg_dat(qRHS, -1, OP_ID, 4 * NUM_SOLUTION_PTS, "double", OP_WRITE));
  }

  // Save the solution
  op_fetch_data_hdf5_file(node_coords, "points.h5");
  op_fetch_data_hdf5_file(x, "points.h5");
  op_fetch_data_hdf5_file(y, "points.h5");

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
  free(J_data);
  free(rx_data);
  free(ry_data);
  free(sx_data);
  free(sy_data);
  free(nx_data);
  free(ny_data);
  free(sJ_data);
  free(q_data);
  free(qRHS_data);
  free(edge2node_data);
  free(edge2cell_data);
  free(bedge2node_data);
  free(bedge2cell_data);
  free(exteriorQ_data);
}
