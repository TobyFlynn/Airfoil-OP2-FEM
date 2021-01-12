#define NUM_SOLUTION_PTS 15

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

using namespace std;

int main(int argc, char **argv) {
  double *coords;
  int *cgnsCells;
  int numNodes, numCells;

  load_mesh("./naca0012.cgns", &numNodes, &numCells, &coords, &cgnsCells);

  // Initialise OP2
  op_init(argc, argv, 2);

  // TODO get input from args
  int iter = 1;

  // Declare memory for data that will be calculated in initialisation kernel
  double *solution_pt_coords_x_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));
  double *solution_pt_coords_y_data = (double *)malloc(NUM_SOLUTION_PTS * numCells * sizeof(double));

  // Declare OP2 sets
  op_set nodes = op_decl_set(numNodes, "nodes");
  op_set cells = op_decl_set(numCells, "cells");

  // Declare OP2 maps
  op_map cell2nodes = op_decl_map(cells, nodes, 3, cgnsCells, "cell2nodes");

  // Declare OP2 datasets
    // Structure: {x, y}
  op_dat node_coords = op_decl_dat(nodes, 2, "double", coords, "node_coords");
    // Structure {x_0, x_1, ..}
  op_dat solution_pt_coords_x = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", solution_pt_coords_x_data, "solution_pt_coords_x");
    // Structure {y_0, y_1, ..}
  op_dat solution_pt_coords_y = op_decl_dat(cells, NUM_SOLUTION_PTS, "double", solution_pt_coords_y_data, "solution_pt_coords_y");
  //op_dat faceNormals = op_decl_dat(cells, 6, "double", NULL, "faceNormals");

  // Declare OP2 constants
  op_decl_const(NUM_SOLUTION_PTS, "double", ones);
  op_decl_const(NUM_SOLUTION_PTS, "double", solution_pts_r);
  op_decl_const(NUM_SOLUTION_PTS, "double", solution_pts_s);

  // Initialisation kernels
  op_par_loop(calc_solution_pt_coords, "calc_solution_pt_coords", cells,
              op_arg_dat(node_coords, 0, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(node_coords, 1, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(node_coords, 2, cell2nodes, 2, "double", OP_READ),
              op_arg_dat(solution_pt_coords_x, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE),
              op_arg_dat(solution_pt_coords_y, -1, OP_ID, NUM_SOLUTION_PTS, "double", OP_WRITE));


  // Run the simulation
  for(int i = 0; i < iter; i++) {

  }

  // Save the solution
  op_fetch_data_hdf5_file(node_coords, "points.h5");
  op_fetch_data_hdf5_file(solution_pt_coords_x, "points.h5");
  op_fetch_data_hdf5_file(solution_pt_coords_y, "points.h5");

  // Clean up OP2
  op_exit();

  free(coords);
  free(cgnsCells);
  free(solution_pt_coords_x_data);
  free(solution_pt_coords_y_data);
}
