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
  int order = 4;
  int numSolutionPts = ((order + 1) * (order + 2)) / 2;

  // Declare OP2 sets
  op_set nodes = op_decl_set(numNodes, "nodes");
  op_set cells = op_decl_set(numCells, "cells");

  // Declare OP2 maps
  op_map cell2nodes = op_decl_map(cells, nodes, 3, cgnsCells, "cell2nodes");

  // Declare OP2 datasets
  op_dat nodeCoords = op_decl_dat(nodes, 2, "double", coords, "nodeCoords");
  //op_dat solutionPtCoords = op_decl_dat(cells, 2 * numSolutionPts, "double", NULL, "solutionPtCoords");
  //op_dat faceNormals = op_decl_dat(cells, 6, "double", NULL, "faceNormals");

  // Initialisation kernel

  // Run the simulation
  for(int i = 0; i < iter; i++) {

  }

  // Save the solution

  // Clean up OP2
  op_exit();

  free(coords);
  free(cgnsCells);
}
