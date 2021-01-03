// Include OP2 stuff
#include "op_seq.h"
// Include CGNS stuff
#include "cgnslib.h"
// Include C++ stuff
#include <string>
#include <iostream>
#include <memory>
#include <vector>

using namespace std;

int main(int argc, char **argv) {
  // Initialise OP2
  op_init(argc, argv, 2);

  // Read CGNS grid
  int file;
  if(cg_open("./naca0012.cgns", CG_MODE_READ, &file)) {
    cg_error_exit();
  }

  int baseIndex = 1;
  int zoneIndex = 1;
  int numNodes;
  char zoneName[33];
  // Get zone name and size
  cg_zone_read(file, baseIndex, zoneIndex, zoneName, &numNodes);
  cout << "Zone Name: " << zoneName << endl;
  cout << "Zone Size: " << numNodes << endl;

  // Get vertices
  vector<double> x(numNodes);
  vector<double> y(numNodes);
  int minVertex = 1;
  int maxVertex = numNodes;
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateX",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, x.data());
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateY",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, y.data());
  // Combine x and y vectors (necessary for OP2)
  vector<double> coords;
  for(int i = 0; i < x.size(); i++) {
    coords.push_back(x[i]);
    coords.push_back(y[i]);
  }

  // Get number of sections
  int numSections;
  cg_nsections(file, baseIndex, zoneIndex, &numSections);
  cout << "Number of sections: " << numSections << endl << endl;

  // Get cell section
  char sectionName[33];
  CGNS_ENUMT(ElementType_t) elementType;
  int elementStart, elementEnd, elementNumBoundary, parentFlag;
  cg_section_read(file, baseIndex, zoneIndex, 1, sectionName, &elementType,
                  &elementStart, &elementEnd, &elementNumBoundary, &parentFlag);
  cout << "Section 1: " << sectionName << endl;
  cout << "Element Type: " << ElementTypeName[elementType] << endl;
  cout << "Start: " << elementStart << " End: " << elementEnd << endl;
  cout << "Element Number Boundary: " << elementNumBoundary;
  cout << " Parent Flag: " << parentFlag << endl;

  // Get cells
  int numCells = elementEnd - elementStart + 1;
  vector<int> cgnsCells(numCells * 3);
  int parentData;
  cg_elements_read(file, baseIndex, zoneIndex, 1, cgnsCells.data(),
                   &parentData);

  cg_close(file);

  // Declare OP2 sets
  op_set nodes = op_decl_set(numNodes, "nodes");
  op_set cells = op_decl_set(numCells, "cells");

  // Declare OP2 maps
  op_map cell2nodes = op_decl_map(cells, nodes, 3, cgnsCells.data(),
                                  "cell2nodes");

  // Declare OP2 datasets
  op_dat nodeCoordinates = op_decl_dat(nodes, 2, "double", coords.data(),
                                       "nodeCoordinates");

  // Clean up OP2
  op_exit();
}
