// Include CGNS stuff
#include "cgnslib.h"
// Include C++ stuff
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;


// Template functions for loading cells into an int array
// (cgsize_t can be long or int depending on how the library was compiled)
template<typename T>
void cgns_load_cells(int file, int baseIndex, int zoneIndex, int *cgnsCells, int numCells);

template<>
void cgns_load_cells<int>(int file, int baseIndex, int zoneIndex, int *cgnsCells, int numCells) {
  cgsize_t parentData;
  cg_elements_read(file, baseIndex, zoneIndex, 1, (cgsize_t *)cgnsCells, &parentData);
  // CGNS starts numbering from 1 but OP2 starts from 0
  transform(cgnsCells, cgnsCells + 3 * numCells, cgnsCells, [](int x) { return x - 1;});
}

template<>
void cgns_load_cells<long>(int file, int baseIndex, int zoneIndex, int *cgnsCells, int numCells) {
  vector<cgsize_t> cells(numCells * 3);
  cgsize_t parentData;
  cg_elements_read(file, baseIndex, zoneIndex, 1, cells.data(), &parentData);
  // CGNS starts numbering from 1 but OP2 starts from 0
  transform(cells.begin(), cells.end(), cgnsCells, [](long x) { return (int)x - 1;});
}

void load_mesh(string filename, int *numNodes, int *numCells, double **coords, int **cgnsCells) {
  // Read CGNS grid
  int file;
  if(cg_open(filename.c_str(), CG_MODE_READ, &file)) {
    cg_error_exit();
  }

  int baseIndex = 1;
  int zoneIndex = 1;
  cgsize_t cg_numNodes;
  char zoneName[33];
  // Get zone name and size
  cg_zone_read(file, baseIndex, zoneIndex, zoneName, &cg_numNodes);
  *numNodes = (int) cg_numNodes;
  cout << "Zone Name: " << zoneName << endl;
  cout << "Zone Size: " << *numNodes << endl;

  // Get vertices
  vector<double> x(*numNodes);
  vector<double> y(*numNodes);
  cgsize_t minVertex = 1;
  cgsize_t maxVertex = cg_numNodes;
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateX",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, x.data());
  cg_coord_read(file, baseIndex, zoneIndex, "CoordinateY",
                CGNS_ENUMV(RealDouble), &minVertex, &maxVertex, y.data());

  *coords = (double *)malloc(2 * *numNodes * sizeof(double));
  for(int i = 0; i < x.size(); i++) {
    (*coords)[2 * i] = x[i];
    (*coords)[2 * i + 1] = y[i];
  }

  // Get number of sections
  int numSections;
  cg_nsections(file, baseIndex, zoneIndex, &numSections);
  cout << "Number of sections: " << numSections << endl << endl;

  // Get cell section
  char sectionName[33];
  CGNS_ENUMT(ElementType_t) elementType;
  cgsize_t elementStart, elementEnd;
  int elementNumBoundary, parentFlag;
  cg_section_read(file, baseIndex, zoneIndex, 1, sectionName, &elementType,
                  &elementStart, &elementEnd, &elementNumBoundary, &parentFlag);
  cout << "Section 1: " << sectionName << endl;
  cout << "Element Type: " << ElementTypeName[elementType] << endl;
  cout << "Start: " << elementStart << " End: " << elementEnd << endl;
  cout << "Element Number Boundary: " << elementNumBoundary;
  cout << " Parent Flag: " << parentFlag << endl;

  // Get cells
  *numCells = elementEnd - elementStart + 1;
  *cgnsCells = (int *)malloc(*numCells * 3 * sizeof(int));
  cgns_load_cells<cgsize_t>(file, baseIndex, zoneIndex, *cgnsCells, *numCells);
  // cout << "Cell data" << endl;
  // for(int i = 0; i < *numCells; i++) {
  //   cout << cgnsCells[3 * i] << " " << cgnsCells[3 * i + 1] << " " << cgnsCells[3 * i + 2] << endl;
  // }

  cg_close(file);
}
