// Include CGNS stuff
#include "cgnslib.h"
// Include C++ stuff
#include <iostream>
#include <fstream>
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

void load_mesh(std::string filename, int *numNodes, int *numCells,
               int *numEdges, int *numBoundaryEdges, double **coords,
               int **cgnsCells, int **edge2node, int **edge2cell,
               int **bedge2node, int **bedge2cell, int **bedge_type) {
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

  // Get edge data
  cg_gopath(file, "/Base/Zone1/Edges");
  char arrayName[33];
  DataType_t arrayDataType;
  int arrayRank;
  cgsize_t arrayDims[2];
  cg_array_info(1, arrayName, &arrayDataType, &arrayRank, arrayDims);
  cout << "Array Name: " << arrayName << endl;
  cout << "Array Dims: " << arrayDims[0] << " " << arrayDims[1] << endl;
  *numEdges = arrayDims[1];
  *edge2node = (int *)malloc(2 * *numEdges * sizeof(int));
  *edge2cell = (int *)malloc(2 * *numEdges * sizeof(int));
  vector<int> edgeData(arrayDims[0] * arrayDims[1]);
  cg_array_read(1, edgeData.data());

  for(int i = 0; i < *numEdges; i++) {
    // - 1 as CGNS counts points from 1 but OP2 counts from 0
    // Cell index do start from one in this data
    (*edge2node)[i * 2]     = edgeData[i * 4] - 1;
    (*edge2node)[i * 2 + 1] = edgeData[i * 4 + 1] - 1;
    (*edge2cell)[i * 2]     = edgeData[i * 4 + 2];
    (*edge2cell)[i * 2 + 1] = edgeData[i * 4 + 3];
  }

  // Get boundary edge data
  cg_gopath(file, "/Base/Zone1/BoundaryEdges");
  char barrayName[33];
  DataType_t barrayDataType;
  int barrayRank;
  cgsize_t barrayDims[2];
  cg_array_info(1, barrayName, &barrayDataType, &barrayRank, barrayDims);
  cout << "Array Name: " << barrayName << endl;
  cout << "Array Dims: " << barrayDims[0] << " " << barrayDims[1] << endl;
  *numBoundaryEdges = barrayDims[1];
  *bedge2node = (int *)malloc(2 * *numBoundaryEdges * sizeof(int));
  *bedge2cell = (int *)malloc(*numBoundaryEdges * sizeof(int));
  *bedge_type = (int *)malloc(*numBoundaryEdges * sizeof(int));
  vector<int> bedgeData(barrayDims[0] * barrayDims[1]);
  cg_array_read(1, bedgeData.data());

  for(int i = 0; i < *numBoundaryEdges; i++) {
    // - 1 as CGNS counts points from 1 but OP2 counts from 0
    // Cell index do start from one in this data
    (*bedge2node)[i * 2]     = bedgeData[i * 3] - 1;
    (*bedge2node)[i * 2 + 1] = bedgeData[i * 3 + 1] - 1;
    (*bedge2cell)[i]         = bedgeData[i * 3 + 2];
    // TODO set boundary type
    // 0 = inflow, 1 = outflow, 2 = wall
  }

  cg_close(file);

  // Write cell data to text file (for python sanity check)
  ofstream cellFile;
  cellFile.open("cells.txt");
  for(int i = 0; i < *numCells; i++) {
    cellFile << (*cgnsCells)[3 * i] << " " << (*cgnsCells)[3 * i + 1] << " " << (*cgnsCells)[3 * i + 2] << endl;
  }
  cellFile.close();
}
