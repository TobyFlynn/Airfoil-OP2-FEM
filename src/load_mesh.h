#ifndef __AIRFOIL_LOADMESH_H
#define __AIRFOIL_LOADMESH_H

#include <string>

void load_mesh(std::string filename, int *numNodes, int *numCells,
               int *numEdges, int *numBoundaryEdges, double **coords,
               int **cgnsCells, int **edge2node, int **edge2cell, int **edgeNum,
               int **bedge2node, int **bedge2cell, int **bedgeNum,
               int **bedge_type);

#endif
