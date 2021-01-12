#ifndef __AIRFOIL_LOADMESH_H
#define __AIRFOIL_LOADMESH_H

#include <string>

void load_mesh(std::string filename, int *numNodes, int *numCells, double **coords, int **cgnsCells);

#endif
