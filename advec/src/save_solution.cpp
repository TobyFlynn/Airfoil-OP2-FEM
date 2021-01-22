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

void save_solution(string filename, int numPts, int numCells, double *q, int *cellMap, int gam) {
  vector<double> velX(numPts);
  vector<double> velY(numPts);
  //vector<double> rho(numPts);
  //vector<double> p(numPts);

  int qNode0Ind = 0;
  int qNode1Ind = 4;
  int qNode2Ind = 14;

  for(int i = 0; i < numCells; i++) {
    int qCellInd = i * 15;
    int node0 = cellMap[i * 3];
    int node1 = cellMap[i * 3 + 1];
    int node2 = cellMap[i * 3 + 2];

    /*rho[node0]  = q[qCellInd + 4 * qNode0Ind];
    velX[node0] = q[qCellInd + 4 * qNode0Ind + 1] / rho[node0];
    velY[node0] = q[qCellInd + 4 * qNode0Ind + 2] / rho[node0];
    p[node0]    = (gam - 1) * (q[qCellInd + 4 * qNode0Ind + 3] - 0.5 * (q[qCellInd + 4 * qNode0Ind + 1] * velX[node0] + q[qCellInd + 4 * qNode0Ind + 2] * velY[node0]));

    rho[node1]  = q[qCellInd + 4 * qNode1Ind];
    velX[node1] = q[qCellInd + 4 * qNode1Ind + 1] / rho[node1];
    velY[node1] = q[qCellInd + 4 * qNode1Ind + 2] / rho[node1];
    p[node1]    = (gam - 1) * (q[qCellInd + 4 * qNode1Ind + 3] - 0.5 * (q[qCellInd + 4 * qNode1Ind + 1] * velX[node1] + q[qCellInd + 4 * qNode1Ind + 2] * velY[node1]));

    rho[node2]  = q[qCellInd + 4 * qNode2Ind];
    velX[node2] = q[qCellInd + 4 * qNode2Ind + 1] / rho[node2];
    velY[node2] = q[qCellInd + 4 * qNode2Ind + 2] / rho[node2];
    p[node2]    = (gam - 1) * (q[qCellInd + 4 * qNode2Ind + 3] - 0.5 * (q[qCellInd + 4 * qNode2Ind + 1] * velX[node2] + q[qCellInd + 4 * qNode2Ind + 2] * velY[node2]));*/
    velX[node0] = q[qCellInd + qNode0Ind];
    // velY[node0] = q[qCellInd + 2 * qNode0Ind + 1];

    velX[node1] = q[qCellInd + qNode1Ind];
    // velY[node1] = q[qCellInd + 2 * qNode1Ind + 1];

    velX[node2] = q[qCellInd + qNode2Ind];
    // velY[node2] = q[qCellInd + 2 * qNode2Ind + 1];
  }

  // Write out CGNS file
  int file;
  if (cg_open(filename.c_str(), CG_MODE_MODIFY, &file)) {
    cg_error_exit();
  }

  int baseIndex = 1;
  int zoneIndex = 1;
  int flowIndex;
  // Create flow solution node
  cg_sol_write(file, baseIndex, zoneIndex, "FlowSolution", CGNS_ENUMV(Vertex), &flowIndex);

  // // Write density
  // int rhoIndex;
  // cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Density", rho.data(), &rhoIndex);
  //
  // // Write pressure
  // int pIndex;
  // cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "Pressure", p.data(), &pIndex);

  // Write velocity x
  int velXIndex;
  cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityX", velX.data(), &velXIndex);

  // Write velocity y
  // int velYIndex;
  // cg_field_write(file, baseIndex, zoneIndex, flowIndex, CGNS_ENUMV(RealDouble), "VelocityY", velY.data(), &velYIndex);

  float exp[5];
  cg_goto(file, baseIndex, "end");
  cg_dataclass_write(CGNS_ENUMV(Dimensional));
  cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Degree));

  // exp[0] = 1.0f; exp[1] = -3.0f; exp[2] = 0.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  // cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", rhoIndex, "end");
  // cg_exponents_write(CGNS_ENUMV(RealSingle), exp);
  //
  // exp[0] = 1.0f; exp[1] = -1.0f; exp[2] = -2.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  // cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", pIndex, "end");
  // cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", velXIndex, "end");
  cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  // exp[0] = 0.0f; exp[1] = 1.0f; exp[2] = -1.0f; exp[3] = 0.0f; exp[4] = 0.0f;
  // cg_goto(file, baseIndex, "Zone_t", zoneIndex, "FlowSolution_t", flowIndex, "DataArray_t", velYIndex, "end");
  // cg_exponents_write(CGNS_ENUMV(RealSingle), exp);

  cg_close(file);
}
