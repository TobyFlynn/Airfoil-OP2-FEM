// Include VTK stuff
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellIterator.h>
#include <vtkIdList.h>
// Include CGNS stuff
#include "cgnslib.h"
// Include C++ stuff
#include <string>
#include <iostream>
#include <map>
#include <memory>
#include <vector>
#include <utility>

using namespace std;

struct Point2D {
  double x;
  double y;
};

struct Cell {
  int points[3];
};

struct Edge {
  int points[2];
  int cells[2];
  int num[2];
};

int main(int argc, char **argv) {
  // Read in VTK file generated by Gmsh
  string fileName = "naca0012.vtk";
  vtkSmartPointer<vtkUnstructuredGrid> grid;
  auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName (fileName.c_str());
  reader->Update();
  grid = reader->GetOutput();

  map<int,unique_ptr<Point2D>> pointMap;
  map<int,unique_ptr<Cell>> cellMap;
  map<pair<int,int>,unique_ptr<Edge>> internalEdgeMap;
  // Maps for each type of boundary
  map<int,unique_ptr<Edge>> wallBoundaryEdgeMap;
  map<int,unique_ptr<Edge>> inflowBoundaryEdgeMap;
  map<int,unique_ptr<Edge>> outflowBoundaryEdgeMap;
  map<int,unique_ptr<Edge>> airfoilBoundaryEdgeMap;

  // Iterate over cells
  vtkSmartPointer<vtkCellIterator> cellIterator = grid->NewCellIterator();
  while(!cellIterator->IsDoneWithTraversal()) {
    // Check that this is a cell (not a line or point)
    // only an issue due to how Gmsh generates the VTK file
    if(cellIterator->GetNumberOfPoints() == 3) {
      vtkSmartPointer<vtkIdList> ids = cellIterator->GetPointIds();
      // Add points to pointMap if it does not already contain the point
      for(int i = 0; i < 3; i++) {
        if(pointMap.count(ids->GetId(i)) == 0) {
          double coords[3];
          grid->GetPoint(ids->GetId(i), coords);
          unique_ptr<Point2D> point = make_unique<Point2D>();
          point->x = coords[0];
          point->y = coords[1];
          pointMap.insert(pair<int,unique_ptr<Point2D>>(ids->GetId(i), move(point)));
        }
      }
      // Add cell to map
      unique_ptr<Cell> cell = make_unique<Cell>();
      // cell->points[0] = ids->GetId(0) + 1;
      // cell->points[1] = ids->GetId(1) + 1;
      // cell->points[2] = ids->GetId(2) + 1;
      // Convert from clockwise to anticlockwise
      cell->points[0] = ids->GetId(0) + 1;
      cell->points[1] = ids->GetId(2) + 1;
      cell->points[2] = ids->GetId(1) + 1;
      cellMap.insert(pair<int,unique_ptr<Cell>>(cellIterator->GetCellId(), move(cell)));
    }
    // Go to next cell
    cellIterator->GoToNextCell();
  }

  // Collate data
  vector<double> x;
  vector<double> y;
  for(auto const &point : pointMap) {
    x.push_back(point.second->x);
    y.push_back(point.second->y);
  }

  vector<cgsize_t> elements;
  for(auto const &elem : cellMap) {
    elements.push_back(elem.second->points[0]);
    elements.push_back(elem.second->points[1]);
    elements.push_back(elem.second->points[2]);
  }

  cout << "Number of points: " << x.size() << endl;
  cout << "VTK Number of points: " << grid->GetNumberOfPoints() << endl;
  cout << "Number of cell: " << elements.size() << endl;

  if(x.size() != grid->GetNumberOfPoints()) {
    cout << "******* Potential error when generating the grid *******" << endl;
    cout << "  Difference between number of points in grid and number in VTK file, indices may be wrong" << endl;
  }

  // Add edges to edge map if not already contained in mapping
  // If already added then update cell field
  // Try both combinations
  for(int i = 0; i < elements.size() / 3; i++) {
    int ind = i * 3;
    // Check that points are anticlockwise
    int p_1 = elements[ind];
    int p_2 = elements[ind + 1];
    int p_3 = elements[ind + 2];
    double val = (x[p_2 - 1] - x[p_1 - 1]) * (y[p_3 - 1] - y[p_2 - 1]) - (y[p_2 - 1] - y[p_1 - 1]) * (x[p_3 - 1] - x[p_2 - 1]);
    if(val < 0) {
      cout << "Switching points" << endl;
      cout << "  Old val: " << val << endl;
      elements[ind] = p_3;
      elements[ind + 2] = p_1;
      p_1 = elements[ind];
      p_2 = elements[ind + 1];
      p_3 = elements[ind + 2];
      val = (x[p_2 - 1] - x[p_1 - 1]) * (y[p_3 - 1] - y[p_2 - 1]) - (y[p_2 - 1] - y[p_1 - 1]) * (x[p_3 - 1] - x[p_2 - 1]);
      cout << "  New val: " << val << endl;
    }
    for(int j = 0; j < 3; j++) {
      int p1, p2;
      if(j == 0) {
        p1 = elements[ind];
        p2 = elements[ind + 1];
      } else if(j == 1) {
        p1 = elements[ind + 1];
        p2 = elements[ind + 2];
      } else {
        p1 = elements[ind + 2];
        p2 = elements[ind];
      }

      pair<int,int> key;
      if(p1 < p2) {
        key = make_pair(p1, p2);
      } else {
        key = make_pair(p2, p1);
      }

      if(internalEdgeMap.count(key) == 0) {
        unique_ptr<Edge> edge = make_unique<Edge>();
        edge->points[0] = key.first;
        edge->points[1] = key.second;
        edge->cells[0] = i;
        edge->cells[1] = -1;
        edge->num[0] = j;
        edge->num[1] = -1;
        internalEdgeMap.insert(make_pair(key, move(edge)));
      } else {
          if(internalEdgeMap.at(key)->cells[1] != -1) {
            cout << "ERROR in edge mapping: " << endl;
            cout << "  Old values: " << internalEdgeMap.at(key)->cells[0] << " " << internalEdgeMap.at(key)->cells[1] << endl;
            cout << "  New Value: " << i << endl;
            cout << "  Edges: " << internalEdgeMap.at(key)->points[0] << " " << internalEdgeMap.at(key)->points[1] << endl;
            cout << "  Key: " << key.first << " " << key.second << endl;
          }
          if(internalEdgeMap.at(key)->points[0] != key.first || internalEdgeMap.at(key)->points[1] != key.second) {
            cout << "ERROR in edge mapping: " << endl;
            cout << "  Prev Edge Nodes: " << internalEdgeMap.at(key)->points[0] << " " << internalEdgeMap.at(key)->points[1] << endl;
            cout << "  Current Nodes: " << p1 << " " << p2 << endl;
            cout << "  Key: " << key.first << " " << key.second << endl;
          }
          internalEdgeMap.at(key)->cells[1] = i;
          internalEdgeMap.at(key)->num[1] = j;
      }
    }
  }

  vector<int> edges;
  vector<int> boundaryEdges;
  for(auto const &edge : internalEdgeMap) {
    if(edge.second->cells[1] == -1) {
      boundaryEdges.push_back(edge.second->points[0]);
      boundaryEdges.push_back(edge.second->points[1]);
      boundaryEdges.push_back(edge.second->cells[0]);
      boundaryEdges.push_back(edge.second->num[0]);
    } else {
      if(edge.second->points[0] == edge.second->points[1])
        cout << "***** ERROR: Edge with identical points *****" << endl;
      if(edge.second->cells[0] == edge.second->cells[1])
        cout << "***** ERROR: Edge with identical cells *****" << endl;
      edges.push_back(edge.second->points[0]);
      edges.push_back(edge.second->points[1]);
      edges.push_back(edge.second->cells[0]);
      edges.push_back(edge.second->cells[1]);
      edges.push_back(edge.second->num[0]);
      edges.push_back(edge.second->num[1]);
    }
  }

  // cout << "Number of edges: " << edges.size() / 4 << endl;

  /*for(int i = 0; i < edges.size() / 4; i++) {
    int ind = i * 4;
    if(edges[ind + 2] == 2403 || edges[ind + 3] == 2403){// && edges[ind + 3] == 1676) {
      cout << "EDGE" << endl;
      cout << "  " << edges[ind] << " " << edges[ind + 1] << endl;
    }
  }*/

  // Write out CGNS file
  int file;
  if (cg_open("naca0012.cgns", CG_MODE_WRITE, &file)) {
    cg_error_exit();
  }

  // Create base
  string baseName = "Base";
  int cellDim = 2;
  int physicalDim = 2;
  int baseIndex;
  cg_base_write(file, baseName.c_str(), cellDim, physicalDim, &baseIndex);
  // Create zone
  string zoneName = "Zone1";
  cgsize_t sizes[3];
  // Number of vertices
  sizes[0] = pointMap.size();
  // Number of cells
  sizes[1] = cellMap.size();
  // Number of boundary vertices (zero if elements not sorted)
  sizes[2] = 0;
  int zoneIndex;
  cg_zone_write(file, baseIndex, zoneName.c_str(), sizes,
                CGNS_ENUMV(Unstructured), &zoneIndex);
  // Write grid coordinates
  int coordIndex;
  cg_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble),
                 "CoordinateX", x.data(), &coordIndex);
  cg_coord_write(file, baseIndex, zoneIndex, CGNS_ENUMV(RealDouble),
                 "CoordinateY", y.data(), &coordIndex);
  // Write elements
  int sectionIndex;
  int start = 1;
  int end = sizes[1];
  cg_section_write(file, baseIndex, zoneIndex, "GridElements", CGNS_ENUMV(TRI_3),
                   start, end, 0, elements.data(), &sectionIndex);
  // Write edges
  // {p1, p2, c1, c2, num1, num2}
  int numEdges = edges.size() / 6;
  cout << "Number of edges: " << numEdges << endl;
  cgsize_t dim[2] = {6, numEdges};
  cg_gopath(file, "/Base/Zone1");
  cg_user_data_write("Edges");
  cg_gopath(file, "/Base/Zone1/Edges");
  cg_array_write("EdgesData", CGNS_ENUMV(Integer), 2, dim, edges.data());

  // Write boundary edges
  // {p1, p2, c1, num}
  int numBoundaryEdges = boundaryEdges.size() / 4;
  cout << "Number of boundary edges: " << numBoundaryEdges << endl;
  cgsize_t boundaryDim[2] = {4, numBoundaryEdges};
  cg_gopath(file, "/Base/Zone1");
  cg_user_data_write("BoundaryEdges");
  cg_gopath(file, "/Base/Zone1/BoundaryEdges");
  cg_array_write("BoundaryEdgesData", CGNS_ENUMV(Integer), 2, boundaryDim, boundaryEdges.data());

  cg_close(file);
}
