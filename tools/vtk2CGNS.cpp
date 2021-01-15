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
      // Check each edge for whether it is a boundary
      for(int j = 0; j < 3; j++) {

      }
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
    for(int j = 0; j < 3; j++) {
      int p1 = elements[ind + j];
      int p2 = elements[ind + ((j + 1) % 3)];
      pair<int,int> key1 = make_pair(p1, p2);
      pair<int,int> key2 = make_pair(p2, p1);
      if(internalEdgeMap.count(key1) == 0 && internalEdgeMap.count(key2) == 0) {
        unique_ptr<Edge> edge = make_unique<Edge>();
        edge->points[0] = p1;
        edge->points[1] = p2;
        edge->cells[0] = i;
        edge->cells[1] = -1;
        internalEdgeMap.insert(make_pair(key1, move(edge)));
      } else {
        if(internalEdgeMap.count(key1) != 0) {
          if(internalEdgeMap.at(key1)->cells[1] != -1) {
            cout << "ERROR in edge mapping: " << endl;
            cout << "  Old values: " << internalEdgeMap.at(key1)->cells[0] << " " << internalEdgeMap.at(key1)->cells[1] << endl;
            cout << "  New Value: " << i << endl;
            cout << "  Edges: " << internalEdgeMap.at(key1)->points[0] << " " << internalEdgeMap.at(key1)->points[1] << endl;
            cout << "  Key: " << key1.first << " " << key1.second << endl;
          }
          internalEdgeMap.at(key1)->cells[1] = i;
        } else {
          if(internalEdgeMap.at(key2)->cells[1] != -1) {
            cout << "ERROR in edge mapping"  << endl;
            cout << "  Old values: " << internalEdgeMap.at(key2)->cells[0] << " " << internalEdgeMap.at(key2)->cells[1] << endl;
            cout << "  New Value: " << i << endl;
            cout << "  Edges: " << internalEdgeMap.at(key2)->points[0] << " " << internalEdgeMap.at(key2)->points[1] << endl;
            cout << "  Key: " << key2.first << " " << key2.second << endl;
          }
          internalEdgeMap.at(key2)->cells[1] = i;
        }
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
    } else {
      edges.push_back(edge.second->points[0]);
      edges.push_back(edge.second->points[1]);
      edges.push_back(edge.second->cells[0]);
      edges.push_back(edge.second->cells[1]);
    }
  }

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
  // {p1, p2, c1, c2}
  int numEdges = edges.size() / 4;
  cgsize_t dim[2] = {4, numEdges};
  cg_gopath(file, "/Base/Zone1");
  cg_user_data_write("Edges");
  cg_gopath(file, "/Base/Zone1/Edges");
  cg_array_write("EdgesData", CGNS_ENUMV(Integer), 2, dim, edges.data());

  // Write boundary edges
  // {p1, p2, c1}
  int numBoundaryEdges = boundaryEdges.size() / 3;
  cgsize_t boundaryDim[2] = {3, numBoundaryEdges};
  cg_gopath(file, "/Base/Zone1");
  cg_user_data_write("BoundaryEdges");
  cg_gopath(file, "/Base/Zone1/BoundaryEdges");
  cg_array_write("BoundaryEdgesData", CGNS_ENUMV(Integer), 2, boundaryDim, boundaryEdges.data());

  cg_close(file);
}
