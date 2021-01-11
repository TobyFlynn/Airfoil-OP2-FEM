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
      cell->points[0] = ids->GetId(0) + 1;
      cell->points[1] = ids->GetId(1) + 1;
      cell->points[2] = ids->GetId(2) + 1;
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
  cout << "Number of cell: " << elements.size() << endl;

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
  cg_section_write(file, baseIndex, zoneIndex, "Elements", CGNS_ENUMV(TRI_3),
                   start, end, 0, elements.data(), &sectionIndex);
  cg_close(file);
}
