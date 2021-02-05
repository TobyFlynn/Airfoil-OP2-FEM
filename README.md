# Airfoil-OP2-FEM
A Discontinuous Galerkin FEM CFD application. Uses SSPRK3 time integration.
Solves the same problem that the Airfoil sample application included in OP2 solves.

Dependencies:
- OP2
- CGNS
- OpenBLAS
- VTK

Directory structure:
- The 'src' directory contains the code for the Airfoil application.
- The 'tools' directory contains an unstructured VTK to CGNS mesh conversion code.
- The 'vortex' directory contains an application that simulates the Euler equations but for a problem with an analytical solution.
- The 'advec' directory contains a 2D advection code.

Build instructions:
```
mkdir build
cd build
cmake ..
make
```

CMakeLists.txt sets the locations of libraries to what I have on my system but you can change these by passing the following flags to the CMake command: `OP2_DIR`, `CGNS_DIR` and `OPENBLAS_DIR`.
