project(Airfoil-OP2-FEM)

# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

# Default OP2 location for me
if(NOT OP2_DIR)
  set(OP2_DIR "/home/u1717021/Code/PhD/OP2-Common/op2/c/build")
endif()

# Load OP2 configuration (Set OP2_DIR to the directory containing OP2Config.cmake)
find_package(OP2 REQUIRED PATHS
  ${OP2-APPS_SOURCE_DIR}/../../op2/c/build
  ${PROJECT_SOURCE_DIR}/../../../op2/c/build
  ${CMAKE_INSTALL_PREFIX}/lib/op2)

# Import compiler flags for all build types
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OP2_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OP2_C_FLAGS}")

include_directories(${OP2_INCLUDE_DIRS})
add_definitions(${OP2_USER_DEFINITIONS})

include(OP2Application)

# CGNS library for me
if(NOT CGNS_DIR)
  set(CGNS_DIR "/home/u1717021/Code/PhD/CGNS/build")
endif()

add_subdirectory(tools)
add_subdirectory(src)
