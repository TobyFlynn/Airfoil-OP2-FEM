op2_application(airfoil_seq LIBS op2_seq SOURCES airfoil.cpp load_mesh.cpp)

# Link CGNS Library
target_include_directories(airfoil_seq PUBLIC ${CGNS_DIR}/include)
target_link_libraries(airfoil_seq ${CGNS_DIR}/lib/libcgns.so)
