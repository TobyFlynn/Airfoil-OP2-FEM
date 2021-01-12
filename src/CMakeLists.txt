set(KERNEL_HEADERS calc_solution_pt_coords.h)

op2_application(airfoil_seq LIBS op2_seq SOURCES airfoil.cpp load_mesh.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(airfoil_seq PUBLIC ${CGNS_DIR}/include ${OPENBLAS_DIR}/include)
target_link_libraries(airfoil_seq ${CGNS_DIR}/lib/libcgns.so -L${OPENBLAS_DIR}/lib -lopenblas)
