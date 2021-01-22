set(KERNEL_HEADERS advec_rhs.h get_neighbour_q.h get_bedge_q.h set_ic.h set_workingQ.h update_Q.h init_grid.h neighbour_zero.h)

op2_application(advec_seq LIBS op2_seq op2_hdf5 SOURCES advec.cpp load_mesh.cpp save_solution.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(advec_seq PUBLIC ${CGNS_DIR}/include ${OPENBLAS_DIR}/include)
target_link_libraries(advec_seq ${CGNS_DIR}/lib/libcgns.so -L${OPENBLAS_DIR}/lib -lopenblas)
