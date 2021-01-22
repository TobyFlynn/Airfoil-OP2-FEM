set(KERNEL_HEADERS init_grid.h euler_rhs.h get_neighbour_q.h get_bedge_q.h set_ic.h set_workingQ.h update_Q.h calc_dt.h l2_error.h filterQ.h neighbour_zero.h)

op2_application(vortex_seq LIBS op2_seq op2_hdf5 SOURCES vortex.cpp load_mesh.cpp save_solution.cpp ${KERNEL_HEADERS})

# Link CGNS Library
target_include_directories(vortex_seq PUBLIC ${CGNS_DIR}/include ${OPENBLAS_DIR}/include)
target_link_libraries(vortex_seq ${CGNS_DIR}/lib/libcgns.so -L${OPENBLAS_DIR}/lib -lopenblas)
