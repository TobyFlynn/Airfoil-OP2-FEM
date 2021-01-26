set(KERNEL_HEADERS init_grid.h euler_rhs.h get_neighbour_q.h get_bedge_q.h set_ic.h set_workingQ.h update_Q.h calc_dt.h neighbour_zero.h)

# op2_application(airfoil_seq LIBS op2_seq SOURCES airfoil.cpp load_mesh.cpp save_solution.cpp ${KERNEL_HEADERS})
#op2_application(airfoil_seq LIBS op2_seq SOURCES airfoil_op.cpp seq/airfoil_seqkernels.cpp load_mesh.cpp save_solution.cpp ${KERNEL_HEADERS})

# Link CGNS Library
#target_include_directories(airfoil_seq PUBLIC ${CGNS_DIR}/include ${OPENBLAS_DIR}/include)
#target_link_libraries(airfoil_seq ${CGNS_DIR}/lib/libcgns.so -L${OPENBLAS_DIR}/lib -lopenblas)

#op2_application(airfoil_openmp LIBS op2_openmp SOURCES airfoil_op.cpp openmp/airfoil_kernels.cpp load_mesh.cpp save_solution.cpp ${KERNEL_HEADERS})

# Link CGNS Library
#target_include_directories(airfoil_openmp PUBLIC ${CGNS_DIR}/include ${OPENBLAS_DIR}/include)
#target_link_libraries(airfoil_openmp ${CGNS_DIR}/lib/libcgns.so -L${OPENBLAS_DIR}/lib -lopenblas)

op2_application(airfoil_cuda LIBS op2_cuda SOURCES airfoil_op.cpp cuda/airfoil_kernels.cu load_mesh.cpp save_solution.cpp ${KERNEL_HEADERS})
target_compile_definitions(airfoil_cuda PUBLIC -DBUILD_CUDA)
# Link CGNS Library
target_include_directories(airfoil_cuda PUBLIC ${CGNS_DIR}/include)
target_link_libraries(airfoil_cuda ${CGNS_DIR}/lib/libcgns.so CUDA::cublas)
