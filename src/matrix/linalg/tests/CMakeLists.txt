include_directories(
        ${PROJECT_SOURCE_DIR}/src/core
        ${PROJECT_SOURCE_DIR}/src/matrix
        ${PROJECT_SOURCE_DIR}/src/matrix/simd
        ${PROJECT_SOURCE_DIR}/src/matrix/linalg
        ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
        ${PROJECT_SOURCE_DIR}/src/matrix/io)

add_executable (t_eig t_eig.cpp)
target_link_libraries (t_eig ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
add_test (eig t_eig)

add_executable (t_svd t_svd.cpp)
target_link_libraries (t_svd ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
add_test (svd t_svd)

add_executable (t_gemv t_gemv.cpp)
target_link_libraries (t_gemv ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
add_test (gemv t_gemv)

add_executable (t_gemm t_gemm.cpp)
target_link_libraries (t_gemm ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
add_test (gemm t_gemm)

add_executable (t_dot t_dot.cpp)
target_link_libraries (t_dot ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
add_test (dot t_dot)

add_executable (t_chol t_chol.cpp)
target_link_libraries (t_chol ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
add_test (chol t_chol)

add_executable (t_inv t_inv.cpp)
target_link_libraries (t_inv ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
add_test (inv t_inv)

add_executable (t_pinv t_pinv.cpp)
target_link_libraries (t_pinv ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
add_test (pinv t_pinv)

