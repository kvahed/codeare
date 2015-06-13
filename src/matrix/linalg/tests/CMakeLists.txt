include_directories(
        ${PROJECT_SOURCE_DIR}/src/core
        ${PROJECT_SOURCE_DIR}/src/matrix
        ${PROJECT_SOURCE_DIR}/src/matrix/simd
        ${PROJECT_SOURCE_DIR}/src/matrix/linalg
        ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
        ${PROJECT_SOURCE_DIR}/src/matrix/io)

list (APPEND COMMON_LIBS ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})

add_executable (t_eig t_eig.cpp)
target_link_libraries (t_eig ${COMMON_LIBS})
add_test (eig t_eig)

add_executable (t_svd t_svd.cpp)
target_link_libraries (t_svd ${COMMON_LIBS})
add_test (svd t_svd)

add_executable (t_pca t_pca.cpp)
target_link_libraries (t_pca ${COMMON_LIBS})
add_test (pca t_pca)

add_executable (t_gemv t_gemv.cpp)
target_link_libraries (t_gemv ${COMMON_LIBS})
add_test (gemv t_gemv)

add_executable (t_gemm t_gemm.cpp)
target_link_libraries (t_gemm ${COMMON_LIBS})
add_test (gemm t_gemm)

add_executable (t_dot t_dot.cpp)
target_link_libraries (t_dot ${COMMON_LIBS})
add_test (dot t_dot)

add_executable (t_chol t_chol.cpp)
target_link_libraries (t_chol ${COMMON_LIBS})
add_test (chol t_chol)

add_executable (t_inv t_inv.cpp)
target_link_libraries (t_inv ${COMMON_LIBS})
add_test (inv t_inv)

add_executable (t_pinv t_pinv.cpp)
target_link_libraries (t_pinv ${COMMON_LIBS})
add_test (pinv t_pinv)

add_executable (t_norm t_norm.cpp)
target_link_libraries (t_norm ${COMMON_LIBS})
add_test (norm t_norm)

