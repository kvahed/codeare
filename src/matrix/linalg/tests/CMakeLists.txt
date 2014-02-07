include_directories(
        ${PROJECT_SOURCE_DIR}/src/core
        ${PROJECT_SOURCE_DIR}/src/matrix
        ${PROJECT_SOURCE_DIR}/src/matrix/simd
        ${PROJECT_SOURCE_DIR}/src/matrix/linalg
        ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
        ${PROJECT_SOURCE_DIR}/src/matrix/io)

add_executable (t_eig t_eig.cpp)
target_link_libraries (t_eig ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})
