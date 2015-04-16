include_directories(
        ${PROJECT_SOURCE_DIR}/src/core
        ${PROJECT_SOURCE_DIR}/src/matrix
        ${PROJECT_SOURCE_DIR}/src/matrix/simd
        ${PROJECT_SOURCE_DIR}/src/matrix/dwt
        ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
        ${PROJECT_SOURCE_DIR}/src/matrix/io)


add_executable (t_dwt t_dwt.cpp)
if (${MSVC})
  target_link_libraries (t_dwt ${FFTW3_LIBRARIES} hdf5 hdf5_cpp)
else()
  target_link_libraries (t_dwt ${FFTW3_LIBRARIES} ${HDF5_LIBRARIES})
endif()

add_test (dwt t_dwt)
