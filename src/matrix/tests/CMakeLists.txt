if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
endif()

include_directories(
        ${PROJECT_SOURCE_DIR}/src/core
        ${PROJECT_SOURCE_DIR}/src/matrix
        ${PROJECT_SOURCE_DIR}/src/matrix/simd
        ${PROJECT_SOURCE_DIR}/src/matrix/ft
        ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
        ${PROJECT_SOURCE_DIR}/src/matrix/io)

add_executable(t_colon t_colon.cpp)
add_test(colon t_colon)

add_executable(t_sort t_sort.cpp)
add_test(sort t_sort)

add_executable(t_tvop t_tvop.cpp)
if (${MSVC})
  target_link_libraries (t_tvop ${FFTW3_LIBRARIES} hdf5 hdf5_cpp)
else()
  target_link_libraries (t_tvop ${FFTW3_LIBRARIES} ${HDF5_LIBRARIES})
endif()
add_test(tvop t_tvop)

