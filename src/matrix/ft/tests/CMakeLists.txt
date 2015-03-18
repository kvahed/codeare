if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
endif()

include_directories(
        ${PROJECT_SOURCE_DIR}/src/core
        ${PROJECT_SOURCE_DIR}/src/matrix
        ${PROJECT_SOURCE_DIR}/src/matrix/simd
        ${PROJECT_SOURCE_DIR}/src/matrix/ft
        ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
        ${PROJECT_SOURCE_DIR}/src/matrix/io)


add_executable(t_dft t_dft.cpp)
target_link_libraries (t_dft ${FFTW3_LIBRARIES})
add_test(dft t_dft)

add_executable(t_fftshift t_fftshift.cpp)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  target_link_libraries (t_fftshift libmx libmat ${FFTW3_LIBRARIES} hdf5 hdf5_cpp)
else()
  target_link_libraries (t_fftshift mx mat ${FFTW3_LIBRARIES} ${HDF5_LIBRARIES})
endif()
add_test(fftshift t_fftshift)
set_tests_properties (fftshift PROPERTIES ENVIRONMENT
  "HDF5_DISABLE_VERSION_CHECK=2") 

