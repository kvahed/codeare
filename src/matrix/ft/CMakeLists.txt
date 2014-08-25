add_subdirectory(tests)

include_directories (${PROJECT_SOURCE_DIR}/src/matrix/io
  ${PROJECT_SOURCE_DIR}/src/core 
  ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic  
  ${PROJECT_SOURCE_DIR}/src/matrix/ft 
  ${PROJECT_SOURCE_DIR}/src/matrix/linalg
  ${FFTW3_INCLUDE_DIR}
  if(${NFFT_FOUND})
    ${NFFT_INCLUDE_DIR}
  endif()
)
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW) # CMake 3.0
endif()

list (APPEND COMLIBS core ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${FFTW3_LIBRARIES} hdf5 hdf5_cpp)
if (${MATLAB_FOUND})
  list (APPEND COMLIBS ${MATLAB_MX_LIBRARY} ${MATLAB_ENG_LIBRARY})
endif ()

list (APPEND SOURCES DFT.hpp DFT.cpp CSENSE.hpp CSENSE.cpp CGRAPPA.cpp)

if (${NFFT_FOUND})
  list (APPEND SOURCES NCSENSE.hpp NCSENSE.cpp NFFT.hpp NFFT.cpp)
  list (APPEND COMLIBS ${NFFT_LIBRARIES})
endif ()

add_library (codeare-ft SHARED ${SOURCES})
target_link_libraries (codeare-ft ${COMLIBS})
