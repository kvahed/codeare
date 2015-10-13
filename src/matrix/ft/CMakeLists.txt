add_subdirectory(tests)

include_directories (${PROJECT_SOURCE_DIR}/src/matrix/io
  ${PROJECT_SOURCE_DIR}/src/core 
  ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic  
  ${PROJECT_SOURCE_DIR}/src/matrix/ft 
  ${PROJECT_SOURCE_DIR}/src/matrix/linalg
  ${PROJECT_SOURCE_DIR}/src/optimisation
  ${FFTW3_INCLUDE_DIR}
  if(${NFFT3_FOUND})
    ${NFFT3_INCLUDE_DIR}
  endif()
)
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW) # CMake 3.0
endif()

list (APPEND COMLIBS core ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${FFTW3_LIBRARIES} hdf5 hdf5_cpp boost_system)
if (${MATLAB_FOUND})
  list (APPEND COMLIBS ${MATLAB_LIBRARIES})
endif ()

list (APPEND SOURCES DFT.cpp CSENSE.cpp CGRAPPA.cpp)

if (${NFFT3_FOUND})
  list (APPEND SOURCES NCSENSE.cpp NFFT.cpp CS_XSENSE.cpp)
  list (APPEND COMLIBS ${NFFT3_LIBRARIES})
endif ()

add_library (codeare-ft SHARED ${SOURCES})
target_link_libraries (codeare-ft ${COMLIBS})

list (APPEND INST_TARGETS codeare-ft)
install (TARGETS ${INST_TARGETS} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 
