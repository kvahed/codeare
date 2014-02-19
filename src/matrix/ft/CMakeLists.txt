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

if (${NFFT_FOUND})
  add_library (codeare-ft SHARED DFT.hpp DFT.cpp CSENSE.hpp CSENSE.cpp CGRAPPA.hpp CGRAPPA.cpp NCSENSE.hpp NCSENSE.cpp)
  target_link_libraries (codeare-ft core ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${FFTW3_LIBRARIES} ${NFFT_LIBRARIES})
else ()
  add_library (codeare-ft SHARED DFT.hpp DFT.cpp CSENSE.hpp CSENSE.cpp )
  target_link_libraries (codeare-ft core ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${FFTW3_LIBRARIES})
endif()
