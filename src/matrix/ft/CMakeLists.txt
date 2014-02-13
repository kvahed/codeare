add_subdirectory(tests)

include_directories (${PROJECT_SOURCE_DIR}/src/matrix/io
  ${PROJECT_SOURCE_DIR}/src/core 
  ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic  
  ${PROJECT_SOURCE_DIR}/src/matrix/ft 
  ${PROJECT_SOURCE_DIR}/src/matrix/linalg
  ${CBLAS_INCLUDE_DIR}
  ${FFTW3_INCLUDE_DIR}
  if(${NFFT_FOUND})
    ${NFFT_INCLUDE_DIR}
  endif()
)

if (${NFFT_FOUND})
  add_library (codeare-ft SHARED DFT.hpp DFT.cpp CSENSE.hpp CSENSE.cpp CGRAPPA.hpp CGRAPPA.cpp NCSENSE.hpp NCSENSE.cpp)
  target_link_libraries (codeare-ft core ${CBLAS_LINKER_FLAGS} ${CBLAS_LIBRARIES} ${FFTW3_LIBRARIES} ${NFFT_LIBRARIES})
else ()
  add_library (codeare-ft SHARED DFT.hpp DFT.cpp CSENSE.hpp CSENSE.cpp )
  target_link_libraries (codeare-ft ${CBLAS_LINKER_FLAGS} ${CBLAS_LIBRARIES} ${FFTW3_LIBRARIES})
endif()
