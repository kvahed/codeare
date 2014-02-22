include_directories (${PROJECT_SOURCE_DIR}/src/matrix/io
  ${PROJECT_SOURCE_DIR}/src/core 
  ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic  
  ${PROJECT_SOURCE_DIR}/src/matrix/dwt
  ${PROJECT_SOURCE_DIR}/src/matrix/ft 
  ${PROJECT_SOURCE_DIR}/src/matrix/io
  ${PROJECT_SOURCE_DIR}/src/matrix/linalg
  ${PROJECT_SOURCE_DIR}/src/matrix/interp
  ${FFTW3_INCLUDE_DIR}
  if (${ITK_FOUND})
    ${ITK_INCLUDE_DIRS}
  endif()
  )

add_library (CompressedSensing MODULE CompressedSensing.hpp CompressedSensing.cpp)
set_target_properties(CompressedSensing PROPERTIES PREFIX "")
target_link_libraries (CompressedSensing tinyxml core
  ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${FFTW3_LIBRARIES}
  ${OPENSSL_LIBRARIES}) 

add_library (KTPoints MODULE KTPoints.hpp KTPoints.cpp)
set_target_properties(KTPoints PROPERTIES PREFIX "")
target_link_libraries (KTPoints tinyxml core ${BLAS_LINKER_FLAGS}
  ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
  ${OPENSSL_LIBRARIES}) 

add_library (SENSE MODULE SENSE.hpp SENSE.cpp)
set_target_properties(SENSE PROPERTIES PREFIX "")
target_link_libraries (SENSE tinyxml core ${BLAS_LINKER_FLAGS}
  ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
  ${FFTW3_LIBRARIES} ${OPENSSL_LIBRARIES})

add_library (DummyRecon MODULE DummyRecon.hpp DummyRecon.cpp)
set_target_properties(DummyRecon PROPERTIES PREFIX "")
target_link_libraries (DummyRecon tinyxml core ${OPENSSL_LIBRARIES})

install (TARGETS CompressedSensing KTPoints SENSE DummyRecon
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 

if (${NFFT_FOUND})
  add_library (NuFFT MODULE NuFFT.hpp NuFFT.cpp)
  set_target_properties(NuFFT PROPERTIES PREFIX "")
  target_link_libraries (NuFFT tinyxml core ${NFFT_LIBRARIES}
    ${OPENSSL_LIBRARIES}) 
  add_library (CGSENSE MODULE CGSENSE.hpp CGSENSE.cpp)
  set_target_properties(CGSENSE PROPERTIES PREFIX "")
  target_link_libraries (CGSENSE tinyxml core ${NFFT_LIBRARIES}
    ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${OPENSSL_LIBRARIES})  
  install (TARGETS NuFFT CGSENSE DESTINATION
    ${CMAKE_INSTALL_PREFIX}/lib)
endif()

if (${ITK_FOUND})
  add_library (RelativeSensitivities MODULE RelativeSensitivities.cpp
    set_target_properties(RelativeSensitivities PROPERTIES PREFIX "")
    RelativeSensitivities.hpp)  
  target_link_libraries (RelativeSensitivities -L${ITK_LIBRARY_DIRS}
    ${ITK_LIBRARIES} ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES}
    ${OPENSSL_LIBRARIES})   
  install (TARGETS CompressedSensing DESTINATION
    ${CMAKE_INSTALL_PREFIX}/lib) 
endif()

