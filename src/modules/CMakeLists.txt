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

list (APPEND COMLIBS tinyxml core ${BLAS_LINKER_FLAGS}
  ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
  ${FFTW3_LIBRARIES} ${OPENSSL_LIBRARIES}) 

add_library (CompressedSensing MODULE CompressedSensing.hpp
  CompressedSensing.cpp) 
set_target_properties(CompressedSensing PROPERTIES PREFIX "")
target_link_libraries (CompressedSensing ${COMLIBS}) 

add_library (KTPoints MODULE KTPoints.hpp KTPoints.cpp)
set_target_properties(KTPoints PROPERTIES PREFIX "")
target_link_libraries (KTPoints ${COMLIBS}) 

add_library (SENSE MODULE SENSE.hpp SENSE.cpp)
set_target_properties(SENSE PROPERTIES PREFIX "")
target_link_libraries (SENSE ${COMLIBS})

add_library (DummyRecon MODULE DummyRecon.hpp DummyRecon.cpp)
set_target_properties(DummyRecon PROPERTIES PREFIX "")
target_link_libraries (DummyRecon ${COMLIBS})

list (APPEND INST_TARGETS CompressedSensing KTPoints SENSE DummyRecon)

if (${NFFT_FOUND})

  add_library (NuFFT MODULE NuFFT.hpp NuFFT.cpp)
  set_target_properties(NuFFT PROPERTIES PREFIX "")
  target_link_libraries (NuFFT ${COMLIBS} ${NFFT_LIBRARIES}) 

  add_library (CGSENSE MODULE CGSENSE.hpp CGSENSE.cpp)
  set_target_properties(CGSENSE PROPERTIES PREFIX "")
  target_link_libraries (CGSENSE ${COMLIBS} ${NFFT_LIBRARIES})   

  list (APPEND INST_TARGETS NuFFT CGSENSE)

endif()

if (${ITK_FOUND})

  add_library (RelativeSensitivities MODULE RelativeSensitivities.cpp
    RelativeSensitivities.hpp)
  set_target_properties(RelativeSensitivities PROPERTIES PREFIX "")
  target_link_libraries (RelativeSensitivities ${COMLIBS}
    ${ITK_LIBRARIES})

  list (APPEND INST_TARGETS RelativeSensitivities)
endif()

install (TARGETS ${INST_TARGETS} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 
