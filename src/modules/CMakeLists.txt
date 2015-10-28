include_directories (${PROJECT_SOURCE_DIR}/src/matrix/io
  ${PROJECT_SOURCE_DIR}/src/core 
  ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic  
  ${PROJECT_SOURCE_DIR}/src/matrix/dwt
  ${PROJECT_SOURCE_DIR}/src/matrix/ft 
  ${PROJECT_SOURCE_DIR}/src/matrix/io
  ${PROJECT_SOURCE_DIR}/src/matrix/linalg
  ${PROJECT_SOURCE_DIR}/src/matrix/interp
  ${PROJECT_SOURCE_DIR}/src/matrix/curves
  ${PROJECT_SOURCE_DIR}/src/optimisation
  ${FFTW3_INCLUDE_DIR}
  if (${ITK_FOUND})
    ${ITK_INCLUDE_DIRS}
  endif()
  )

if (${MACOSX})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flat_namespace") 
endif()

list (APPEND COMLIBS tinyxml core ${BLAS_LINKER_FLAGS}
  ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
  ${FFTW3_LIBRARIES} ${OPENSSL_LIBRARIES} ${Boost_CHRONO_LIBRARY}
  ${Boost_SYSTEM_LIBRARY} ${Boost_CHRONO_LIBRARY}
  ${Boost_TIMER_LIBRARY}  hdf5 hdf5_cpp)
if (${MATLAB_FOUND})
  list (APPEND COMLIBS ${MATLAB_LIBRARIES})
endif ()

add_library (KTPoints MODULE KTPoints.hpp KTPoints.cpp)
set_target_properties(KTPoints PROPERTIES PREFIX "")
target_link_libraries (KTPoints ${COMLIBS}) 

add_library (SENSE MODULE SENSE.hpp SENSE.cpp)
set_target_properties(SENSE PROPERTIES PREFIX "")
target_link_libraries (SENSE ${COMLIBS})

add_library (GRAPPA MODULE GRAPPA.hpp GRAPPA.cpp)
set_target_properties(GRAPPA PROPERTIES PREFIX "")
target_link_libraries (GRAPPA ${COMLIBS})

add_library (CoilCompression MODULE CoilCompression.hpp CoilCompression.cpp)
set_target_properties(CoilCompression PROPERTIES PREFIX "")
target_link_libraries (CoilCompression ${COMLIBS})

add_library (DummyRecon MODULE DummyRecon.hpp DummyRecon.cpp)
set_target_properties(DummyRecon PROPERTIES PREFIX "")
target_link_libraries (DummyRecon ${COMLIBS})

list (APPEND INST_TARGETS CompressedSensing KTPoints SENSE GRAPPA CoilCompression DummyRecon)

if (${NFFT3_FOUND})
  list (APPEND COMLIBS)

  add_library (NuFFT MODULE NuFFT.hpp NuFFT.cpp)
  set_target_properties(NuFFT PROPERTIES PREFIX "")
  target_link_libraries (NuFFT ${COMLIBS} ${NFFT3_LIBRARIES}) 

  add_library (CGSENSE MODULE CGSENSE.hpp CGSENSE.cpp)
  set_target_properties(CGSENSE PROPERTIES PREFIX "")
  target_link_libraries (CGSENSE ${COMLIBS} ${NFFT3_LIBRARIES})   

  add_library (GRASP MODULE GRASP.hpp GRASP.cpp)
  set_target_properties(GRASP PROPERTIES PREFIX "")
  target_link_libraries (GRASP ${COMLIBS} ${NFFT3_LIBRARIES})   
  
  add_library (EstimateSensitivities MODULE EstimateSensitivities.hpp
    EstimateSensitivities.cpp) 
  set_target_properties(EstimateSensitivities PROPERTIES PREFIX "")
  target_link_libraries (EstimateSensitivities ${COMLIBS}
    ${NFFT3_LIBRARIES} ${GSL_LIBRARIES})
  
  add_library (EstimateSensitivitiesRadialVIBE MODULE
    EstimateSensitivitiesRadialVIBE.hpp
    EstimateSensitivitiesRadialVIBE.cpp)  
  set_target_properties(EstimateSensitivitiesRadialVIBE PROPERTIES PREFIX "")
  target_link_libraries (EstimateSensitivitiesRadialVIBE ${COMLIBS}
    nfft3_threads nfft3f_threads ${GSL_LIBRARIES})
  list (APPEND INST_TARGETS EstimateSensitivitiesRadialVIBE)
  
  if (${HAS_CXX11_CONDITIONAL}) 
    add_library (XDGRASP MODULE XDGRASP.hpp XDGRASP.cpp)
    set_target_properties(XDGRASP PROPERTIES PREFIX "")
    target_link_libraries (XDGRASP ${COMLIBS} nfft3 nfft3f)
    list (APPEND INST_TARGETS XDGRASP)

    add_library (MotionDetectionXDGRASPLiver MODULE
      MotionDetectionXDGRASPLiver.hpp MotionDetectionXDGRASPLiver.cpp)  
    set_target_properties(MotionDetectionXDGRASPLiver PROPERTIES
      PREFIX "") 
    target_link_libraries (MotionDetectionXDGRASPLiver ${COMLIBS}
      ${GSL_LIBRARIES})
    list (APPEND INST_TARGETS MotionDetectionXDGRASPLiver)
  endif()

  list (APPEND INST_TARGETS NuFFT CGSENSE GRASP EstimateSensitivities)

endif()

add_library (CompressedSensing MODULE CompressedSensing.hpp
  CompressedSensing.cpp) 
set_target_properties(CompressedSensing PROPERTIES PREFIX "")
target_link_libraries (CompressedSensing ${COMLIBS}) 

if (${ITK_FOUND})
  add_library (RelativeSensitivities MODULE RelativeSensitivities.cpp
    RelativeSensitivities.hpp)
  set_target_properties(RelativeSensitivities PROPERTIES PREFIX "")
  target_link_libraries (RelativeSensitivities ${COMLIBS}
    ${ITK_LIBRARIES})
  list (APPEND INST_TARGETS RelativeSensitivities)
endif()

install (TARGETS ${INST_TARGETS} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 
