if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
endif()


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
target_link_libraries (CompressedSensing tinyxml core
  ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${FFTW3_LIBRARIES})

add_library (KTPoints MODULE KTPoints.hpp KTPoints.cpp)
target_link_libraries (KTPoints tinyxml core ${BLAS_LINKER_FLAGS}
  ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} )

add_library (SENSE MODULE SENSE.hpp SENSE.cpp)
target_link_libraries (SENSE tinyxml core ${BLAS_LINKER_FLAGS}
  ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
  ${FFTW3_LIBRARIES} )

add_library (DummyRecon MODULE DummyRecon.hpp DummyRecon.cpp)
target_link_libraries (DummyRecon tinyxml core )

if (${NFFT_FOUND})
  add_library (NuFFT MODULE NuFFT.hpp NuFFT.cpp)
  target_link_libraries (NuFFT tinyxml core ${NFFT_LIBRARIES})
  
  add_library (CGSENSE MODULE CGSENSE.hpp CGSENSE.cpp)
  target_link_libraries (CGSENSE tinyxml core ${NFFT_LIBRARIES}
    ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES}) 
endif()

if (${ITK_FOUND})
  add_library (RelativeSensitivities MODULE RelativeSensitivities.cpp
  RelativeSensitivities.hpp)  
  target_link_libraries (RelativeSensitivities -L${ITK_LIBRARY_DIRS}
  ${ITK_LIBRARIES} ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})  
endif()

