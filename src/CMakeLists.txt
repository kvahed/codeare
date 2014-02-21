include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/src/matrix
  ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
  ${PROJECT_SOURCE_DIR}/src/matrix/dwt
  ${PROJECT_SOURCE_DIR}/src/matrix/ft
  ${PROJECT_SOURCE_DIR}/src/matrix/io
  ${PROJECT_SOURCE_DIR}/src/matrix/interp
  ${PROJECT_SOURCE_DIR}/src/matrix/linalg
  ${PROJECT_SOURCE_DIR}/src/matrix/simd
  ${PROJECT_SOURCE_DIR}/src/tinyxml
  ${PROJECT_SOURCE_DIR}/src/misc
  ${PROJECT_SOURCE_DIR}/src/core
  ${PROJECT_SOURCE_DIR}/src/tests
  ${PROJECT_SOURCE_DIR}/src/modules
  ${GSL_INCLUDE_DIR}
  if (${OMNIORB4_FOUND})
    ${OMNIORB4_INCLUDE_DIR}
  endif()
  )

add_subdirectory(tinyxml)
add_subdirectory(core)
add_subdirectory(matrix)
add_subdirectory(mongoose)
add_subdirectory(modules)

add_executable(lclient LocalConnector.hpp LocalConnector.cpp
  testclt.cpp options.cpp) 
target_compile_definitions(lclient PRIVATE LOCAL)
target_link_libraries (lclient core tinyxml ${GSL_LIBRARIES}
  ${OPENSSL_LIBRARIES} hdf5 hdf5_cpp boost_thread boost_system)
if (${MATLAB_FOUND})
  target_link_libraries (lclient ${MATLAB_LIBRARIES})
endif()
if(${ISMRMRD_FOUND})
  target_link_libraries (lclient ${ISMRMRD_LIBRARIES})
endif()
if (${NIFTI_FOUND})
  target_link_libraries (lclient ${NIFTI_NIFTIIO_LIBRARY} ${NIFTI_ZNZ_LIBRARY} z)
endif()

#add_executable(rclient LocalConnector.hpp RemoteConnector.cpp
#  testclt.cpp options.cpp) 
#target_compile_definitions(rclient PRIVATE REMOTE)
#target_link_libraries (rclient core tinyxml ${GSL_LIBRARIES}
#  ${OPENSSL_LIBRARIES} hdf5 hdf5_cpp ${OMNIORB4_LIBRARIES}) 

#add_executable(codeared LocalConnector.hpp RemoteConnector.cpp
#  testclt.cpp options.cpp) 
#target_link_libraries (codeared core tinyxml ${GSL_LIBRARIES}
#  ${OPENSSL_LIBRARIES} hdf5 hdf5_cpp ${OMNIORB4_LIBRARIES}) 
