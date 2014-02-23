include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/src/core
  ${PROJECT_SOURCE_DIR}/src/matrix
  ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
  ${PROJECT_SOURCE_DIR}/src/matrix/dwt
  ${PROJECT_SOURCE_DIR}/src/matrix/ft
  ${PROJECT_SOURCE_DIR}/src/matrix/io
  ${PROJECT_SOURCE_DIR}/src/matrix/interp
  ${PROJECT_SOURCE_DIR}/src/matrix/linalg
  ${PROJECT_SOURCE_DIR}/src/matrix/simd
  ${PROJECT_SOURCE_DIR}/src/misc
  ${PROJECT_SOURCE_DIR}/src/mongoose
  ${PROJECT_SOURCE_DIR}/src/tests
  ${PROJECT_SOURCE_DIR}/src/tinyxml
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

### lclient ###
add_executable(lclient LocalConnector.hpp LocalConnector.cpp
  testclt.cpp options.cpp) 
target_compile_definitions(lclient PRIVATE LOCAL)
target_link_libraries (lclient core tinyxml ${GSL_LIBRARIES}
  ${OPENSSL_LIBRARIES} hdf5 hdf5_cpp ${BOOST_THREAD_LIBRARY}
  ${BOOST_SYSTEM_LIBRARY}) 
if (${MATLAB_FOUND})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    target_link_libraries (lclient libmx libmat)
  else()           
    target_link_libraries (lclient ${MATLAB_LIBRARIES})
  endif()
endif()
if(${ISMRMRD_FOUND})
  target_link_libraries (lclient ${ISMRMRD_LIBRARIES})
endif()
if (${NIFTI_FOUND})
  target_link_libraries (lclient ${NIFTI_NIFTIIO_LIBRARY}
  ${NIFTI_ZNZ_LIBRARY} z) 
endif()
if(${MATLAB_FOUND})
  include_directories (${MATLAB_INCLUDE_DIRS})
endif()

install (TARGETS lclient DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

if (${OMNIORB4_FOUND})

  add_custom_target (RRSModuleSK.cc 
    COMMAND ${OMNIORB4_IDL_COMPILER} -bcxx RRSModule.idl 
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src)

  ### rclient ###
  add_executable(rclient RemoteConnector.hpp RemoteConnector.cpp
    testclt.cpp RRSModule.hh RRSModuleSK.cc options.cpp) 
  add_dependencies(rclient RRSModuleSK.cc)
  target_compile_definitions(rclient PRIVATE REMOTE)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    target_compile_definitions(rclient PRIVATE __WIN32__ __x86__ _WIN32_WINNT=0x0400)
  endif()
  target_link_libraries (rclient core tinyxml ${GSL_LIBRARIES}
    ${OPENSSL_LIBRARIES} hdf5 hdf5_cpp ${BOOST_THREAD_LIBRARY}
    ${BOOST_SYSTEM_LIBRARY}) 

  ### codeared ###
  add_executable(codeared codeared.cpp RRSModuleSK.cc RRSModule.hh
  ReconServant.hpp ReconServant.cpp options.cpp)  
  add_dependencies(codeared RRSModuleSK.cc)
  target_compile_definitions(codeared PRIVATE __WIN32__ __x86__
  _WIN32_WINNT=0x0400)
  target_link_libraries (codeared core mongoose tinyxml
    ${GSL_LIBRARIES} ${OPENSSL_LIBRARIES} hdf5 hdf5_cpp
    ${BOOST_THREAD_LIBRARY} ${BOOST_SYSTEM_LIBRARY}) 

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    target_compile_definitions(rclient PRIVATE REMOTE __WIN32__
    __x86__ _WIN32_WINNT=0x0400)
    target_compile_definitions(codeared PRIVATE __WIN32__ __x86__
    _WIN32_WINNT=0x0400)
    target_link_libraries (rclient ws2_32 mswsock advapi32
    ${OMNIORB4_LIBRARY_omniORB4} ${OMNIORB4_LIBRARY_omniDynamic4}
    ${OMNIORB4_LIBRARY_omnithread}) 
    target_link_libraries (codeared ws2_32 mswsock advapi32
    ${OMNIORB4_LIBRARY_omniORB4} ${OMNIORB4_LIBRARY_omniDynamic4}
    ${OMNIORB4_LIBRARY_omnithread}) 
  else()
    target_link_libraries (rclient ${OMNIORB4_LIBRARIES})
    target_link_libraries (codeared ${OMNIORB4_LIBRARIES})
  endif()
  if (${MATLAB_FOUND})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
      target_link_libraries (rclient libmx libmat)
      target_link_libraries (codeared libmx libmat)
    else()           
      target_link_libraries (rclient ${MATLAB_LIBRARIES})
      target_link_libraries (codeared ${MATLAB_LIBRARIES})
    endif()
  endif()
  if(${ISMRMRD_FOUND})
    target_link_libraries (rclient ${ISMRMRD_LIBRARIES})
    target_link_libraries (codeared ${ISMRMRD_LIBRARIES})
  endif()
  if (${NIFTI_FOUND})
    target_link_libraries (rclient ${NIFTI_NIFTIIO_LIBRARY}
      ${NIFTI_ZNZ_LIBRARY} z) 
    target_link_libraries (codeared ${NIFTI_NIFTIIO_LIBRARY}
      ${NIFTI_ZNZ_LIBRARY} z) 
  endif()

  install (TARGETS rclient codeared DESTINATION
    ${CMAKE_INSTALL_PREFIX}/bin) 

endif()

