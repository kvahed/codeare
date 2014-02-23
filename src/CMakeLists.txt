include_directories (
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

add_subdirectory (core)
add_subdirectory (matrix)
add_subdirectory (modules)
add_subdirectory (mongoose)
add_subdirectory (tinyxml)

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.hpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.hpp" @ONLY)
list (APPEND SOURCES GitSHA1.hpp)
list (APPEND SOURCES options.h options.cpp)

list (APPEND COMLIBS tinyxml core ${GSL_LIBRARIES}
  ${OPENSSL_LIBRARIES} hdf5 hdf5_cpp boost_thread boost_system) 
if (${MATLAB_FOUND})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    list (APPEND COMLIBS libmx libmat)           
  else ()
    list (APPEND COMLIBS ${MATLAB_LIBRARIES})           
  endif ()
endif ()
if (${ISMRMRD_FOUND})
  list (APPEND COMLIBS ${ISMRMRD_LIBRARIES})
endif ()
if (${NIFTI_FOUND})
  list (APPEND COMLIBS ${NIFTI_NIFTIIO_LIBRARY} 
    ${NIFTI_ZNZ_LIBRARY} z)
endif ()
if (${MATLAB_FOUND})
  include_directories (${MATLAB_INCLUDE_DIRS})
endif ()

### lclient ###
add_executable (lclient ${SOURCES} LocalConnector.hpp
  LocalConnector.cpp testclt.cpp) 
target_compile_definitions (lclient PRIVATE LOCAL)
target_link_libraries (lclient ${COMLIBS})
install (TARGETS lclient DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

if (${OMNIORB4_FOUND})

  add_custom_target (RRSModuleSK.cc 
    COMMAND ${OMNIORB4_IDL_COMPILER} -bcxx RRSModule.idl 
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src)

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    list (APPEND OORB_DEFINES __WIN32__ __x86__ _WIN32_WINNT=0x0400)
  endif()

  ### rclient ###
  add_executable (rclient ${SOURCES} RemoteConnector.hpp
    RemoteConnector.cpp testclt.cpp RRSModule.hh RRSModuleSK.cc) 
  target_compile_definitions (rclient PRIVATE REMOTE ${OORB_DEFINES})
  target_link_libraries (rclient ${COMLIBS} ${OMNIORB4_LIBRARIES})

  ### codeared ###
  add_executable (codeared ${SOURCES} codeared.cpp RRSModule.hh
    RRSModuleSK.cc ReconServant.hpp ReconServant.cpp)
  target_compile_definitions (codeared PRIVATE ${OORB_DEFINES})
  target_link_libraries (codeared mongoose ${COMLIBS}
    ${OMNIORB4_LIBRARIES})  
  
  install (TARGETS rclient codeared DESTINATION
    ${CMAKE_INSTALL_PREFIX}/bin) 

endif()

