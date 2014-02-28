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

set (prefix ${CMAKE_INSTALL_PREFIX})
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.hpp" @ONLY) 
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/Loader.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Loader.hpp" @ONLY) 
list (APPEND SOURCES GitSHA1.hpp)
list (APPEND SOURCES options.h options.cpp core/Workspace.hpp
  core/Workspace.cpp) 

if (WINDOWS)
  list (APPEND COMLIBS hdf5 hdf5_cpp)
else ()
  list (APPEND COMLIBS ${HDF5_LIBRARIES})
endif()

list (APPEND COMLIBS core tinyxml ${GSL_LIBRARIES}
  ${OPENSSL_LIBRARIES} ${Boost_THREAD_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}) 

if (${MATLAB_FOUND})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    list (APPEND COMLIBS libmx libmat)           
  else ()
    list (APPEND COMLIBS mx mat)
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

  set (OMNIORB_GENERATED ${CMAKE_CURRENT_SOURCE_DIR}/RRSModule.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/RRSModuleSK.cc) 

  add_custom_command (
    OUTPUT ${OMNIORB_GENERATED} 
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${OMNIORB4_IDL_COMPILER} -bcxx
    ${CMAKE_CURRENT_SOURCE_DIR}/RRSModule.idl) 
  
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    list (APPEND OORB_DEFINES __WIN32__ __x86__ _WIN32_WINNT=0x0400)
    list (APPEND OMNIORB4_LIBRARIES ws2_32 mswsock advapi32)
  endif()
  
  ### rclient ###
  add_executable (rclient ${OMNIORB_GENERATED} ${SOURCES}
    RemoteConnector.hpp RemoteConnector.cpp testclt.cpp) 
  target_compile_definitions (rclient PRIVATE REMOTE ${OORB_DEFINES})
  target_link_libraries (rclient ${COMLIBS} ${OMNIORB4_LIBRARIES})
  
  ### codeared ###
  add_executable (codeared ${OMNIORB_GENERATED} ${SOURCES}
    codeared.cpp ReconServant.hpp ReconServant.cpp)
  target_compile_definitions (codeared PRIVATE REMOTE ${OORB_DEFINES})
  target_link_libraries (codeared mongoose ${COMLIBS}
    ${OMNIORB4_LIBRARIES})  
  
  install (TARGETS rclient codeared DESTINATION
    ${CMAKE_INSTALL_PREFIX}/bin) 
  
endif()

if (${LINUX})
  set (RT_LIBS_ENV LD_LIBRARY_PATH)
elseif (${MACOSX})
  set (RT_LIBS_ENV DYLD_LIBRARY_PATH)
elseif (${WINDOWS})
  set (RT_LIBS_ENV PATH)
endif()

if (${MATLAB_FOUND})
  set (RT_LIBS_DIR "${MATLAB_LIBRARY_DIRS}")
endif()

if (${ISMRMRD_FOUND})
  set (RT_LIBS_DIR "${RT_LIBS_DIR}:${ISMRMRD_INCLUDE_DIR}/../lib")
endif()

set (RT_LIBS_DIR "${RT_LIBS_DIR}:${CMAKE_INSTALL_PREFIX}/lib")

add_test (sense lclient -t SENSE -b ${PROJECT_SOURCE_DIR}/share/sense
  -c 3fold.xml)
set_tests_properties (sense PROPERTIES ENVIRONMENT
  "${RT_LIBS_ENV}=${RT_LIBS_DIR}") 

add_test (cgsense lclient -t CGSENSE -b
  ${PROJECT_SOURCE_DIR}/share/cgsense -c 3fold.xml)
set_tests_properties (cgsense PROPERTIES ENVIRONMENT
  "${RT_LIBS_ENV}=${RT_LIBS_DIR}") 

add_test (ktpoints lclient -t KTPoints -b
  ${PROJECT_SOURCE_DIR}/share/ktpoints -c 7t-8ch-8p.xml)
set_tests_properties (ktpoints PROPERTIES ENVIRONMENT
  "${RT_LIBS_ENV}=${RT_LIBS_DIR}")

add_test (compressedsensing lclient -t CompressedSensing -b
  ${PROJECT_SOURCE_DIR}/share/compressedsensing -c
  phantom_denoising.xml) 
set_tests_properties (compressedsensing PROPERTIES ENVIRONMENT
  "${RT_LIBS_ENV}=${RT_LIBS_DIR}")

add_test (nufft lclient -t NuFFT -b ${PROJECT_SOURCE_DIR}/share/nufft/
  -c 16shot.xml)
set_tests_properties (nufft PROPERTIES ENVIRONMENT
  "${RT_LIBS_ENV}=${RT_LIBS_DIR}") 

add_test (vdspiral lclient -t VDSpiral -b
  ${PROJECT_SOURCE_DIR}/share/gradients -c config_vdcs.xml)
set_tests_properties (vdspiral PROPERTIES ENVIRONMENT
  "${RT_LIBS_ENV}=${RT_LIBS_DIR}")  
