if (POLICY CMP0021)
  cmake_policy(SET CMP0021 OLD)
endif()


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
  ${PROJECT_SOURCE_DIR}/src/tests
  ${PROJECT_SOURCE_DIR}/src/tinyxml
  ${PROJECT_SOURCE_DIR}/src/modules
  if (${GSL_FOUND})
	${GSL_INCLUDE_DIR}
  endif()
  if (${OMNIORB4_FOUND})
	${PROJECT_SOURCE_DIR}/src/mongoose
    ${OMNIORB4_INCLUDE_DIR}
  endif()
  )

add_subdirectory (core)
add_subdirectory (matrix)
add_subdirectory (modules)
add_subdirectory (tinyxml)

if (${OMNIORB4_FOUND})
  add_subdirectory (mongoose)
endif()

configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/config.h")

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
  set (LD_ENV "LD_LIBRARY_PATH=$LD_LIBRARY_PATH")
elseif (${MACOSX})
  set (LD_ENV "DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH")
elseif (${WINDOWS})
  set (LD_ENV "PATH=%PATH%")
endif()

if (${MATLAB_FOUND})
  if (${WINDOWS}) 
    set (LD_ENV "${LD_ENV}\;${MATLAB_LIBRARY_DIRS}")
  else ()
    set (LD_ENV "${LD_ENV}:${MATLAB_LIBRARY_DIRS}")
  endif ()
endif()

if (${ISMRMRD_FOUND})
  if (${WINDOWS}) 
    set (LD_ENV "${LD_ENV}\;${ISMRMRD_INCLUDE_DIR}/../lib")
  else ()
    set (LD_ENV "${LD_ENV}:${ISMRMRD_INCLUDE_DIR}/../lib")
  endif ()
endif()

if (${WINDOWS}) 
  set (LD_ENV
    "${LD_ENV}\;${CMAKE_CURRENT_BINARY_DIR}/modules")
  set (LD_ENV
    "${LD_ENV}\;${CMAKE_CURRENT_BINARY_DIR}/modules/Release")
  set (LD_ENV "${LD_ENV}\;${CMAKE_CURRENT_BINARY_DIR}/Release")
else ()
  set (LD_ENV
    "${LD_ENV}:${CMAKE_CURRENT_BINARY_DIR}/modules")
endif ()

if (${WINDOWS})
  string(REPLACE "/" "\\" LD_ENV "${LD_ENV}")
endif()

# Test macro -------------------------------------------------------- 
include (TestMacro)

set (TEST_CALL lclient -t SENSE -b
  ${PROJECT_SOURCE_DIR}/share/sense -c 3fold.xml)  
MP_TESTS ("sense" "${TEST_CALL}")

set (TEST_CALL lclient -t KTPoints -b
  ${PROJECT_SOURCE_DIR}/share/ktpoints -c 7t-8ch-8p.xml)  
MP_TESTS ("ktpoints" "${TEST_CALL}")

set (TEST_CALL lclient -t CompressedSensing -b
  ${PROJECT_SOURCE_DIR}/share/compressedsensing -c phantom_denoising.xml) 
MP_TESTS ("compressedsensing" "${TEST_CALL}")

if (${GSL_FOUND})
   set (TEST_CALL lclient -t VDSpiral -b
     ${PROJECT_SOURCE_DIR}/share/gradients -c config_vdcs.xml) 
   MP_TESTS ("vdspiral" "${TEST_CALL}")
endif()

if (${NFFT_FOUND})

  set (TEST_CALL lclient -t NuFFT -b ${PROJECT_SOURCE_DIR}/share/nufft/ -c
    16shot.xml) 
  MP_TESTS ("nufft" "${TEST_CALL}")

  set (TEST_CALL lclient -t CGSENSE -b
    ${PROJECT_SOURCE_DIR}/share/cgsense -c 3fold.xml) 
  MP_TESTS ("cgsense" "${TEST_CALL}")

endif()

