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
  ${PROJECT_SOURCE_DIR}/src/mongoose
  )

add_subdirectory (core)
add_subdirectory (matrix)
add_subdirectory (modules)
add_subdirectory (tinyxml)
add_subdirectory (optimisation)

if (${OMNIORB4_FOUND})
  add_subdirectory (mongoose)
endif()

configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/config.h")

set (prefix ${CMAKE_INSTALL_PREFIX})
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.hpp" @ONLY) 
list (APPEND SOURCES GitSHA1.hpp)
list (APPEND SOURCES options.h options.cpp) 

if (WINDOWS)
  list (APPEND COMLIBS hdf5 hdf5_cpp)
else ()
  list (APPEND COMLIBS ${HDF5_LIBRARIES} ${OPENSSL_LIBRARIES})
endif()

list (APPEND COMLIBS core tinyxml ${GSL_LIBRARIES}
  ${OPENSSL_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_CHRONO_LIBRARY}
  ${Boost_TIMER_LIBRARY})  

if (${Matlab_FOUND})
    list (APPEND COMLIBS ${Matlab_MX_LIBRARY} ${Matlab_MAT_LIBRARY})
endif ()
if (${ISMRMRD_FOUND})
  list (APPEND COMLIBS ${ISMRMRD_LIBRARIES})
endif ()
if (${NIFTI_FOUND})
  list (APPEND COMLIBS ${NIFTI_NIFTIIO_LIBRARY} 
    ${NIFTI_ZNZ_LIBRARY} z)
endif ()
if (${Matlab_FOUND})
  include_directories (${Matlab_INCLUDE_DIRS})
endif ()

set (NETWORKING_SOURCES "")

if (${OMNIORB4_FOUND})
  
  set (OMNIORB_GENERATED ${CMAKE_CURRENT_SOURCE_DIR}/RRSModule.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/RRSModuleSK.cc) 
  
  add_custom_command (
    OUTPUT ${OMNIORB_GENERATED} 
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${OMNIORB4_IDL_COMPILER} -bcxx
      ${CMAKE_CURRENT_SOURCE_DIR}/RRSModule.idl
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/RRSModule.idl) 
  
  if (${MSVC})
    list (APPEND OORB_DEFINES __WIN32__ __x86__ _WIN32_WINNT=0x0400)
    list (APPEND OMNIORB4_LIBRARIES ws2_32 mswsock advapi32)
  else()
    list (APPEND OORB_DEFINES "-DHAVE_OMNIORB4")
  endif()
  
  list (APPEND COMLIBS ${OMNIORB4_LIBRARIES})
  list (APPEND NETWORKING_SOURCES "RemoteConnector.hpp" "RemoteConnector.cpp")
  
  ### codeared ###
  add_executable (codeared ${OMNIORB_GENERATED} ${SOURCES}
    codeared.cpp codeared.hpp ReconServant.hpp ReconServant.cpp)
  target_compile_definitions (codeared PRIVATE ${OORB_DEFINES})
  target_link_libraries (codeared mongoose ${COMLIBS} ${Boost_THREAD_LIBRARY})  
  
  install (TARGETS codeared DESTINATION
	${CMAKE_INSTALL_PREFIX}/bin) 
  
endif()  

### codeare ###
add_executable (codeare ${OMNIORB_GENERATED} ${SOURCES}
  Connector.hpp LocalConnector.hpp LocalConnector.cpp
  ${NETWORKING_SOURCES} codeare.hpp codeare.cpp)  
target_compile_definitions (codeare PRIVATE "${OORB_DEFINES}")
target_link_libraries (codeare ${COMLIBS})

install (TARGETS codeare DESTINATION
  ${CMAKE_INSTALL_PREFIX}/bin)

if (${LINUX})
  set (LD_ENV $ENV{LD_LIBRARY_PATH})
elseif (${MACOSX})
  set (LD_ENV $ENV{DYLD_LIBRARY_PATH})
elseif (${WINDOWS})
  set (LD_ENV "PATH=%PATH%")
endif()

if (${Matlab_FOUND})
  if (${WINDOWS}) 
    set (LD_ENV "${LD_ENV}\;${Matlab_ROOT}/bin/win64")
  else ()
    set (LD_ENV "${LD_ENV}:${Matlab_LIBRARY_DIRS}")
  endif ()
endif()

if (${WINDOWS}) 
  set (LD_ENV "${LD_ENV}\;${Boost_LIBRARY_DIRS}")
  set (LD_ENV "${LD_ENV}\;${OPENSSL_ROOT_DIR}")
  set (LD_ENV "${LD_ENV}\;${FFTW3_INCLUDE_DIR}")
  set (LD_ENV "${LD_ENV}\;${NFFT_INCLUDE_DIR}/../bin")
else ()
  set (LD_ENV "${LD_ENV}:${Boost_LIBRARY_DIRS}")
endif ()

if (${ISMRMRD_FOUND})
  if (${WINDOWS}) 
    set (LD_ENV "${LD_ENV}\;${ISMRMRD_INCLUDE_DIR}/../lib")
  else ()
    set (LD_ENV "${LD_ENV}:${ISMRMRD_INCLUDE_DIR}/../lib")
  endif ()
endif()

if (${OMNIORB4_FOUND})
  if (${WINDOWS}) 
    set (LD_ENV "${LD_ENV}\;${OMNIORB4_BIN_DIR_WIN32}")
  endif()
endif()

if (${WINDOWS}) 
  list (APPEND LD_ENV ${HDF5_TOOLS_DIR})
endif()

if (${WINDOWS}) 
  set (LD_ENV
    "${LD_ENV}\;${CMAKE_CURRENT_BINARY_DIR}/modules")
  set (LD_ENV
    "${LD_ENV}\;${CMAKE_CURRENT_BINARY_DIR}/modules/Release")
  set (LD_ENV "${LD_ENV}\;${CMAKE_CURRENT_BINARY_DIR}/Release")
  set (LD_ENV "${LD_ENV}\;${BOOST_LIBRARYDIR}")
else ()
  set (LD_ENV
    "${LD_ENV}:${CMAKE_CURRENT_BINARY_DIR}/modules")
  set (LD_ENV "${LD_ENV}:${BOOST_LIBRARYDIR}")
endif ()

if (${WINDOWS})
  string(REPLACE "/" "\\" LD_ENV "${LD_ENV}")
endif()

# Test macro -------------------------------------------------------- 
include (TestMacro)

set (TEST_CALL codeare -b ${PROJECT_SOURCE_DIR}/share/sense -c 3fold.xml)  
MP_TESTS ("csense" "${TEST_CALL}")

set (TEST_CALL codeare -b ${PROJECT_SOURCE_DIR}/share/grappa -c 4fold.xml)  
MP_TESTS ("grappa" "${TEST_CALL}")

set (TEST_CALL codeare -b ${PROJECT_SOURCE_DIR}/share/ktpoints -c 7t-8ch-8p.xml)  
MP_TESTS ("ktpoints" "${TEST_CALL}")

set (TEST_CALL codeare -b ${PROJECT_SOURCE_DIR}/share/compressedsensing -c phantom_denoising.xml) 
MP_TESTS ("denoise_cs" "${TEST_CALL}")

#if (${GSL_FOUND})
#   set (TEST_CALL codeare -b ${PROJECT_SOURCE_DIR}/share/gradients -c config_vdcs.xml) 
#   MP_TESTS ("vdspiral" "${TEST_CALL}")
#endif()

if (${NFFT3_FOUND})

  set (TEST_CALL codeare -b ${PROJECT_SOURCE_DIR}/share/nufft -c 16shot.xml)
  MP_TESTS ("nufft" "${TEST_CALL}")

  set (TEST_CALL codeare -b ${PROJECT_SOURCE_DIR}/share/cgsense -c 3fold.xml) 
  MP_TESTS ("cgsense" "${TEST_CALL}")

  set (TEST_CALL codeare -b ${PROJECT_SOURCE_DIR}/share/compressedsensing -c nufft_cs_lbfgs.xml) 
  MP_TESTS ("nufft_cs" "${TEST_CALL}")

  set (TEST_CALL codeare -b ${PROJECT_SOURCE_DIR}/share/compressedsensing -c cgsense_r3_cs_nlcg.xml) 
  MP_TESTS ("cgsense_cs" "${TEST_CALL}")

endif()

