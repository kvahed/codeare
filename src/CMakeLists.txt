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

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.hpp" @ONLY) 
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
    message (STATUS ${MATLAB_LIBRARY_DIRS})
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

  add_custom_target (RRSModuleSK.cc 
    COMMAND ${OMNIORB4_IDL_COMPILER} -bcxx RRSModule.idl 
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src)

  add_custom_target (RRSModule.hh 
    COMMAND ${OMNIORB4_IDL_COMPILER} -bcxx RRSModule.idl 
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src)

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    list (APPEND OORB_DEFINES __WIN32__ __x86__ _WIN32_WINNT=0x0400)
    list (APPEND OMNIORB4_LIBRARIES ws2_32 mswsock advapi32)
  endif()

  ### rclient ###
  add_executable (rclient ${SOURCES} RemoteConnector.hpp
    RemoteConnector.cpp testclt.cpp RRSModule.hh RRSModuleSK.cc) 
  target_compile_definitions (rclient PRIVATE REMOTE ${OORB_DEFINES})
  target_link_libraries (rclient ${COMLIBS} ${OMNIORB4_LIBRARIES})

  ### codeared ###
  add_executable (codeared ${SOURCES} codeared.cpp RRSModule.hh
    RRSModuleSK.cc ReconServant.hpp ReconServant.cpp)
  target_compile_definitions (codeared PRIVATE REMOTE ${OORB_DEFINES})
  target_link_libraries (codeared mongoose ${COMLIBS}
    ${OMNIORB4_LIBRARIES})  

  install (TARGETS rclient codeared DESTINATION
    ${CMAKE_INSTALL_PREFIX}/bin) 

endif()

if (${LINUX})
  set (RT_ENV "LD_LIBRARY_PATH=$LD_LIBRARY_PATH")
elseif (${MACOSX})
  set (RT_ENV "DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH")
elseif (${WINDOWS})
  set (RT_ENV "PATH=%PATH%")
endif()

if (${MATLAB_FOUND})
  set (RT_ENV "${RT_ENV}\;${MATLAB_LIBRARY_DIRS}")
endif()

if (${ISMRMRD_FOUND})
  set (RT_ENV "${RT_ENV}:${ISMRMRD_INCLUDE_DIR}/../lib")
endif()

if (${WINDOWS})
  set (RT_ENV "${RT_ENV}\;${CMAKE_CURRENT_BINARY_DIR}/modules/Release")
else ()
  set (RT_ENV "${RT_ENV}:${CMAKE_INSTALL_PREFIX}/lib")
endif ()

message (STATUS ${RT_ENV})

add_test (NAME sense 
	WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/Release"
	COMMAND lclient -t SENSE -b ${PROJECT_SOURCE_DIR}/share/sense -c 3fold.xml)
set_tests_properties (sense PROPERTIES ENVIRONMENT ${RT_ENV})

if (${NFFT_FOUND})
  add_test (nufft lclient -t NuFFT -b ${PROJECT_SOURCE_DIR}/share/nufft/
    -c 16shot.xml)
  set_tests_properties (nufft PROPERTIES ENVIRONMENT ${RT_ENV}) 

  add_test (cgsense lclient -t CGSENSE -b
    ${PROJECT_SOURCE_DIR}/share/cgsense -c 3fold.xml)
  set_tests_properties (cgsense PROPERTIES ENVIRONMENT ${RT_ENV})
endif()

add_test (NAME ktpoints 
	WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/Release"
	COMMAND lclient -t KTPoints -b
  ${PROJECT_SOURCE_DIR}/share/ktpoints -c 7t-8ch-8p.xml)
set_tests_properties (ktpoints PROPERTIES ENVIRONMENT ${RT_ENV})

add_test (NAME compressedsensing
	WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/Release"
	COMMAND lclient -t CompressedSensing -b
	${PROJECT_SOURCE_DIR}/share/compressedsensing -c
	phantom_denoising.xml) 
set_tests_properties (compressedsensing PROPERTIES ENVIRONMENT ${RT_ENV})

add_test (NAME vdspiral
	WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/Release"
	COMMAND lclient -t VDSpiral -b
	${PROJECT_SOURCE_DIR}/share/gradients -c config_vdcs.xml)	
set_tests_properties (vdspiral PROPERTIES ENVIRONMENT ${RT_ENV}) 
