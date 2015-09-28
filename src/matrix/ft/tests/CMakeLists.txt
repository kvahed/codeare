if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
endif()

include_directories(
        ${PROJECT_SOURCE_DIR}/src/core
        ${PROJECT_SOURCE_DIR}/src/matrix
        ${PROJECT_SOURCE_DIR}/src/matrix/simd
        ${PROJECT_SOURCE_DIR}/src/matrix/ft
        ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
        ${PROJECT_SOURCE_DIR}/src/matrix/io)

if (${WINDOWS})
  set (LD_ENV "PATH=%PATH%")
  list (APPEND LD_ENV ${HDF5_TOOLS_DIR})
  list (APPEND LD_ENV Release)
  list (APPEND LD_ENV ${FFTW3_INCLUDE_DIR})
  list (APPEND LD_ENV "${NFFT_INCLUDE_DIR}/../bin")
endif ()

if(${MATLAB_FOUND})
  if (${WINDOWS}) 
    set (LD_ENV "${LD_ENV}\;${MATLAB_ROOT}/bin/win64") 
  else ()
    set (LD_ENV "${LD_ENV}:${MATLAB_LIBRARY_DIRS}")
  endif ()
endif()

add_executable(t_dft t_dft.cpp)
target_link_libraries (t_dft ${FFTW3_LIBRARIES})

include (TestMacro)

# if(${MATLAB_FOUND})
#   add_executable(t_fftshift t_fftshift.cpp)
#   if (${MSVC})
#     target_link_libraries (t_fftshift libmx libmat ${FFTW3_LIBRARIES} hdf5 hdf5_cpp)
#   else()
#     target_link_libraries (t_fftshift mx mat ${FFTW3_LIBRARIES} ${HDF5_LIBRARIES})
#   endif()
#   set (TEST_CALL t_fftshift)  
#   MP_TESTS ("fftshift" "${TEST_CALL}")
# endif()

set (TEST_CALL t_dft)  
MP_TESTS ("dft" "${TEST_CALL}")

