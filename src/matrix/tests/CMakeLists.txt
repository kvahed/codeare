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

add_executable(t_colon t_colon.cpp)
add_test(colon t_colon)

add_executable(t_sort t_sort.cpp)
add_test(sort t_sort)

set (COMLIBS ${FFTW3_LIBRARIES})
if (${MSVC})
  list (APPEND COMLIBS hdf5 hdf5_cpp)
  if(${MATLAB_FOUND})
    list (APPEND COMLIBS libmx libmat)
  endif()
else()
  list (APPEND COMLIBS ${HDF5_LIBRARIES})
  if(${MATLAB_FOUND})
<<<<<<< HEAD
    list (APPEND COMLIBS mx mat)
=======
<<<<<<< HEAD
    list (APPEND COMLIBS ${MATLAB_LIBRARY_DIRS})
=======
    list (APPEND COMLIBS mx mat)
>>>>>>> b450193f5e531dabcc0a09ff5224acc1ca2fbbfc
>>>>>>> c264ad4da08f67481e56e8e9e9cae9c4638a29fa
  endif()
endif()

add_executable(t_tvop t_tvop.cpp)
target_link_libraries (t_tvop ${COMLIBS})

include (TestMacro)

if (${WINDOWS})
  set (LD_ENV "PATH=%PATH%")
  list (APPEND LD_ENV ${HDF5_TOOLS_DIR})
  list (APPEND LD_ENV Release)
endif()

if(${MATLAB_FOUND})
  if (${WINDOWS}) 
    set (LD_ENV "${LD_ENV}\;${MATLAB_ROOT}/bin/win64") 
  else ()
    set (LD_ENV "${LD_ENV}:${MATLAB_LIBRARY_DIRS}")
  endif ()
endif()

set (TEST_CALL t_tvop)  
MP_TESTS ("tvop" "${TEST_CALL}")


