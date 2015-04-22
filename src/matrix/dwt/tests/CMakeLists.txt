include_directories(
        ${PROJECT_SOURCE_DIR}/src/core
        ${PROJECT_SOURCE_DIR}/src/matrix
        ${PROJECT_SOURCE_DIR}/src/matrix/simd
        ${PROJECT_SOURCE_DIR}/src/matrix/dwt
        ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
        ${PROJECT_SOURCE_DIR}/src/matrix/io)


add_executable (t_dwt t_dwt.cpp)

if (${MSVC})
  set (COMLIBS hdf5 hdf5_cpp)
  if(${MATLAB_FOUND})
    list (APPEND COMLIBS libmx libmat)
  endif()
else()
  set (COMLIBS ${HDF5_LIBRARIES})
  if(${MATLAB_FOUND})
<<<<<<< HEAD
    list (APPEND COMLIBS  mx mat)
=======
<<<<<<< HEAD
    list (APPEND COMLIBS ${MATLAB_LIBRARY_DIRS})
=======
    list (APPEND COMLIBS  mx mat)
>>>>>>> b450193f5e531dabcc0a09ff5224acc1ca2fbbfc
>>>>>>> c264ad4da08f67481e56e8e9e9cae9c4638a29fa
  endif()
endif()

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


target_link_libraries (t_dwt ${COMLIBS})

set (TEST_CALL t_dwt)  
MP_TESTS ("dwt" "${TEST_CALL}")

