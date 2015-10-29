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
  if(${Matlab_FOUND})
    list (APPEND COMLIBS  ${Matlab_MX_LIBRARIES} ${Matlab_MAT_LIBRARIES})
  endif()
else()
  set (COMLIBS ${HDF5_LIBRARIES})
  if(${Matlab_FOUND})
    list (APPEND COMLIBS)
  endif()
endif()

if (${WINDOWS})
  set (LD_ENV "PATH=%PATH%")
  list (APPEND LD_ENV ${HDF5_TOOLS_DIR})
  list (APPEND LD_ENV Release)
endif()

if(${Matlab_FOUND})
  if (${WINDOWS}) 
    set (LD_ENV "${LD_ENV}\;${Matlab_ROOT}/bin/win64") 
  else ()
    set (LD_ENV "${LD_ENV}:${Matlab_LIBRARY_DIRS}")
  endif ()
endif()


target_link_libraries (t_dwt ${COMLIBS})

set (TEST_CALL t_dwt)  
MP_TESTS ("dwt" "${TEST_CALL}")

