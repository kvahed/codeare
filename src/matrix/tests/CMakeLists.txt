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

set (COMLIBS ${FFTW3_LIBRARIES})
if (${MSVC})
  list (APPEND COMLIBS hdf5 hdf5_cpp)
  if(${MATLAB_FOUND})
    list (APPEND COMLIBS libmx libmat)
  endif()
else()
  list (APPEND COMLIBS ${HDF5_LIBRARIES})
  if(${MATLAB_FOUND})
    list (APPEND COMLIBS mx mat)
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

add_executable(t_emul t_emul.cpp)
add_test(emul t_emul)

add_executable(t_ediv t_ediv.cpp)
add_test(ediv t_ediv)

add_executable(t_eadd t_eadd.cpp)
add_test(eadd t_eadd)

add_executable(t_esub t_esub.cpp)
add_test(esub t_esub)

add_executable(t_cast t_cast.cpp)
add_test(cast t_cast)

add_executable(t_sort t_sort.cpp)
add_test(sort t_sort)

add_executable(t_diff t_diff.cpp)
add_test(diff t_diff)

add_executable(t_flip t_flip.cpp)
add_test(flip t_flip)

add_executable(t_friends t_friends.cpp)
add_test(friends t_friends)

if (${HAS_CXX11_CONDITIONAL}) 
add_executable(t_colon t_colon.cpp)
add_test(colon t_colon)
add_executable(t_view t_view.cpp)
add_test(view t_view)
endif()


