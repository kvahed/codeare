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
    list (APPEND COMLIBS ${Matlab_MX_LIBRARIES} ${Matlab_MAT_LIBRARIES})
  endif()
else()
  list (APPEND COMLIBS ${HDF5_LIBRARIES})
  if(${MATLAB_FOUND})
    list (APPEND COMLIBS ${Matlab_MX_LIBRARIES} ${Matlab_MAT_LIBRARIES})
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

add_executable(t_econj t_econj.cpp)
add_test(conj t_econj)

add_executable(t_cov t_cov.cpp)
add_test(cov t_cov)
target_link_libraries (t_cov ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})

add_executable(t_emul t_emul.cpp)
add_test(emul t_emul)

add_executable(t_ediv t_ediv.cpp)
add_test(ediv t_ediv)

add_executable(t_eadd t_eadd.cpp)
add_test(eadd t_eadd)

add_executable(t_esub t_esub.cpp)
add_test(esub t_esub)

add_executable(t_sum t_sum.cpp)
add_test(sum t_sum)

add_executable(t_mean t_mean.cpp)
add_test(mean t_mean)

add_executable(t_max t_max.cpp)
add_test(max t_max)

add_executable(t_min t_min.cpp)
add_test(min t_min)
#add_executable(t_speed t_speed.cpp)
#add_test(speed t_speed)
#target_link_libraries (t_speed ${Boost_CHRONO_LIBRARY}
#  ${Boost_SYSTEM_LIBRARY} ${Boost_CHRONO_LIBRARY}
#  ${Boost_TIMER_LIBRARY})

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

add_executable(t_squeeze t_squeeze.cpp)
add_test(squeeze t_squeeze)

add_executable(t_issymmetric t_issymmetric.cpp)
add_test(issymmetric t_issymmetric)
target_link_libraries(t_issymmetric ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(t_ishermitian t_ishermitian.cpp)
add_test(ishermitian t_ishermitian)
target_link_libraries(t_ishermitian ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(t_issame t_issame.cpp)
add_test(issame t_issame)
target_link_libraries(t_issame)

if (${HAS_CXX11_CONDITIONAL}) 
add_executable(t_colon t_colon.cpp)
add_test(colon t_colon)
add_executable(t_view t_view.cpp)
add_test(view t_view)
endif()


