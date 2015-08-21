if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW) # CMake 3.0
endif()

include_directories (${PROJECT_SOURCE_DIR}/src/optimisation
  ${PROJECT_SOURCE_DIR}/src/matrix/io)

list (APPEND OPMTIMISATION_SOURCE Linear.hpp CGLS.hpp CGLS.cpp
  NonLinear.hpp NLCG.hpp NLCG.cpp SplitBregman.hpp SplitBregman.cpp
  LBFGS.hpp LBFGS.cpp lbfgs.h arithmetic_ansi.h lbfgs.h
  arithmetic_sse_double.h arithmetic_sse_float.h) 

SET_PROPERTY(SOURCE lbfgs.c PROPERTY COMPILE_FLAGS "-DHAVE_SSE2")
if (${HAVE_EMMINTRIN_H})
  SET_PROPERTY(SOURCE lbfgs.c PROPERTY COMPILE_FLAGS "-DHAVE_SSE2 -DHAVE_EMMINTRIN_H")
endif()

add_library (codeare-optimisation SHARED ${OPMTIMISATION_SOURCE})
target_link_libraries (codeare-optimisation  ${BLAS_LINKER_FLAGS}
  ${BLAS_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
list (APPEND INST_TARGETS codeare-optimisation)
install (TARGETS ${INST_TARGETS} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 
