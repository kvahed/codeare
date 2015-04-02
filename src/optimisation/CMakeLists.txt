include_directories (${PROJECT_SOURCE_DIR}/src/optimisation
  ${PROJECT_SOURCE_DIR}/src/matrix/io)

list (APPEND OPMTIMISATION_SOURCE NonLinear.hpp
  NLCG.hpp NLCG.cpp SplitBregman.hpp SplitBregman.cpp LBFGS.hpp
  LBFGS.cpp lbfgs.h arithmetic_ansi.h arithmetic_sse_double.h
  arithmetic_sse_float.h lbfgs.h lbfgs.c)

SET_PROPERTY(SOURCE foo.cpp PROPERTY COMPILE_FLAGS HAVE_SSE2)

add_library (codeare-optimisation SHARED ${OPMTIMISATION_SOURCE})
list (APPEND INST_TARGETS codeare-optimisation)
install (TARGETS ${INST_TARGETS} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 
