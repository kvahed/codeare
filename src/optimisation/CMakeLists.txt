include_directories (${PROJECT_SOURCE_DIR}/src/optimisation
  ${PROJECT_SOURCE_DIR}/src/matrix/io)

list (APPEND OPMTIMISATION_SOURCE NonLinear.hpp
  NLCG.hpp NLCG.cpp SplitBregman.hpp SplitBregman.cpp)

if (${LBFGS_FOUND})
  list (APPEND OPMTIMISATION_SOURCE LBFGS.hpp LBFGS.cpp)
endif()

add_library (optimisation ${OPMTIMISATION_SOURCE})
