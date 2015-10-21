include_directories (
  ${PROJECT_SOURCE_DIR}/src/core
  ${PROJECT_SOURCE_DIR}/src/matrix
  ${PROJECT_SOURCE_DIR}/src/matrix/curves
  ${PROJECT_SOURCE_DIR}/src/matrix/arithmetics
  )

add_executable (t_lmaxima t_lmaxima.cpp)
add_test (lmaxima t_lmaxima)
