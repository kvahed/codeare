include_directories (${PROJECT_SOURCE_DIR}/src/optimisation
  ${PROJECT_SOURCE_DIR}/src/matrix/io)

list (APPEND OPMTIMISATION_SOURCE NonLinear.hpp
  NLCG.hpp NLCG.cpp SplitBregman.hpp SplitBregman.cpp)

if (${LBFGS_FOUND})
  list (APPEND OPMTIMISATION_SOURCE LBFGS.hpp LBFGS.cpp)
endif()

add_library (codeare-optimisation SHARED ${OPMTIMISATION_SOURCE})
if (${LBFGS_FOUND})
  target_link_libraries (codeare-optimisation ${LBFGS_LIBRARIES})  
endif()

list (APPEND INST_TARGETS codeare-optimisation)
install (TARGETS ${INST_TARGETS} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 
