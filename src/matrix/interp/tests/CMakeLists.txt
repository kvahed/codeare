include_directories (
  ${PROJECT_SOURCE_DIR}/src/core
  ${PROJECT_SOURCE_DIR}/src/matrix
  ${PROJECT_SOURCE_DIR}/src/matrix/simd
  ${PROJECT_SOURCE_DIR}/src/matrix/interp
  ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic
  ${PROJECT_SOURCE_DIR}/src/matrix/io
  )

add_executable (t_medfilt2 t_medfilt2.cpp)
add_test (medfilt2 t_medfilt2)

if (${ITK_FOUND})
#  add_executable (t_resample t_resample.cpp)
#  set_target_properties (t_resample 
#    PROPERTIES COMPILE_FLAGS "-DHAVE_INSIGHT -fpermissive") 
#  target_link_libraries (t_resample ${ITK_LIBRARIES})
#  add_test (resample t_resample)
  #set_tests_properties (resample PROPERTIES ENVIRONMENT ${LD_ENV})
endif()

if (${GSL_FOUND})
  add_executable (t_ppval t_ppval.cpp)
  target_link_libraries (t_ppval ${GSL_LIBRARIES})
  add_test (ppval t_ppval)
  add_executable (t_interp1 t_interp1.cpp)
  target_link_libraries (t_interp1 ${GSL_LIBRARIES})
  add_test (interp1 t_interp1)
endif()

