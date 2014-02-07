if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
endif()

include_directories (${PROJECT_SOURCE_DIR}/src/matrix/io ${PROJECT_SOURCE_DIR}/src/core ${PROJECT_SOURCE_DIR}/src/matrix/ft ${PROJECT_SOURCE_DIR}/src/matrix/arithmetic ${PROJECT_SOURCE_DIR}/src/matrix/dwt)

add_library (DummyRecon SHARED DummyRecon.hpp DummyRecon.cpp)
target_link_libraries (DummyRecon tinyxml core)

#add_library (CompressedSensing SHARED CompressedSensing.hpp CompressedSensing.cpp)
#target_link_libraries (CompressedSensing tinyxml core)
