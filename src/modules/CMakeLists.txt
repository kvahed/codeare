include_directories(${PROJECT_SOURCE_DIR}/src/matrix/io)

add_library (DummyRecon SHARED DummyRecon.hpp DummyRecon.cpp)
target_link_libraries (DummyRecon tinyxml)

#add_library (CompressedSensing SHARED ../Workspace.cpp CompressedSensing.hpp CompressedSensing.cpp)
#target_link_libraries (CompressedSensing tinyxml)
