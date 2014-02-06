link_directories(${PROJECT_BINARY_DIR}/src/tinyxml)
add_library (DummyRecon SHARED ../Workspace.hpp ../Workspace.hpp DummyRecon.hpp DummyRecon.cpp)
target_link_libraries (tinyxml)
