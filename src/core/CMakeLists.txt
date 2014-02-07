include_directories (${PROJECT_SOURCE_DIR}/src/core ${PROJECT_SOURCE_DIR}/src/matrix/io)
add_library (core Params.hpp Queue.hpp Queue.cpp ReconContext.hpp ReconContext.cpp Toolbox.hpp Toolbox.cpp Workspace.hpp Workspace.cpp)
