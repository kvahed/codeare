include_directories (${PROJECT_SOURCE_DIR}/src/core
${PROJECT_SOURCE_DIR}/src/matrix/io)
list (APPEND CORE_SOURCE Params.hpp Queue.hpp Queue.cpp
ReconContext.hpp ReconContext.cpp Toolbox.hpp Toolbox.cpp
Workspace.hpp Workspace.cpp)  
add_library (core ${CORE_SOURCE})
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  target_link_libraries (core dl)
endif()
