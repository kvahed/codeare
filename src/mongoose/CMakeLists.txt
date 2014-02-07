include_directories (${PROJECT_SOURCE_DIR}/src/core)
set (HEADER_FILES ${PROJECT_SOURCE_DIR}/src/core/Workspace.hpp)
add_library (mongoose mongoose.c mongoose.h MongooseService.cpp MongooseService.hpp ${HEADER_FILES})
