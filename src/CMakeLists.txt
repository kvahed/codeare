include_directories(
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/src/matrix
        ${PROJECT_SOURCE_DIR}/src/matrix/simd
        ${PROJECT_SOURCE_DIR}/src/tinyxml
        ${PROJECT_SOURCE_DIR}/src/misc)
add_subdirectory(tinyxml)
add_subdirectory(mongoose)
add_subdirectory(modules)

#add_executable(lclient LocalConnector.hpp LocalConnector.cpp testclt.cpp)
