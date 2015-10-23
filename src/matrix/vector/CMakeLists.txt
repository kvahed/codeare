include_directories(
        ${PROJECT_SOURCE_DIR}/src/matrix
)

add_executable(t_find t_find.cpp)
add_test(find t_find)
