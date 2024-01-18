cmake_minimum_required(VERSION 3.16)

project(Lab4)

find_package(fmt REQUIRED)

add_executable(${PROJECT_NAME} src/main.cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt)