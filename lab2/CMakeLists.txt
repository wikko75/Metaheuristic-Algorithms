cmake_minimum_required(VERSION 3.16)
project(lab2 VERSION 0.1.0 LANGUAGES CXX)

include(CTest)
enable_testing()

file(GLOB GRAPH_SRC "../Graph/src/*.cpp")
file(GLOB UTILS_SRC "../Utils/src/*.cpp")
file(GLOB MAIN_SRC "src/*.cpp")

add_executable(${PROJECT_NAME} ${UTILS_SRC} ${GRAPH_SRC} ${MAIN_SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
target_include_directories(${PROJECT_NAME} PRIVATE ../Utils/inc ../Graph/inc inc )