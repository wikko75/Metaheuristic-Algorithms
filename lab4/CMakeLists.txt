cmake_minimum_required(VERSION 3.16)

project(lab4)

find_package(fmt REQUIRED)

file(GLOB UTILS_SRC "../Utils/src/*.cpp")
file(GLOB GRAPH_SRC "../Graph/src/*.cpp")
file(GLOB MAIN_SRC "src/*.cpp")

set(UTILS_INC "../Utils/inc")
set(GRAPH_INC "../Graph/inc")
set(MAIN_INC "inc")

add_executable(${PROJECT_NAME} ${GRAPH_SRC} ${UTILS_SRC} ${MAIN_SRC} )

target_include_directories(${PROJECT_NAME} PRIVATE ${GRAPH_INC} ${UTILS_INC} ${MAIN_INC})

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -O3)

target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt)