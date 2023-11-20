cmake_minimum_required(VERSION 3.16)
project(Lab1 VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

add_executable(${PROJECT_NAME} src/main.cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
