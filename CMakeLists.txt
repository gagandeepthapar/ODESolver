cmake_minimum_required(VERSION 3.26)
project(DOUBLE)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)
set (CMAKE_EXPORT_COMPILE_COMMANDS 1)


add_executable(main main.cpp odesolvers.cpp)
add_executable(test test.cpp odesolvers.cpp)

