cmake_minimum_required(VERSION 3.26)
project(ODESolver)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)
set (CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Add the path where you have your numpy headers installed.
include_directories(/opt/homebrew/lib/python3.11/site-packages/numpy/core/include)
find_library(PYTHON_LIB
             NAMES libpython3.11.dylib
             PATHS /opt/homebrew/Cellar/python@3.11/3.11.3/Frameworks/Python.framework/Versions/3.11/lib
             NO_DEFAULT_PATH)

find_package(Matplot++ REQUIRED)

# add_executable(main.out main.cpp odesolvers.cpp)
add_executable(example.out example.cpp ode_solver.cpp)
target_link_libraries(example.out ${PYTHON_LIB})
target_link_libraries(example.out Matplot++::matplot)
