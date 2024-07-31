enable_language(CXX)
cmake_minimum_required(VERSION 3.10)

project(BIT LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(/usr/local/include/matplotlibcpp)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
link_libraries(${PYTHON_LIBRARIES})

add_executable(main main.cpp BITStarforssl.cpp)
target_link_libraries(main ${PYTHON_LIBRARIES})