enable_language(CXX)
cmake_minimum_required(VERSION 3.10)

project(BIT LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

include_directories(/usr/local/include/matplotlibcpp)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
link_libraries(${PYTHON_LIBRARIES})

add_executable(Testplt Testplt.cpp)
target_link_libraries(Testplt ${PYTHON_LIBRARIES})
