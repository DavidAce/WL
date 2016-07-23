cmake_minimum_required(VERSION 3.5)
project(WL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
set(SOURCE_FILES main.cpp)
add_executable(WL ${SOURCE_FILES})