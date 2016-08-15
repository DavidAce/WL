cmake_minimum_required(VERSION 3.5)
project(WL)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -std=c++11 -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -std=c++11 -funroll-loops -march=native -frename-registers -Wno-deprecated-declarations -DNDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -funroll-loops -march=native -frename-registers -Wno-deprecated-declarations -DNDEBUG")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_CXX_STANDARD 11)


# Output directory
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")
#set(CMAKE_CURRENT_BINARY_DIR "build")
#SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE})

SET(SOURCE_FILES main.cpp source/class_worker.cpp source/class_worker.h source/class_master.cpp source/class_master.h source/class_model.cpp source/class_model.h source/constants.h source/randomNumbers.h source/randomNumbers.cpp source/wanglandau.cpp source/wanglandau.cpp source/wanglandau.h source/counters_timers.h source/class_data.cpp source/class_data.h source/math_algorithms.h source/math_algorithms.cpp source/MPI_algorithms.cpp source/MPI_algorithms.h source/class_profiling.cpp source/class_profiling.h)
add_executable(WL ${SOURCE_FILES})
