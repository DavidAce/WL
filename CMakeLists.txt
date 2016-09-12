set(CMAKE_C_COMPILER "mpiicc")
set(CMAKE_CXX_COMPILER "mpiicpc")
cmake_minimum_required(VERSION 2.8)
project(WL)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
#find_package(Eigen3 REQUIRED)
#include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(eigen3)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
link_libraries(-lstdc++fs)

# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

##################################
# Check for c++14 compiler support
##################################
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++1y" COMPILER_SUPPORTS_CXX1Y)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
elseif(COMPILER_SUPPORTS_CXX1Y)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
else()
    message("Compiler ${CMAKE_CXX_COMPILER} may not have C++14 support.")
endif()
SET(CMAKE_CXX_STANDARD 14)

###################################
#Change flags depending on compiler
###################################
message("Compiling with: ${CMAKE_CXX_COMPILER_ID}")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    SET(CMAKE_CXX_COMPILER mpicxx)
    set(CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-deprecated-declarations")
    set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops -march=native -frename-registers -Wno-deprecated-declarations -DNDEBUG")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    SET(CMAKE_CXX_COMPILER mpiicpc)
    set(CMAKE_CXX_FLAGS_DEBUG "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O3 -ip -xavx -funroll-loops -DNDEBUG")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # using Visual Studio C++
endif()


file(GLOB_RECURSE WL_SOURCES "source/*.cpp")
file(GLOB_RECURSE WL_HEADERS "source/*.h")

file(GLOB_RECURSE EMC_SOURCES "EMC-Lib/*.cpp")
file(GLOB_RECURSE EMC_HEADERS "EMC-Lib/*.hpp")

set (WL_INCLUDE_DIRS "")
foreach (_headerFile ${WL_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND WL_INCLUDE_DIRS ${_dir})
endforeach()
foreach (_headerFile ${EMC_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND WL_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES WL_INCLUDE_DIRS)

include_directories(${WL_INCLUDE_DIRS})
add_executable (WL main.cpp ${WL_SOURCES} ${EMC_SOURCES})

