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

# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

##################################
# Check for c++14 compiler support
##################################
include(CheckCXXCompilerFlag)
set(ENABLE_CXXFLAGS_TO_CHECK
        -std=gnu++14
        -std=c++14
        -std=gnu++1y
        -std=c++1y
        -std=gnu++11
        -std=c++11
        -std=gnu++0x
        -std=c++0x)

foreach(flag ${ENABLE_CXXFLAGS_TO_CHECK})
    string(REPLACE "-std=" "_" flag_var ${flag})
    string(REPLACE "+" "x" flag_var ${flag_var})
    check_cxx_compiler_flag("${flag}" COMPILER_HAS_CXX_FLAG${flag_var})
    if(COMPILER_HAS_CXX_FLAG${flag_var})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
        break()
    endif()
endforeach()

#
#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
#CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#CHECK_FLAG_SUPPORTED(-std=c++14)
##CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
##CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#if(COMPILER_SUPPORTS_CXX14)
#    set(CMAKE_CXX_STANDARD 14)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++14")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
#elseif(COMPILER_SUPPORTS_CXX11)
#    set(CMAKE_CXX_STANDARD 11)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#else()
#    set(CMAKE_CXX_STANDARD 14)
#    message("WARNING: Compiler ${CMAKE_CXX_COMPILER} may not have C++11/C++14 support.")
#endif()



###################################
#Change flags depending on host
###################################
site_name(this_host)
if ("${this_host}" MATCHES "triolith")
    message("Host: ${this_host} (Intel) ")
    set(CMAKE_CXX_COMPILER mpiicpc)
    set(CMAKE_CXX_FLAGS_DEBUG "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -cxxlib=/software/apps/gcc/6.1.0/build01/ -gxx-name=g++ -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O3 -cxxlib=/software/apps/gcc/6.1.0/build01/ -gxx-name=g++ -ip -xavx -funroll-loops -DNDEBUG")
else()
    message("Host: ${this_host} (GNU)")
    set(CMAKE_CXX_COMPILER mpicxx)
    set(CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-deprecated-declarations")
    set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops -march=native -frename-registers -Wno-deprecated-declarations -DNDEBUG")
endif()
message("Compiling with: ${CMAKE_CXX_COMPILER_ID}")



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

