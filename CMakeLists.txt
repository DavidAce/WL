cmake_minimum_required(VERSION 2.8)
project(WL)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
#find_package(Eigen3 REQUIRED)
#include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(eigen3)

# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

##################################
# Check for c++14 compiler support
##################################
include(CheckCXXCompilerFlag)
set(ENABLE_CXXFLAGS_TO_CHECK
        -std=c++14
        -std=c++1y
        -std=gnu++1y
        -std=gnu++14
        -std=c++11
        -std=gnu++11
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


###################################
#Change flags depending on host
###################################
site_name(this_host)
if ("${this_host}" MATCHES "triolith")
    message("Host: ${this_host} (Intel) ")
    set(CMAKE_CXX_COMPILER mpiicpc)
    set(CMAKE_CXX_FLAGS_DEBUG "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -g -funroll-loops -DNDEBUG -cxxlib=/software/apps/gcc/5.3.0/build01/ -gxx-name=g++ -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -g -O3 -cxxlib=/software/apps/gcc/5.3.0/build01/ -gxx-name=g++ -ip -xavx -funroll-loops -DNDEBUG")
elseif("${this_host}" MATCHES "beskow")
    message("Host: ${this_host} (Intel) ")
    set(CMAKE_CXX_COMPILER CC)
    set(CMAKE_CXX_FLAGS_DEBUG "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -g -shared -ansi-alias -no-prec-div -ipo -xCORE-AVX2 -DEIGEN_NO_DEBUG -DNDEBUG  -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -g -O3 -shared -ansi-alias -no-prec-div -ipo -xCORE-AVX2 -DEIGEN_NO_DEBUG -DNDEBUG")
else()
    message("Host: ${this_host} (GNU)")
    set(CMAKE_CXX_COMPILER mpicxx)
    set(CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -Wno-deprecated-declarations")
    set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -g -O3 -funroll-loops -march=native -frename-registers -Wno-deprecated-declarations -DNDEBUG")
endif()
message("Compiling with: ${CMAKE_CXX_COMPILER_ID}")

#target_link_libraries(WL shared)


file(GLOB_RECURSE WL_SOURCES "source/*.cpp")
file(GLOB_RECURSE WL_HEADERS "source/*.h")

file(GLOB_RECURSE EMC_SOURCES "EMC-Lib/source/*.cpp")
file(GLOB_RECURSE EMC_HEADERS "EMC-Lib/source/*.hpp")

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

