cmake_minimum_required(VERSION 3.15)
project(WL
        VERSION 1.1.0
        DESCRIPTION "Wang-Landau algorithm"
        HOMEPAGE_URL "https://github.com/DavidAce/WL")

include(cmake/Deprecation.cmake) # Warn if user passes deprecated settings

# Set options this build
option(WL_INSTALL                         "Enable installation"                                          OFF)
option(WL_ENABLE_ASAN                     "Enable runtime address sanitizer -fsanitize=address"          OFF)
option(WL_CMAKE_DEBUG                     "Extra information during CMake configuration"                 OFF)
option(COMPILER_ENABLE_ASAN               "Enable runtime address sanitizer -fsanitize=address"          OFF)
option(COMPILER_ENABLE_USAN               "Enable undefined behavior sanitizer -fsanitize=undefined"     OFF)
option(CMAKE_POSITION_INDEPENDENT_CODE    "Use -fPIC when compiling shared libraries"                    ON)
option(THREADS_PREFER_PTHREAD_FLAG        "Prefer -pthread flag over -lpthread or similar"               ON)
option(CMAKE_INTERPROCEDURAL_OPTIMIZATION "Enable link-time-optimization"                                OFF)

# Make an "enum" for valid download methods
set(WL_PACKAGE_MANAGERS_VALID find cmake cpm conan)
set(WL_PACKAGE_MANAGER cmake CACHE STRING "Package manager for external dependencies")
set_property(CACHE WL_PACKAGE_MANAGER PROPERTY STRINGS ${WL_PACKAGE_MANAGERS_VALID})
if (NOT WL_PACKAGE_MANAGER IN_LIST WL_PACKAGE_MANAGERS_VALID)
    message(FATAL_ERROR "WL_PACKAGE_MANAGER must be one of ${WL_PACKAGE_MANAGERS_VALID}")
endif ()


################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git version numbers.###
################################################################
include(cmake/gitversion.cmake)


### Print extra info during CMake configure
include(cmake/PrintBuildInfo.cmake)


###  Add targets to collect common settings
add_library(wl-flags INTERFACE)
add_library(wl-deps INTERFACE)


###  Apply compiler flags
include(cmake/CompilerFlags.cmake)

### Find or install all the dependencies
include(cmake/SetupDependencies.cmake)

###########################################
###  Add all source files recursively   ###
###########################################
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC source)

# Link all the things
target_link_libraries(${PROJECT_NAME} PUBLIC wl-deps wl-flags)

# Print summary of CMake configuration
if (WL_CMAKE_DEBUG)
    # Print summary of CMake configuration
    include(cmake/PrintTargetInfo.cmake)
    print_and_write_project_summary(WL)
endif()


