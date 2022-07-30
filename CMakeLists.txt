cmake_minimum_required(VERSION 3.15)
project(WL
        VERSION 1.0.0
        DESCRIPTION "Wang-Landau algorithm"
        HOMEPAGE_URL "https://github.com/DavidAce/WL")

include(cmake/DeprecationWarnings.cmake) # Warn if user passes deprecated settings

# Set options this build
option(WL_INSTALL                       "Enable installation"                                          OFF)
option(WL_ENABLE_LTO                    "Enable link time optimization"                                OFF)
option(WL_ENABLE_ASAN                   "Enable runtime address sanitizer -fsanitize=address"          OFF)
option(WL_PREFIX_ADD_PKGNAME            "Install dependencies into CMAKE_INSTALL_PREFIX/<libname>"     OFF)
option(CMAKE_POSITION_INDEPENDENT_CODE  "Use -fPIC when compiling shared libraries"                    ON)


# Make an "enum" for valid download methods
set(WL_PACKAGE_MANAGERS_VALID find cmake conan)
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
add_library(flags INTERFACE)
add_library(deps INTERFACE)


###  Apply compiler flags
include(cmake/CompilerFlags.cmake)

### Find or install all the dependencies
include(cmake/SetupDependencies.cmake)

###########################################
###  Add all source files recursively   ###
###########################################
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")
file(GLOB_RECURSE EMC_SOURCES "EMC-Lib/source/*.cpp")
file(GLOB_RECURSE EMC_HEADERS "EMC-Lib/source/*.hpp")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS} ${EMC_SOURCES} ${EMC_HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC EMC-Lib)
target_include_directories(${PROJECT_NAME} PUBLIC source)
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_WHAT_YOU_USE TRUE)

# Link all the things
target_link_libraries(${PROJECT_NAME} PUBLIC deps flags)

# Print summary
include(cmake/PrintTargetInfo.cmake)
print_project_summary(${PROJECT_NAME})



###################
### INFORMATION ###
###################
message(STATUS "")
message(STATUS "=============================================================================")
message(STATUS "")
message(STATUS "      To clean CMake cache, run './build.sh -c'")
message(STATUS "      To launch this program,  './run.sh'")
message(STATUS "")
message(STATUS "=============================================================================")
message(STATUS "")



