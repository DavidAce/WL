cmake_minimum_required(VERSION 3.15)

# Append search paths for find_package and find_library calls
list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)

# Make sure find_library prefers static/shared library depending on BUILD_SHARED_LIBS
# This is important when finding dependencies such as zlib which provides both shared and static libraries.
# Note that we do not force this cache variable, so users can override it
if(BUILD_SHARED_LIBS)
    # This is order is the default
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX};${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE STRING "Prefer finding shared libraries")
else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX};${CMAKE_SHARED_LIBRARY_SUFFIX} CACHE STRING "Prefer finding static libraries")
endif()

# Transform CMAKE_INSTALL_PREFIX to full path
if(DEFINED CMAKE_INSTALL_PREFIX AND NOT IS_ABSOLUTE CMAKE_INSTALL_PREFIX)
    get_filename_component(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}
                           ABSOLUTE BASE_DIR ${CMAKE_BINARY_DIR} CACHE FORCE)
    message(DEBUG "Setting absolute path CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
endif()

# Setup build and install directories for dependencies
if(NOT WL_PKG_BUILD_DIR)
    set(WL_PKG_BUILD_DIR ${CMAKE_BINARY_DIR}/pkg-build)
endif()
if(NOT WL_PKG_INSTALL_DIR)
    if(DEFINED CMAKE_INSTALL_PREFIX)
        set(WL_PKG_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}) # Install to the same location as WL by default
    else()
        set(WL_PKG_INSTALL_DIR ${CMAKE_BINARY_DIR}/pkg-install)
    endif()
endif()
if(WL_PACKAGE_MANAGER MATCHES "cmake")
    set(PKG_INSTALL_DIR_DEFAULT ${WL_PKG_INSTALL_DIR} CACHE INTERNAL "Install directory for automatically installed dependency packages")
    set(PKG_BUILD_DIR_DEFAULT ${WL_PKG_BUILD_DIR} CACHE INTERNAL "Build directory for automatically installed dependency packages")
endif()

set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${PKG_INSTALL_DIR_DEFAULT};${CMAKE_INSTALL_PREFIX}")
list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" CACHE INTERNAL "Paths for find_package lookup")

if(WL_PACKAGE_MANAGER MATCHES "conan")
    # Paths to search for conan installation.
    list(APPEND WL_CONAN_HINTS
         ${CONAN_PREFIX}
         $ENV{CONAN_PREFIX}
         ${CONDA_PREFIX}
         $ENV{CONDA_PREFIX}
         $ENV{HOME}/anaconda3
         $ENV{HOME}/anaconda
         $ENV{HOME}/miniconda3
         $ENV{HOME}/miniconda
         $ENV{HOME}/.conda
         )
    list(APPEND WL_CONAN_PATH_SUFFIXES bin envs/dmrg/bin)
    list(REMOVE_DUPLICATES WL_CONAN_HINTS)
    list(REMOVE_DUPLICATES WL_CONAN_PATH_SUFFIXES)
    mark_as_advanced(WL_CONAN_HINTS)
    mark_as_advanced(WL_CONAN_PATH_SUFFIXES)
endif()
