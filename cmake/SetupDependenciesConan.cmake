cmake_minimum_required(VERSION 3.15)

if(WL_PACKAGE_MANAGER MATCHES "conan")
    ##################################################################
    ### Install dependencies from conanfile.py                     ###
    ##################################################################
    unset(CONAN_COMMAND CACHE)
    find_program (CONAN_COMMAND conan PATHS ${WL_CONAN_HINTS} PATH_SUFFIXES ${WL_CONAN_PATH_SUFFIXES} REQUIRED)

    # Download cmake-conan integrator
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan/conan.cmake")
        message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
        file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/0.18.1/conan.cmake"
             "${CMAKE_BINARY_DIR}/conan/conan.cmake"
             EXPECTED_HASH MD5=81d5eab13a49f43527e35a90bfac6960
             TLS_VERIFY ON)
    endif()

    if(BUILD_SHARED_LIBS)
        list(APPEND CONAN_OPTIONS OPTIONS "*:shared=True")
    else()
        list(APPEND CONAN_OPTIONS OPTIONS "*:shared=False")
    endif()

    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        list(APPEND CONAN_OPTIONS OPTIONS "ceres-solver:use_glog=False")
    endif()

    # Copy the current compiler flags to conan
    string(TOUPPER "${CMAKE_BUILD_TYPE}" CONAN_BUILD_TYPE)
    set(CONAN_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CONAN_BUILD_TYPE}}")
    set(CONAN_CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CONAN_BUILD_TYPE}}")
    set(CONAN_LDFLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    conan_cmake_autodetect(CONAN_AUTODETECT)
    conan_cmake_install(
            CONAN_COMMAND ${CONAN_COMMAND}
            BUILD ${CONAN_BUILD} missing outdated cascade
            UPDATE
            GENERATOR CMakeDeps
            SETTINGS ${CONAN_AUTODETECT}
            INSTALL_FOLDER ${CMAKE_BINARY_DIR}/conan
            ENV CC=${CMAKE_C_COMPILER} # Fixes issue with CMake not detecting the right compiler when not building from scratch
            ENV CXX=${CMAKE_CXX_COMPILER} # Fixes issue with CMake not detecting the right compiler when not building from scratch
            ENV CXXFLAGS=${CONAN_CXXFLAGS}
            ENV CFLAGS=${CONAN_CFLAGS}
            ENV LDFLAGS=${CONAN_LDFLAGS}
            ${CONAN_OPTIONS}
            PATH_OR_REFERENCE ${CMAKE_SOURCE_DIR}
    )

    ##################################################################
    ### Find all the things!                                       ###
    ##################################################################
    list(PREPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR}/conan)
    list(PREPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/conan)
    list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
    # Use CONFIG to avoid MODULE mode. This is recommended for the cmake_find_package_multi generator
    find_package(Eigen3 3.4 REQUIRED CONFIG)
    find_package(spdlog 1.10.0 REQUIRED CONFIG)
    find_package(fmt 8.1.1 REQUIRED CONFIG)
    find_package(h5pp 1.10.0 REQUIRED CONFIG)
    find_package(CLI11 2.2.0 REQUIRED CONFIG)
    target_link_libraries(wl-deps INTERFACE Eigen3::Eigen spdlog::spdlog fmt::fmt h5pp::h5pp CLI11::CLI11)
endif()
