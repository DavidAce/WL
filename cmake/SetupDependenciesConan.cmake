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
        file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/release/0.17/conan.cmake"
                "${CMAKE_BINARY_DIR}/conan/conan.cmake"
                EXPECTED_HASH MD5=52a255a933397fdce3d0937f9c737e98
                TLS_VERIFY ON)
    endif()

    include(${CMAKE_BINARY_DIR}/conan/conan.cmake)
    conan_cmake_autodetect(CONAN_AUTODETECT)
    conan_cmake_install(
            CONAN_COMMAND ${CONAN_COMMAND}
            BUILD missing outdated cascade
            GENERATOR CMakeDeps
            SETTINGS ${CONAN_AUTODETECT}
            INSTALL_FOLDER ${CMAKE_BINARY_DIR}/conan
            PATH_OR_REFERENCE ${CMAKE_SOURCE_DIR}
    )


    ##################################################################
    ### Find all the things!                                       ###
    ##################################################################
    if(NOT CONAN_CMAKE_SILENT_OUTPUT)
        set(CONAN_CMAKE_SILENT_OUTPUT OFF) # Default is off
    endif()
    list(PREPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/conan)
    list(PREPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR}/conan)
    # Use CONFIG to avoid MODULE mode. This is recommended for the cmake_find_package_multi generator
    find_package(Eigen3 3.4 REQUIRED CONFIG)
    find_package(spdlog 1.10.0 REQUIRED CONFIG)
    find_package(fmt 8.1.1 REQUIRED CONFIG)
    find_package(h5pp 1.10.0 REQUIRED CONFIG)
    find_package(CLI11::CLI11 REQUIRED CONFIG)
    target_link_libraries(wl-deps INTERFACE Eigen3::Eigen spdlog::spdlog fmt::fmt h5pp::h5pp CLI11::CLI11)
endif()
