cmake_minimum_required(VERSION 3.15)
project(external-h5pp)
if(INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif()
include(GNUInstallDirs)

option(BUILD_SHARED_LIBS               "Builds shared libraries"                                     OFF)
option(H5PP_BUILD_EXAMPLES             "Builds examples"                                             OFF)
option(H5PP_ENABLE_TESTS               "Enable testing"                                              OFF)
option(H5PP_PREFIX_ADD_PKGNAME         "Install dependencies into CMAKE_INSTALL_PREFIX/<libname>"    OFF)
option(H5PP_ENABLE_EIGEN3              "Enables Eigen3 linear algebra library"                       ON)
option(H5PP_ENABLE_FMT                 "Enables the {fmt} formatting library"                        ON)
option(H5PP_ENABLE_SPDLOG              "Enables Spdlog for logging h5pp internal info to stdout"     ON)
option(H5PP_ENABLE_MPI                 "Enables use of MPI (work in progress)"                       OFF)
option(H5PP_ENABLE_ASAN                "Enable runtime address sanitizer -fsanitize=address"         OFF)
option(CMAKE_POSITION_INDEPENDENT_CODE "Use -fPIC when compiling shared libraries"                   ON)


include(ExternalProject)
ExternalProject_Add(external_h5pp
        URL https://github.com/DavidAce/h5pp/archive/refs/tags/v1.10.0.tar.gz
        URL_MD5 e1d21f9516df086184f794f9bce3bf4a

        BUILD_ALWAYS TRUE
        PREFIX          ${CMAKE_BINARY_DIR}
        INSTALL_DIR     ${CMAKE_INSTALL_PREFIX}
        #Give Ninja access to the terminal.
        USES_TERMINAL_DOWNLOAD TRUE
        USES_TERMINAL_UPDATE TRUE
        USES_TERMINAL_CONFIGURE TRUE
        USES_TERMINAL_BUILD TRUE
        USES_TERMINAL_INSTALL TRUE
        USES_TERMINAL_TEST TRUE

        CMAKE_ARGS
        ${INIT_CACHE_ARG}
        -DCMAKE_POLICY_DEFAULT_CMP0074=NEW
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_INSTALL_RPATH:PATH=<INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR}
        -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
        -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON
        -DCMAKE_INSTALL_MESSAGE=LAZY
        -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=${CMAKE_POSITION_INDEPENDENT_CODE}

        # h5pp flags
        -DH5PP_BUILD_EXAMPLES:BOOL=${H5PP_BUILD_EXAMPLES}
        -DH5PP_ENABLE_TESTS:BOOL=${H5PP_ENABLE_TESTS}
        -DH5PP_PREFIX_ADD_PKGNAME:BOOL=${H5PP_PREFIX_ADD_PKGNAME}
        -DH5PP_ENABLE_EIGEN3:BOOL=${H5PP_ENABLE_EIGEN3}
        -DH5PP_ENABLE_FMT:BOOL=${H5PP_ENABLE_FMT}
        -DH5PP_ENABLE_SPDLOG:BOOL=${H5PP_ENABLE_SPDLOG}
        -DH5PP_ENABLE_MPI:BOOL=${H5PP_ENABLE_MPI}
        -DH5PP_ENABLE_ASAN:BOOL=${H5PP_ENABLE_ASAN}
        -DH5PP_PACKAGE_MANAGER:STRING=${H5PP_PACKAGE_MANAGER}
        -Dfmt_ROOT:PATH=${fmt_ROOT}
        -Dspdlog_ROOT:PATH=${spdlog_ROOT}
        -DEigen3_ROOT:PATH=${Eigen3_ROOT}
        )