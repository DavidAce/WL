
if(WL_PACKAGE_MANAGER MATCHES "cpm")
    include(cmake/CPM.cmake)
    CPMAddPackage(
            NAME eigen
            GITLAB_REPOSITORY libeigen/eigen
            GIT_TAG 3.4.0
            OPTIONS
            "EIGEN_TEST_CXX11 OFF"
            "BUILD_TESTING OFF"
            "EIGEN_BUILD_DOC OFF"
            "EIGEN_BUILD_PKGCONFIG OFF"
    )
    target_link_libraries(wl-deps INTERFACE eigen)

    CPMAddPackage(NAME fmt
            GITHUB_REPOSITORY fmtlib/fmt
            GIT_TAG 8.1.1
            OPTIONS "FMT_INSTALL ON" )
    target_link_libraries(wl-deps INTERFACE fmt)


    CPMAddPackage(
            NAME spdlog
            GITHUB_REPOSITORY gabime/spdlog
            VERSION 1.10.0
            OPTIONS
            "SPDLOG_INSTALL ON"
            "SPDLOG_BUILD_SHARED ${BUILD_SHARED_LIBS}"
            "SPDLOG_FMT_EXTERNAL ON"
            "SPDLOG_FMT_EXTERNAL_HO OFF"
            "SPDLOG_ENABLE_PCH ON"
            "SPDLOG_BUILD_EXAMPLE OFF"
            "SPDLOG_BUILD_EXAMPLE_HO OFF"
            "SPDLOG_BUILD_TESTS OFF"
            "SPDLOG_BUILD_TESTS_HO OFF"
    )
    target_link_libraries(wl-deps INTERFACE spdlog)

    CPMAddPackage(
            NAME CLI11
            GITHUB_REPOSITORY CLIUtils/CLI11
            VERSION 2.2.0
            OPTIONS
            "CLI11_BUILD_EXAMPLES OFF"
            "CLI11_BUILD_TESTS OFF"
            "CLI11_BUILD_DOC OFF")
    target_link_libraries(wl-deps INTERFACE CLI11)


    CPMAddPackage(
            NAME h5pp
            GITHUB_REPOSITORY DavidAce/h5pp
            VERSION 1.10.0
            OPTIONS
            "H5PP_PACKAGE_MANAGER cpm"
            "H5PP_ENABLE_EIGEN3 ON"
            "H5PP_ENABLE_SPDLOG ON"
            "H5PP_ENABLE_FMT    ON"
    )
    find_library(ZLIB_LIBRARY
            NAMES
            z zlib zdll zlib1 zlibstatic # Release names
            zd zlibd zdlld zlibd1 zlib1d zlibstaticd # Debug names
            HINTS ${CMAKE_INSTALL_PREFIX} NO_DEFAULT_PATH )
    find_package(ZLIB REQUIRED)
    find_package(SZIP REQUIRED HINTS ${CMAKE_INSTALL_PREFIX} NO_DEFAULT_PATH)
    find_package(HDF5 1.13.1 REQUIRED CONFIG HINTS ${CMAKE_INSTALL_PREFIX} NO_DEFAULT_PATH)
    if(NOT TARGET HDF5::HDF5)
        if(BUILD_SHARED_LIBS)
            set(HDF5_TGTSFX shared)
        else()
            set(HDF5_TGTSFX static)
        endif()
        add_library(HDF5::HDF5 IMPORTED INTERFACE)
        target_link_libraries(HDF5::HDF5 INTERFACE hdf5::hdf5_hl-${HDF5_TGTSFX} hdf5::hdf5-${HDF5_TGTSFX})
    endif()
    target_link_libraries(wl-deps INTERFACE h5pp)

endif()