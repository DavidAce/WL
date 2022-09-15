
if(WL_PACKAGE_MANAGER MATCHES "cmake")
    include(cmake/InstallPackage.cmake)
    list(PREPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR}/${WL_PKG_INSTALL_DIR})
    list(PREPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/${WL_PKG_INSTALL_DIR})
    list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)

    # Install all the depeendencies
    install_package(fmt VERSION 8.1.1)
    install_package(spdlog VERSION 1.10.0
            DEPENDS fmt::fmt
            CMAKE_ARGS
            -DSPDLOG_FMT_EXTERNAL:BOOL=ON
            -Dfmt_ROOT:PATH=${WL_PKG_INSTALL_DIR})
    install_package(Eigen3 VERSION 3.4.0 TARGET_NAME Eigen3::Eigen)
    install_package(cli11 VERSION 2.1.1 TARGET_NAME CLI11::CLI11 FIND_NAME CLI11)
    install_package(h5pp VERSION 1.10.1
                    CMAKE_ARGS
                    -DH5PP_PACKAGE_MANAGER=${WL_PACKAGE_MANAGER}
                    -DCMAKE_PREFIX_PATH:PATH=${WL_PKG_INSTALL_DIR} )

    # Link to dependencies
    target_link_libraries(wl-deps INTERFACE fmt::fmt spdlog::spdlog Eigen3::Eigen CLI11::CLI11 h5pp::h5pp)

endif()
