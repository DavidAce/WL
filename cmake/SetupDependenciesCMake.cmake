cmake_minimum_required(VERSION 3.15)

macro(unset_NOT_FOUND pkg)
    if(NOT ${pkg}_FOUND)
        unset(${pkg}_FOUND)
        unset(${pkg}_FOUND CACHE)
    endif()
endmacro()


if(WL_PACKAGE_MANAGER MATCHES "cmake")
    unset_NOT_FOUND(fmt)
    unset_NOT_FOUND(spdlog)
    unset_NOT_FOUND(Eigen3)

    include(cmake/InstallPackage.cmake)
    if(WL_PREFIX_ADD_PKGNAME)
        set(INSTALL_PREFIX_PKGNAME INSTALL_PREFIX_PKGNAME)
    endif()


    # Install all the depeendencies
    install_package(fmt VERSION 8.1.1 ${INSTALL_PREFIX_PKGNAME})
    install_package(spdlog VERSION 1.10.0 ${INSTALL_PREFIX_PKGNAME}
            DEPENDS fmt::fmt
            CMAKE_ARGS
            -DSPDLOG_FMT_EXTERNAL:BOOL=ON
            -Dfmt_ROOT:PATH=${WL_DEPS_INSTALL_DIR})
    install_package(Eigen3 VERSION 3.4.0 TARGET_NAME Eigen3::Eigen ${INSTALL_PREFIX_PKGNAME})

    # Link to dependencies
    target_link_libraries(deps INTERFACE ${PKG_fmt_TARGET} ${PKG_spdlog_TARGET} ${PKG_Eigen3_TARGET})

endif()
