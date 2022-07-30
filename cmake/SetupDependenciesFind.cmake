if(WL_PACKAGE_MANAGER MATCHES "find")
    if(WL_PACKAGE_MANAGER STREQUAL "find")
        set(REQUIRED REQUIRED)
    endif()

    if(NOT "${CMAKE_BINARY_DIR}" IN_LIST CMAKE_MODULE_PATH)
        list(PREPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
        list(PREPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR}/conan)
    endif()

    # Start finding the dependencies
    if(WL_ENABLE_EIGEN3 AND NOT Eigen3_FOUND )
        find_package(Eigen3 3.3 ${REQUIRED})
        if(Eigen3_FOUND)
            target_link_libraries(deps INTERFACE Eigen3::Eigen)
        endif()
    endif()

    if(WL_ENABLE_FMT AND NOT fmt_FOUND)
        find_package(fmt 6.1.2 ${REQUIRED})
        if(fmt_FOUND)
            target_link_libraries(deps INTERFACE fmt::fmt)
        endif()
    endif()
    if(WL_ENABLE_SPDLOG AND NOT spdlog_FOUND)
        find_package(spdlog 1.5.0 ${REQUIRED})
        if(spdlog_FOUND AND TARGET spdlog AND NOT TARGET spdlog::spdlog)
            add_library(spdlog::spdlog ALIAS spdlog)
        endif()
        if(spdlog_FOUND)
            target_link_libraries(deps INTERFACE spdlog::spdlog)
        endif()
    endif()
endif()