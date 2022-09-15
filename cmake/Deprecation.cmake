
if(WL_DOWNLOAD_METHOD)
    message(FATAL_ERROR "The variable [WL_DOWNLOAD_METHOD] has been deprecated. Replace by:\n"
            "WL_PACKAGE_MANAGER:STRING=[find|cmake|find-or-cmake|conan]")
endif()

if(WL_DEPS_IN_SUBDIR)
    message(FATAL_ERROR "The option [WL_DEPS_IN_SUBDIR] has been deprecated. Replace by:\n"
            "WL_PREFIX_ADD_PKGNAME:BOOL=[TRUE|FALSE]")
endif()

if(WL_PRINT_INFO)
    message(FATAL_ERROR "The option [WL_PRINT_INFO] has been deprecated. Replace by:\n"
            "CMake CLI option --loglevel=[TRACE|DEBUG|VERBOSE|STATUS...] or\n"
            "set CMAKE_MESSAGE_LOG_LEVEL=[TRACE|DEBUG|VERBOSE|STATUS...]"
            )
endif()

if(WL_MICROARCH)
    message(FATAL_ERROR "The option [WL_MICROARCH] has been deprecated.\n")
endif()