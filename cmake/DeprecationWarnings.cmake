
if(WL_DOWNLOAD_METHOD)
    message(WARNING "The variable [WL_DOWNLOAD_METHOD] has been deprecated\n"
            "Use the following variable instead:\n"
            "\t WL_PACKAGE_MANAGER:STRING=[find|cmake|find-or-cmake|conan]")
    set(WL_PACKAGE_MANAGER ${WL_DOWNLOAD_METHOD} CACHE STRING "")
endif()

if(WL_DEPS_IN_SUBDIR)
    message(WARNING "The option [WL_DEPS_IN_SUBDIR] has been deprecated\n"
            "Use the following variable instead:\n"
            "\t WL_PREFIX_ADD_PKGNAME:BOOL=[TRUE|FALSE]")
    set(WL_PREFIX_ADD_PKGNAME ${WL_DEPS_IN_SUBDIR})
endif()

if(WL_PRINT_INFO)
    message(WARNING "The option [WL_PRINT_INFO] has been deprecated\n"
            "Use the built-in CMake CLI option --loglevel=[TRACE|DEBUG|VERBOSE|STATUS...] instead")
endif()