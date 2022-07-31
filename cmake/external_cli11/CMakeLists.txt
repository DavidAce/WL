cmake_minimum_required(VERSION 3.15)
project(external-cli11)
if(INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif()
include(ExternalProject)
ExternalProject_Add(external_cli11
        URL  https://github.com/CLIUtils/CLI11/archive/refs/tags/v2.2.0.tar.gz
        URL_MD5 f0b70187eadf05a23257ac1d17c9a0e3
        PREFIX      ${CMAKE_BINARY_DIR}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}

        #Give Ninja access to the terminal.
        USES_TERMINAL_DOWNLOAD TRUE
        USES_TERMINAL_UPDATE TRUE
        USES_TERMINAL_CONFIGURE TRUE
        USES_TERMINAL_BUILD TRUE
        USES_TERMINAL_INSTALL TRUE
        USES_TERMINAL_TEST TRUE

        BUILD_ALWAYS TRUE
        CMAKE_ARGS
        ${INIT_CACHE_ARG}
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_INSTALL_MESSAGE=NEVER
        # cli11 flags
        -DCLI11_BUILD_EXAMPLES:BOOL=OFF
        -DCLI11_BUILD_TESTS:BOOL=OFF
        -DCLI11_BUILD_DOC:BOOL=OFF
        )