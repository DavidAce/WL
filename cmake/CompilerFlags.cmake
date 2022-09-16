cmake_minimum_required(VERSION 3.20)
set(PROJECT_UNAME WL)

if (NOT TARGET wl-flags)
    add_library(wl-flags INTERFACE)
endif ()


###  Add optional RELEASE/DEBUG compile to flags
target_compile_options(wl-flags INTERFACE $<$<AND:$<CONFIG:DEBUG>,$<CXX_COMPILER_ID:Clang>>: -fstandalone-debug>)
target_compile_options(wl-flags INTERFACE $<$<AND:$<CONFIG:RELWITHDEBINFO>,$<CXX_COMPILER_ID:Clang>>: -fstandalone-debug>)
target_compile_options(wl-flags INTERFACE
                       $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>:/W4>
                       $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:-Wall -Wextra -Wpedantic -Wconversion -Wunused>
                       )


###  Enable c++17 support
target_compile_features(wl-flags INTERFACE cxx_std_17)

# Settings for sanitizers
if(COMPILER_ENABLE_ASAN)
    target_compile_options(wl-flags INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-fsanitize=address -fno-omit-frame-pointer>)
    target_link_libraries(wl-flags INTERFACE -fsanitize=address)
endif()
if(COMPILER_ENABLE_USAN)
    target_compile_options(wl-flags INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-fsanitize=undefined,leak,pointer-compare,pointer-subtract,alignment,bounds -fno-omit-frame-pointer>)
    target_link_libraries(wl-flags INTERFACE -fsanitize=undefined,leak,pointer-compare,pointer-subtract,alignment,bounds)
endif()

### Enable link time optimization
if(CMAKE_INTERPROCEDURAL_OPTIMIZATION)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT lto_supported OUTPUT lto_error)
    if(lto_supported)
        message(STATUS "LTO enabled")
    else()
        message(FATAL_ERROR "LTO is not supported: ${lto_error}")
    endif()
endif()
