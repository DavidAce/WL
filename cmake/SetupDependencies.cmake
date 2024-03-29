
if(NOT TARGET wl-deps)
    add_library(wl-deps INTERFACE)
endif()
include(cmake/SetupPaths.cmake) # Setup paths that find_package should search

# Link MPI
set(MPI_CXX_LINK_FLAGS "-lopen-pal")
find_package(MPI COMPONENTS C CXX REQUIRED)
target_link_libraries(wl-deps INTERFACE MPI::MPI_CXX MPI::MPI_C  ${MPI_CXX_LINK_FLAGS})

# Link std::filesystem
find_package(Filesystem COMPONENTS Final Experimental REQUIRED)
target_link_libraries(wl-deps INTERFACE std::filesystem)



# Install dependencies
include(cmake/SetupDependenciesCMake.cmake)
include(cmake/SetupDependenciesCPM.cmake)
include(cmake/SetupDependenciesConan.cmake)
include(cmake/SetupDependenciesFind.cmake)


