#include "algorithm/bootstrap.h"
#include "algorithm/simulation.h"
#include <algorithm/class_WL_worker.h>
#include <iostream>
#include <mpi.h>
#include <params/nmspc_WL_constants.h>
#include <general/nmspc_logger.h>
#include <CLI/CLI.hpp>

template<typename T>
constexpr auto sv2enum(std::string_view item) {
    if constexpr (std::is_same_v<T, spdlog::level::level_enum>){
        if(item == "trace") return spdlog::level::level_enum::trace;
        if(item == "debug") return spdlog::level::level_enum::debug;
        if(item == "info") return spdlog::level::level_enum::info;
        if(item == "warn") return spdlog::level::level_enum::warn;
        if(item == "error") return spdlog::level::level_enum::err;
        if(item == "critical") return spdlog::level::level_enum::critical;
        if(item == "off") return spdlog::level::level_enum::off;    }
    throw std::runtime_error("sv2enum given invalid string item: " + std::string(item));
}


template<typename T, auto num>
using enumarray_t = std::array<std::pair<std::string, T>, num>;

template<typename T, typename... Args>
constexpr auto mapStr2Enum(Args... names) {
    constexpr auto num     = sizeof...(names);
    auto           pairgen = [](const std::string &name) -> std::pair<std::string, T> {
        return {name, sv2enum<T>(name)};
    };
    return enumarray_t<T, num>{pairgen(names)...};
}

// MWE: https://godbolt.org/z/jddxod53d
int parse(int argc, char **argv) {
    using namespace spdlog;
    auto s2e_log     = mapStr2Enum<spdlog::level::level_enum>("trace", "debug", "info", "warn", "error", "critical","off");

    CLI::App app;
    app.description("WL: Parallel Wang-Landau 1/t algorithm");
    app.get_formatter()->column_width(90);
    app.option_defaults()->always_capture_default();
    app.allow_extras(false);
    /* clang-format off */
    app.add_option("-L, --lattice-size"                , constants::L                   , "Linear size of the lattice");
    app.add_option("--minimum-lnf"                     , constants::minimum_lnf         , "Convergence threshold for the modification factor lnf");
    app.add_option("--num-teams"                       , constants::num_teams           , "Number of energy sub-windows (if num_teams < MPI world size, take the latter)");
    app.add_option("-v,--log,--verbosity,--loglevel"   , constants::loglevel            , "Log level of WL")->transform(CLI::CheckedTransformer(s2e_log, CLI::ignore_case))->type_name("ENUM");
    app.add_flag("--timestamp,--logstamp"              , constants::logstamp            , "Prepend timestamp to log messages");
    /* clang-format on */

    app.parse(argc, argv);
    return 0;
}


int main(int argc, char *argv[]) {
    parse(argc, argv);


    // Initialize the MPI environment
    MPI_Init(nullptr, nullptr);
    int world_ID, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_ID);   // Establish thread number of this worker
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Get total number of threads
    constants::num_teams = std::min(world_size, constants::num_teams);
    class_worker worker(world_ID, world_size);

    wl::setLogger(fmt::format("WL-{:2}", world_ID), constants::loglevel, constants::logstamp);

    do_simulations(worker);
    do_bootstrap(worker);
    do_sampling(worker);
    if(world_ID == 0) { std::cout << "Finished successfully" << std::endl; }

    MPI_Finalize();
    return 0;
}