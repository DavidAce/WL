#pragma once
#include <memory>
#include <fmt/compile.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace wl {
    inline std::shared_ptr<spdlog::logger> log;
    inline void enableTimestamp() {
        if(not log) return;
        log->trace("Enabled timestamp");
        log->set_pattern("[%Y-%m-%d %H:%M:%S][%n]%^[%=8l]%$ %v");
    }
    inline void disableTimestamp() {
        if(not log) return;
        log->trace("Disabled timestamp");
        log->set_pattern("[%n]%^[%=8l]%$ %v");
    }

    inline void setLogger(const std::string &name, spdlog::level::level_enum level = spdlog::level::info, bool timestamp = false) {
        if(spdlog::get(name) == nullptr)
            log = spdlog::stdout_color_mt(name, spdlog::color_mode::automatic);
        else
            log = spdlog::get(name);
        log->set_pattern("[%n]%^[%=8l]%$ %v"); // Disabled timestamp is the default
        log->set_level(level);
        if(timestamp) enableTimestamp();
    }

}