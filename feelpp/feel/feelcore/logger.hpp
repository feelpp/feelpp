/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2025-01-11

 Copyright (C) 2025 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 
 @file logger.hpp
 @brief Spdlog-based logging system with glog-compatible interface
 
 This provides a glog-compatible logging interface using spdlog as the backend.
 
 MPI-Aware Logging:
   Control which MPI ranks produce log output via --log.mpi option:
   - "none": No ranks log (all use null sink)
   - "master" (default): Only rank 0 logs
   - "all": All ranks log to separate per-rank log files
   
 Usage:
   LOG(INFO) << "Information message";
   LOG(WARNING) << "Warning message";
   LOG(ERROR) << "Error message";
   VLOG(1) << "Verbose message at level 1";
   VLOG(2) << "More verbose message at level 2";
   CHECK(condition) << "Assertion message";
*/
#pragma once

// Include fmt formatters FIRST to ensure all formatters are available
#include <feel/feelcore/fmt.hpp>
#include <boost/mpi.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/null_sink.h>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <type_traits>

namespace Feel
{
namespace Logger
{
    /**
     * @brief Convert glog severity level string to spdlog level
     */
    inline spdlog::level::level_enum to_level(const char* t)
    {
        if (std::strcmp(t, "INFO") == 0)    return spdlog::level::info;
        if (std::strcmp(t, "WARNING") == 0) return spdlog::level::warn;
        if (std::strcmp(t, "ERROR") == 0)   return spdlog::level::err;
        if (std::strcmp(t, "FATAL") == 0)   return spdlog::level::critical;
        if (std::strcmp(t, "DEBUG") == 0)   return spdlog::level::debug;
        return spdlog::level::info;
    }
    
    /**
     * @brief Check if logging is enabled (not set to off)
     * 
     * This can be used to avoid expensive log message construction
     * when logging is disabled for this rank.
     */
    inline bool is_enabled()
    {
        return spdlog::default_logger() && 
               spdlog::get_level() != spdlog::level::off;
    }

    /**
     * @brief Get logger for console output only
     * 
     * Creates a logger that writes to stderr with colors, respecting --log.mpi setting.
     * When MPI is initialized:
     *   - log.mpi=all: All ranks print to console
     *   - log.mpi=master (default): Only rank 0 prints to console
     *   - log.mpi=none: No ranks print to console (null sink)
     * 
     * This is useful for user-facing messages that should be visible
     * on the console based on MPI rank filtering.
     * 
     * Usage: Feel::Logger::console()->info("Message: {}", value);
     */
    inline std::shared_ptr<spdlog::logger> console()
    {
        // Need to forward-declare Environment to avoid circular dependency
        // We'll check if MPI is initialized and get rank/mode dynamically
        static std::shared_ptr<spdlog::logger> logger;
        static bool initialized = false;
        
        // Lazy initialization - check MPI state when first called
        if (!initialized)
        {
            bool should_log = true;
            
            if (boost::mpi::environment::initialized())
            {
                boost::mpi::communicator world;
                should_log = (world.rank() == 0);
            }
            
            if (should_log)
            {
                auto console_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
                logger = std::make_shared<spdlog::logger>("console", console_sink);
                logger->set_level(spdlog::level::info);
                logger->set_pattern("%v"); // Simple pattern (just the message)
            }
            else
            {
                // Use null sink for ranks that shouldn't log
                auto null_sink = std::make_shared<spdlog::sinks::null_sink_st>();
                logger = std::make_shared<spdlog::logger>("console", null_sink);
                logger->set_level(spdlog::level::off);
            }
            
            initialized = true;
        }
        
        return logger;
    }

    /**
     * @brief Get logger for file output only
     * 
     * Uses the default logger which is configured to write to files.
     * This is equivalent to using LOG() macros directly.
     * 
     * Usage: Feel::Logger::file()->info("Message: {}", value);
     */
    inline std::shared_ptr<spdlog::logger> file()
    {
        return spdlog::default_logger();
    }

    /**
     * @brief Get null logger (no output)
     * 
     * Discards all log messages. Useful when you need a logger
     * interface but want to suppress output.
     * 
     * Usage: Feel::Logger::null()->info("Message: {}", value);
     */
    inline std::shared_ptr<spdlog::logger> null()
    {
        static auto logger = []() {
            auto null_sink = std::make_shared<spdlog::sinks::null_sink_st>();
            auto l = std::make_shared<spdlog::logger>("null", null_sink);
            l->set_level(spdlog::level::off);
            return l;
        }();
        return logger;
    }

    /**
     * @brief Get the default logger
     * 
     * Returns the default logger which respects --log.mpi and --log.console settings.
     * This is what LOG() macro uses internally.
     * 
     * Usage: Feel::Logger::logger()->info("Message: {}", value);
     */
    inline std::shared_ptr<spdlog::logger> logger()
    {
        return spdlog::default_logger();
    }

    /**
     * @brief Stream-based logging wrapper for spdlog
     * 
     * This provides glog-compatible << streaming syntax that accumulates
     * messages in an ostringstream and logs them via spdlog on destruction.
     * Uses pure ostream formatting to avoid fmt template instantiation issues
     * with boost::ublas types.
     */
    struct Stream
    {
        explicit Stream(spdlog::level::level_enum lvl, bool fatal = false)
            : level(lvl), fatal(fatal)
        {
        }

        ~Stream()
        {
            spdlog::log(level, "{}", oss.str());
            if (fatal)
                std::abort();
        }

        // Unified formatting approach:
        // Use fmt::format for all types, relying on:
        // 1. Built-in fmt formatters for standard types
        // 2. Module-specific formatters (feelcore/fmt.hpp, feelalg/fmt.hpp, etc.)
        // 3. fmt::streamed as fallback for types with ostream operators
        template <typename T>
        Stream& operator<<(T const& v)
        {
            if constexpr (fmt::is_formattable<T>::value)
            {
                // Type has a fmt formatter (either built-in or user-defined)
                oss << fmt::format("{}", v);
            }
            else
            {
                // Fallback to ostream operators via fmt::streamed
                oss << fmt::format("{}", fmt::streamed(v));
            }
            return *this;
        }

        // Support for stream manipulators like std::endl
        Stream& operator<<(std::ostream& (*manip)(std::ostream&))
        {
            oss << manip;
            return *this;
        }

        spdlog::level::level_enum level;
        bool fatal;
        std::ostringstream oss;
    };

    /**
     * @brief Get/set verbosity level for VLOG filtering
     * 
     * This is mapped from Environment::logVerbosityLevel() by Environment::startLogging()
     */
    inline int& verbosity()
    {
        static int v = 0;
        return v;
    }

    /**
     * @brief Set the global log level
     * @param v Verbosity level (0=info, 1=debug, 2+=trace)
     */
    inline void setLevel(int v)
    {
        if (v >= 2) 
            spdlog::set_level(spdlog::level::trace);
        else if (v >= 1) 
            spdlog::set_level(spdlog::level::debug);
        else 
            spdlog::set_level(spdlog::level::info);
    }

    /**
     * @brief Set automatic flush on specified level
     * Ensures logs are flushed to disk at the specified level or higher
     */
    inline void flushOn(int level = 0)
    {
        spdlog::flush_on(level >= 2 ? spdlog::level::trace : 
                        level >= 1 ? spdlog::level::debug : 
                        spdlog::level::info);
    }

    /**
     * @brief Create a null sink logger (discards all output)
     * @return Shared pointer to null logger
     */
    inline std::shared_ptr<spdlog::logger> createNullLogger(const std::string& name = "null")
    {
        auto null_sink = std::make_shared<spdlog::sinks::null_sink_st>();
        auto logger = std::make_shared<spdlog::logger>(name, null_sink);
        logger->set_level(spdlog::level::off);
        return logger;
    }

    /**
     * @brief Create a file sink logger
     * @param filename Log file path
     * @param truncate Whether to truncate existing file
     * @return Shared pointer to file logger
     */
    inline std::shared_ptr<spdlog::logger> createFileLogger(const std::string& name, 
                                                            const std::string& filename, 
                                                            bool truncate = false)
    {
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename, truncate);
        return std::make_shared<spdlog::logger>(name, file_sink);
    }

    /**
     * @brief Create a multi-sink logger (file + console)
     * @param name Logger name
     * @param filename Log file path
     * @param truncate Whether to truncate existing file
     * @return Shared pointer to multi-sink logger
     */
    inline std::shared_ptr<spdlog::logger> createMultiLogger(const std::string& name,
                                                             const std::string& filename,
                                                             bool truncate = false)
    {
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename, truncate);
        auto err_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
        return std::make_shared<spdlog::logger>(name, spdlog::sinks_init_list{file_sink, err_sink});
    }

    /**
     * @brief Set the default logger for the application
     * @param logger Shared pointer to logger instance
     */
    inline void setDefaultLogger(std::shared_ptr<spdlog::logger> logger)
    {
        spdlog::set_default_logger(logger);
    }

    /**
     * @brief Set the log message pattern format
     * @param pattern Format pattern string (spdlog format)
     */
    inline void setPattern(const std::string& pattern)
    {
        spdlog::set_pattern(pattern);
    }

    /**
     * @brief Disable all logging (set level to off)
     */
    inline void disable()
    {
        spdlog::set_level(spdlog::level::off);
    }

    /**
     * @brief Flush the default logger
     * Forces all buffered log messages to be written
     */
    inline void flush()
    {
        if (spdlog::default_logger())
            spdlog::default_logger()->flush();
    }

    /**
     * @brief Shutdown the logging system
     * Flushes and releases all loggers
     */
    inline void shutdown()
    {
        spdlog::shutdown();
    }

    /**
     * @brief Log level enumeration
     * 
     * Provides abstraction for log levels without exposing spdlog
     */
    enum class Level
    {
        trace = 0,
        debug = 1,
        info = 2,
        warn = 3,
        err = 4,
        critical = 5,
        off = 6
    };

} // namespace Logger
} // namespace Feel

// Severity level constants for glog compatibility
namespace Feel {
    enum LogSeverity { INFO = 0, WARNING = 1, ERROR = 2, FATAL = 3 };
    
    // Make constants available in Feel namespace
    using ::Feel::INFO;
    using ::Feel::WARNING;
    using ::Feel::ERROR;
    using ::Feel::FATAL;
}

// glog-compatible macros using spdlog backend
#define LOG(SEV) ::Feel::Logger::Stream(::Feel::Logger::to_level(#SEV), std::strcmp(#SEV,"FATAL")==0)
#define LOG_IF(SEV, condition) if (condition) LOG(SEV)

// Helper macros for LOG_FIRST_N to properly expand __LINE__
#define FEELPP_LOG_COUNTER_CONCAT_IMPL(a, b) a##b
#define FEELPP_LOG_COUNTER_CONCAT(a, b) FEELPP_LOG_COUNTER_CONCAT_IMPL(a, b)
#define LOG_FIRST_N(SEV, N) \
    static int FEELPP_LOG_COUNTER_CONCAT(_log_counter_, __LINE__) = 0; \
    if (FEELPP_LOG_COUNTER_CONCAT(_log_counter_, __LINE__)++ < N) LOG(SEV)

#define VLOG(N)  if ((N) <= ::Feel::Logger::verbosity()) ::Feel::Logger::Stream(spdlog::level::debug)
#define DVLOG(N) VLOG(N)  // Debug VLOG - same as VLOG in this implementation
#define VLOG_IS_ON(N) ((N) <= ::Feel::Logger::verbosity())
#define DVLOG_IF(N, condition) if ((condition) && (N) <= ::Feel::Logger::verbosity()) ::Feel::Logger::Stream(spdlog::level::debug)

// Debug logging - only in debug builds
#ifndef NDEBUG
#define DLOG(SEV) LOG(SEV)
#define DLOG_IF(SEV, condition) LOG_IF(SEV, condition)
#else
#define DLOG(SEV) if (false) LOG(SEV)
#define DLOG_IF(SEV, condition) if (false) LOG(SEV)
#endif

// CHECK macros that abort on failure
#define CHECK(X)      if(!(X)) ::Feel::Logger::Stream(spdlog::level::critical, true) << "CHECK(" #X ") failed: "
#define CHECK_EQ(a,b) CHECK((a)==(b)) << (a) << " vs " << (b)
#define CHECK_NE(a,b) CHECK((a)!=(b)) << (a) << " vs " << (b)
#define CHECK_LE(a,b) CHECK((a)<=(b)) << (a) << " vs " << (b)
#define CHECK_LT(a,b) CHECK((a)< (b)) << (a) << " vs " << (b)
#define CHECK_GE(a,b) CHECK((a)>=(b)) << (a) << " vs " << (b)
#define CHECK_GT(a,b) CHECK((a)> (b)) << (a) << " vs " << (b)

// Debug checks - only enabled in debug builds
#ifndef NDEBUG
#define DCHECK(X) CHECK(X)
#define DCHECK_EQ(a,b) CHECK_EQ(a,b)
#define DCHECK_NE(a,b) CHECK_NE(a,b)
#define DCHECK_LE(a,b) CHECK_LE(a,b)
#define DCHECK_LT(a,b) CHECK_LT(a,b)
#define DCHECK_GE(a,b) CHECK_GE(a,b)
#define DCHECK_GT(a,b) CHECK_GT(a,b)
#else
#define DCHECK(X) if(false) ::Feel::Logger::Stream(spdlog::level::debug)
#define DCHECK_EQ(a,b) if(false) ::Feel::Logger::Stream(spdlog::level::debug)
#define DCHECK_NE(a,b) if(false) ::Feel::Logger::Stream(spdlog::level::debug)
#define DCHECK_LE(a,b) if(false) ::Feel::Logger::Stream(spdlog::level::debug)
#define DCHECK_LT(a,b) if(false) ::Feel::Logger::Stream(spdlog::level::debug)
#define DCHECK_GE(a,b) if(false) ::Feel::Logger::Stream(spdlog::level::debug)
#define DCHECK_GT(a,b) if(false) ::Feel::Logger::Stream(spdlog::level::debug)
#endif
