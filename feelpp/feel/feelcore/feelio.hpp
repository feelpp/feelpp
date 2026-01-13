/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 19 Feb 2016

 Copyright (C) 2016 Feel++ Consortium

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
 */
#ifndef FEELPP_FEELIO_HPP
#define FEELPP_FEELIO_HPP 1

#include <iostream>
#include <fstream>
#include <sstream>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/logger.hpp>
#include <feel/feelcore/fmt.hpp>


namespace Feel {

/**
 * Output Stream that outputs only on master rank of a worldcomm
 *
 * The first application is to use the instantiation of MasterStream cout, cerr
 * or clog to default output on master rank process of a Feel++ application
 * @code
 * Environment env(...); // initialize Feel++ environment
 * // cout only on master rank process
 * cout << "Hello World from process " << Environment::rank()  << std::endl;
 * @endcode
 *
 * Behavior based on MPI state:
 * - Before MPI initialization: outputs on all processes (typically just one)
 * - After MPI initialization: outputs only on master rank
 * - After MPI finalization: outputs on all processes (to avoid MPI crashes)
 * - No MPI available: always outputs
 *
 * Logging integration:
 * - When logging is enabled (via spdlog or glog), output is also sent to LOG(INFO/WARNING/ERROR)
 * - This provides automatic logging of all Feel::cout/cerr output
 * - Logging happens on flush (std::endl, std::flush, or explicit flush())
 *
 * This ensures MasterStream is safe to use throughout the entire program lifecycle,
 * including during global object construction/destruction.
 */
class FEELPP_EXPORT MasterStream
{
public:
    using LogLevel = Logger::Level;

    /**
     * Construct a MasterStream from a std::ostream and optionally a WorldComm.
     * If no WorldComm is provided, it will be retrieved from Environment when needed.
     * 
     * @param _out the underlying output stream (std::cout, std::cerr, etc.)
     * @param _wc optional WorldComm pointer (defaults to nullptr, retrieved lazily)
     * @param _log_level the log level for logging backend (info, warn, err)
     */
    MasterStream(std::ostream& _out, std::shared_ptr<WorldComm> _wc = nullptr, LogLevel _log_level = Logger::Level::info )
        :
        out(_out),
        wc(_wc),
        log_level(_log_level),
        buffer()
        {}

    void attachWorldComm( std::shared_ptr<WorldComm>& _wc )
        {
            wc = _wc;
        }

private:
    /**
     * Check if we should output based on current MPI state
     * 
     * @return true if output should proceed, false otherwise
     * 
     * Logic:
     * - If MPI not initialized yet: output (return true)
     * - If MPI finalized: output to avoid crash (return true) 
     * - If no WorldComm available: output (return true)
     * - If WorldComm available and MPI active: check isMasterRank()
     */
    inline bool shouldOutput() const noexcept
        {
            // MPI not initialized or already finalized: always output
            if ( !Environment::initialized() || Environment::finalized() )
                return true;
            
            // Get WorldComm - if not set, try to get from Environment
            auto wc_ptr = wc;
            if ( !wc_ptr )
            {
                try 
                {
                    wc_ptr = Environment::worldCommPtr();
                }
                catch (const std::exception&)
                {
                    // If we can't get WorldComm, output anyway
                    return true;
                }
            }
            
            // WorldComm available and MPI is active: check master rank
            if ( wc_ptr )
            {
                try
                {
                    return wc_ptr->isMasterRank();
                }
                catch (const std::exception&)
                {
                    // If isMasterRank throws, output anyway to be safe
                    return true;
                }
            }
            
            // Fallback: output
            return true;
        }

    /**
     * Flush the internal buffer to both output stream and logging backend
     */
    void flushBuffer() const
        {
            if ( !shouldOutput() )
            {
                buffer.str("");  // Clear buffer but don't output
                buffer.clear();
                return;
            }

            std::string content = buffer.str();
            if ( content.empty() )
                return;

            // Output to stream
            out << content;
            out.flush();

            // Also log if logging is active
            try
            {
                // Remove trailing newline for LOG macro (it adds its own)
                if ( !content.empty() && content.back() == '\n' )
                    content.pop_back();
                
                // Log based on severity level
                if ( log_level == Logger::Level::info )
                    LOG(INFO) << content;
                else if ( log_level == Logger::Level::warn )
                    LOG(WARNING) << content;
                else if ( log_level == Logger::Level::err )
                    LOG(ERROR) << content;
            }
            catch (const std::exception&)
            {
                // Silently ignore logging errors to avoid cascading failures
            }

            // Clear buffer
            buffer.str("");
            buffer.clear();
        }

public:
    /**
     * Output operator for generic types
     * Buffers output and sends to both stream and logging backend on flush
     */
    template<typename T>
    const MasterStream& operator<<(const T& v) const
    {
        // Use fmt for better formatting of containers
        if constexpr (
            std::is_same_v<T, std::vector<double>> ||
            std::is_same_v<T, std::vector<float>> ||
            std::is_same_v<T, std::vector<int>> ||
            std::is_same_v<T, std::vector<std::string>>
        ) {
            buffer << fmt::format("[{}]", fmt::join(v, ", "));
        } else {
            buffer << v;
        }
        return *this;
    }
    
    /**
     * Output operator for stream manipulators (std::endl, std::flush, etc.)
     * Flushes buffer to both output stream and logging backend
     */
    MasterStream const& operator<<(std::ostream& (*F)(std::ostream&)) const
        {
            // Check if this is a flush/endl manipulator
            if ( F == static_cast<std::ostream& (*)(std::ostream&)>(std::endl) ||
                 F == static_cast<std::ostream& (*)(std::ostream&)>(std::flush) )
            {
                // Apply manipulator to buffer first
                F(buffer);
                // Then flush everything
                flushBuffer();
            }
            else
            {
                // Other manipulators just applied to buffer
                F(buffer);
            }
            return *this;
        }

    /**
     * Explicit flush method
     */
    void flush() const
        {
            flushBuffer();
        }

    /**
     * provide an interface to \c str() for ostringstream
     *
     * @return the current buffer content
     */
    std::string str() const
        {
            return buffer.str();
        }

    /**
     * @return the shared_ptr WorldComm
     */
    std::shared_ptr<WorldComm> worldCommPtr() const { return wc; }

    /**
     * Set the log level for this stream
     */
    void setLogLevel(LogLevel level)
        {
            log_level = level;
        }

    /**
     * variadic template function to provide open() from std::ofstream
     * Only opens on master rank when MPI is active
     */
    template<typename... Args>
    void open(Args... args)
        {
            auto* o = dynamic_cast<std::ofstream*>( &out );
            if ( o && shouldOutput() )
                o->open( args... );
        }

    /**
     * provides close() from std::ofstream
     * Only closes on master rank when MPI is active
     */
    void close()
        {
            flushBuffer();  // Flush before closing
            auto* o = dynamic_cast<std::ofstream*>( &out );
            if ( o && shouldOutput() )
                o->close();
        }

protected:
    std::ostream& out;                      ///< Underlying output stream
    std::shared_ptr<WorldComm> wc;          ///< WorldComm for MPI rank checking
    LogLevel log_level;                     ///< Log severity level
    mutable std::ostringstream buffer;      ///< Internal buffer for line-based output
};


extern MasterStream cout;
extern MasterStream cerr;
extern MasterStream clog;

}
#endif
