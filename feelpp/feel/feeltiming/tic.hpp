/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  SPDX-FileCopyrightText: 2012 Université Joseph Fourier (Grenoble I)
  SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
  SPDX-License-Identifier: LGPL-2.1-or-later
*/

#if !defined(FEELPP_TIMING_TIC_HPP)
#define FEELPP_TIMING_TIC_HPP 1

#include <feel/feelcore/environment.hpp>

/** 
 * @defgroup Timing
 * @ingroup Feelpp
 */
#include <feel/feeltiming/timer.hpp>
#include <feel/feeltiming/now.hpp>

namespace Feel
{

namespace details
{
struct FEELPP_NO_EXPORT SecondBasedTimer
{
    static void print( std::string const& msg, const std::pair<double,int>& val )
    {
        if ( Environment::isMasterRank() )
        {
            int cols = (val.second < 15) ? val.second : 15;
            if ( !msg.empty() )
                std::cout << std::setw(1+cols) << "[" << msg << "] Time : " << val.first << "s\n";
            else
                std::cout << std::setw(7+cols) << "Time : " << val.first << "s\n";
        }
    }
    static inline time_point  time()
    {
        return Feel::details::now();
    }
};

counter<time_point,SecondBasedTimer> const sec_timer = {};
}  // details
} // Feel

namespace Feel
{
//! display
const inline bool display = true;

//! no display
const inline bool no_display = false;

namespace time
{

//! 
/**
 * @brief Record internal time at its execution. To be used with toc.
 * 
 * \code {.cpp}
 * tic()
 * ...
 * // some code here
 * ...
 * std::cout << fmt::format("time spent in block in seconds : {}", toc() );
 * \endcode
 * 
 */
inline void tic()
{
    Feel::details::sec_timer.tic();
}



/**
 * @brief toc returns the time elapsed since the last tic
 * 
 * @param msh identifier of the timer
 * @param _display if true (default) the time is displayed
 * @param uiname user interface name
 * @return double the time elapsed since the last tic
 */
inline double  toc( std::string const& msg = "",
                    bool _display = display, 
                    std::string const& uiname = "" )
{
    auto t = Feel::details::sec_timer.toc( msg, _display );
    Environment::addTimer( msg, t, uiname );
    return t.first;
}

/**
 * @brief Gather and print the latest local toc() samples across MPI ranks.
 *
 * Call after the measured work: this function communicates, while tic() and
 * toc() remain local. Every rank in @p comm must pass identical @p labels in
 * the same order. The structured reports remain accessible through
 * Environment::timerRankReports() and Environment::timerRankReportsJson().
 *
 * @param labels Timer labels to summarize.
 * @param comm Communicator whose ranks contributed the samples.
 * @param showRankValues Print durations in rank order as well as summaries.
 */
inline void reportTocRankStatistics( std::vector<std::string> const& labels,
                                     mpi::communicator const& comm,
                                     bool showRankValues = false )
{
    Environment::gatherTimerRankStatistics( labels, comm );
    if ( comm.rank() == 0 )
        Environment::printTimerRankReports( std::cout, showRankValues );
}
} // time
} // Feel

namespace Feel
{
// Convenience namespace injection from time:: into Feel::
using time::tic;
using time::toc;
using time::reportTocRankStatistics;
}

#endif /* FEELPP_TIMING_TIC_HPP */
