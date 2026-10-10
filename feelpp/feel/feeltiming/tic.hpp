/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  SPDX-FileCopyrightText: 2012 Université Joseph Fourier (Grenoble I)
  SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
  SPDX-License-Identifier: LGPL-2.1-or-later
*/

#if !defined(FEELPP_TIMING_TIC_HPP)
#define FEELPP_TIMING_TIC_HPP 1

#include <feel/feelcore/environment.hpp>
#include <concepts>
#include <type_traits>
#include <utility>

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
 * @param msg Identifier of the timer.
 * @param _display if true (default) the time is displayed
 * @param uiname User interface name stored with the timer.
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
 * @brief Conditionally record consecutive phases into a profiling struct.
 *
 * Constructing with a null profile disables timing. Each lap records one
 * duration under its label, adds it to the selected profile member, and starts
 * the next phase. finish() records the last phase and disables the timer.
 * This stopwatch has its own start time, so it does not change the tic()/toc()
 * stack used by enclosing timers.
 *
 * @code{.cpp}
 * ProfileLapTimer timer{ profile }; // profile may be null
 * // Work on the metadata phase.
 * timer.lap( "metadata", &Profile::metadataSeconds );
 * // Work on the insertion phase.
 * timer.finish( "insert", &Profile::insertSeconds );
 * @endcode
 *
 * @tparam Profile Profiling struct with floating-point duration members.
 */
template <class Profile>
    requires std::is_class_v<Profile>
class ProfileLapTimer
{
public:
    /// Type of profiling struct updated by this timer.
    using profile_type = Profile;

    /**
     * @brief Start timing when @p profile is nonnull.
     * @param profile Destination for phase totals, or null to disable timing.
     */
    explicit ProfileLapTimer( Profile* profile )
        : M_Profile( profile ),
          M_Start( profile ? Feel::details::now() : Feel::time_point{} ),
          M_Level( profile ? static_cast<int>( Feel::details::sec_timer.times().size() ) : 0 )
    {}

    /** @brief Prevent two timers from recording the same phase. */
    ProfileLapTimer( ProfileLapTimer const& ) = delete;
    /** @brief Prevent two timers from recording the same phase. */
    ProfileLapTimer& operator=( ProfileLapTimer const& ) = delete;

    /**
     * @brief Record the current phase and start timing the next one.
     * @tparam ValueType Floating-point type of the profile duration member.
     * @param label TimerTable label for this phase.
     * @param member Profile duration member to increment.
     */
    template <std::floating_point ValueType>
    void lap( std::string const& label, ValueType Profile::* member )
    {
        record( label, member );
    }

    /**
     * @brief Record the final phase and disable this timer.
     * @tparam ValueType Floating-point type of the profile duration member.
     * @param label TimerTable label for this phase.
     * @param member Profile duration member to increment.
     */
    template <std::floating_point ValueType>
    void finish( std::string const& label, ValueType Profile::* member )
    {
        if ( !M_Profile )
            return;
        record( label, member );
        M_Profile = nullptr;
    }

private:
    /**
     * @brief Store one phase duration and reset the start time.
     * @tparam ValueType Floating-point type of the profile duration member.
     * @param label TimerTable label for this phase.
     * @param member Profile duration member to increment.
     */
    template <std::floating_point ValueType>
    void record( std::string const& label, ValueType Profile::* member )
    {
        if ( !M_Profile )
            return;
        double seconds = std::chrono::duration<double>( Feel::details::now() - M_Start ).count();
        Environment::addTimer( label, { seconds, M_Level }, "" );
        M_Profile->*member += static_cast<ValueType>( seconds );
        M_Start = Feel::details::now();
    }

    /// Profile to update, or null when timing is disabled or complete.
    Profile* M_Profile = nullptr;
    /// Start of the current phase.
    Feel::time_point M_Start{};
    /// Nesting level of the enclosing tic()/toc() stack.
    int M_Level = 0;
};

/**
 * @brief Gather pooled and per-rank distributions of selected ranks' toc() calls.
 *
 * Every rank in @p comm must pass the same labels in the same order. The
 * quantiles are approximate and are computed from bounded local sketches.
 *
 * @param labels Timer labels to profile.
 * @param comm Communicator containing the sampled ranks.
 * @param showRanks Print each rank's profile after the pooled profile.
 * @param selectedRanks Rank IDs to include; defaults to {0}, while {} includes all ranks.
 */
inline void reportTocProfileStatistics( std::vector<std::string> const& labels,
                                        mpi::communicator const& comm,
                                        bool showRanks = false,
                                        std::vector<int> selectedRanks = { 0 } )
{
    Environment::gatherTimerProfileStatistics( labels, comm, std::move( selectedRanks ) );
    if ( comm.rank() == 0 )
        Environment::printTimerProfileReports( std::cout, showRanks );
}

/**
 * @brief Discover all labels on selected ranks and print their call profiles.
 *
 * Every rank in @p comm must call this with the same rank selection. The
 * default selects communicator rank 0; pass {} to include all ranks and the
 * union of their labels.
 *
 * @param comm Communicator containing the sampled ranks.
 * @param showRanks Print each selected rank's profile after its pooled profile.
 * @param selectedRanks Rank IDs to include; defaults to {0}, while {} includes all ranks.
 */
inline void reportTocProfileStatistics( mpi::communicator const& comm,
                                        bool showRanks = false,
                                        std::vector<int> selectedRanks = { 0 } )
{
    Environment::gatherTimerProfileStatistics( comm, std::move( selectedRanks ) );
    if ( comm.rank() == 0 )
        Environment::printTimerProfileReports( std::cout, showRanks );
}
} // time
} // Feel

namespace Feel
{
// Convenience namespace injection from time:: into Feel::
using time::tic;
using time::toc;
using time::ProfileLapTimer;
using time::reportTocProfileStatistics;
}

#endif /* FEELPP_TIMING_TIC_HPP */
