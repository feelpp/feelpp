
// Copyright (c) 2005 Christopher Diggins
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// Disclaimer: Not a Boost library

#ifndef BOOST_PROFILER_MANAGER_HPP_INCLUDED
#define BOOST_PROFILER_MANAGER_HPP_INCLUDED

#include <iostream>
#include <map>
#include <string>
#include <utility>

#include "high_resolution_timer.hpp"

namespace boost {
namespace prof {

  using namespace std;

  // this type is used to keep running totals for computing averages
  template<typename sum_t>
  struct counted_sum : pair<int, sum_t>
  {
    typedef pair<int, sum_t> inherited_t;
    counted_sum() : inherited_t(0, 0) { }
    counted_sum(int x, sum_t y) : inherited_t(x, y) { }
    void operator+=(const sum_t& x) {
      this->first++;
      this->second += x;
    }
  };

  // profile statistics are stored in this type
  template<typename id_t, typename duration_t>
  struct default_stats_policy
  {
    typedef std::map<id_t, counted_sum<duration_t> > stats_map;
      typedef typename stats_map::const_iterator const_iterator;
            typedef typename stats_map::iterator iterator;

    static stats_map stats;

    static void on_profiled(const id_t& id, const duration_t& dur, bool underflow, bool overflow) {
      if (underflow) {
        stats[id] = counted_sum<duration_t>(-1, 0);
      }
      else
      if (overflow) {
        stats[id] = counted_sum<duration_t>(-2, 0);
      }
      else {
        stats[id] += dur;
      }
    }

    // this is a default function for generating reports
    static void generate_report(ostream& o = cerr) {
      o << "profile id" << '\t'
        << "total elapsed" << '\t'
        << "entry count" << '\t'
        << "average" << '\t'
        << "percentage" << endl;

      duration_t dTotalTotal = 0.0;
      int nTotalCount = 0;

      for (const_iterator i=stats.begin(); i != stats.end(); i++)
      {
        int nCount = i->second.first;
        nTotalCount += nCount;
        duration_t dTotal = i->second.second;
        dTotalTotal += dTotal;
      }
      for (const_iterator i=stats.begin(); i != stats.end(); i++)
      {
        int nCount = i->second.first;
        duration_t dTotal = i->second.second;
        double dAvg = dTotal / nCount;
        double dPercent = (dTotalTotal / dTotal) * 100;
        o << fixed
          << i->first << '\t'
          << dTotal << '\t'
          << nCount << '\t'
          << dAvg << '\t'
          << dPercent << endl;
      }

      o << "total entries" << nTotalCount << endl;
      o << "total time elapsed" << dTotalTotal << endl;
      o << "overall average" << dTotalTotal / nTotalCount << endl;
    }
  };

  // static initialization for stats data
  template<typename id_t, typename duration_t>
  std::map<id_t, counted_sum<duration_t> >
  default_stats_policy<id_t, duration_t>::stats;


// this logging policy causes no data to be logged
  struct empty_logging_policy
  {
    template<typename id_t>
    static void on_start(const id_t& /*id*/) { };
    template<typename id_t>
    static void on_restart(const id_t& /*id*/) { };
    template<typename id_t>
    static void on_resume(const id_t& /*id*/) { };
    template<typename id_t, typename duration_t>
    static void on_pause(const id_t& /*id*/, const duration_t& /*dur*/, bool /*underflow*/, bool /*overflow*/) { };
    template<typename id_t, typename duration_t>
    static void on_stop(const id_t& /*id*/, const duration_t& /*dur*/, bool /*underflow*/, bool /*overflow*/) { };
  };

  // this is the default logging policy
  struct default_logging_policy
  {
    template<typename id_t>
    static void on_start(const id_t& id) {
      cerr << "starting profile " << id << endl;
    }
      template<typename id_t>
    static void on_restart(const id_t& id) {
    }
    template<typename id_t>
    static void on_resume(const id_t& id) {
      cerr << "resuming profile " << id << endl;
    }
    template<typename id_t, typename duration_t>
    static void on_pause(const id_t& id, const duration_t& d, bool underflow, bool overflow) {
      cerr << "pausing profile " << id;
      cerr << " time elapsed = " << d;
      if (underflow) cerr << " underflow occurred";
      if (overflow) cerr << " overflow occurred";
      cerr << endl;
    }
    template<typename id_t, typename duration_t>
    static void on_stop(const id_t& id, const duration_t& d, bool underflow, bool overflow) {
      cerr << "stopping profile " << id;
      cerr << " time elapsed = " << d;
      if (underflow) cerr << " underflow occurred";
      if (overflow) cerr << " overflow occurred";
      cerr << endl;
    }
  };

  // this is used for managing logging responsibilities and managing stats
  // for all of the profilers
  template
  <
    typename id_t = string,
    typename duration_t = double,
    typename timer_t = high_resolution_timer,
    typename logging_policy_t = default_logging_policy,
    typename stats_policy_t = default_stats_policy<id_t, duration_t>
  >
  struct basic_profile_manager
  {
    typedef id_t id_type;
    typedef duration_t duration_type;
    typedef timer_t timer_type;

    typedef logging_policy_t logging_policy;
    typedef stats_policy_t stats_policy;

    static void on_start(const id_t& id) {
      logging_policy::on_start(id);
    }
    static void on_restart(const id_t& id) {
      logging_policy::on_restart(id);
    }
    static void on_resume(const id_t& id) {
      logging_policy::on_resume(id);
    }
    static void on_pause(const id_t& id, const duration_t& dur, bool underflow, bool overflow) {
      logging_policy::on_pause(id, dur, underflow, overflow);
      stats_policy::stats[id] += dur;
    }
    static void on_stop(const id_t& id, const duration_t& dur, bool underflow, bool overflow) {
      logging_policy::on_stop(id, dur, underflow, overflow);
      stats_policy::on_profiled(id, dur, underflow, overflow);
    }
    static void generate_report() {
      stats_policy::generate_report();
    }
  };


  typedef basic_profile_manager<> profiler_manager;
} }

#endif
