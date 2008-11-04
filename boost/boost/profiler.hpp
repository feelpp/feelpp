
// Copyright (c) 2005 Christopher Diggins
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// Disclaimer: Not a Boost library

#ifndef BOOST_PROFILER_HPP_INCLUDED
#define BOOST_PROFILER_HPP_INCLUDED

#include "profile_manager.hpp"

namespace boost {
namespace prof {

  using namespace std;

  #ifndef PROFILING_OFF
    template<typename policy_t>
    class basic_profiler
    {
      typedef typename policy_t::timer_type     timer_type;
      typedef typename policy_t::duration_type  duration_type;
      typedef typename policy_t::id_type        id_type;

      public:
        basic_profiler()
          : underflow(false), overflow(false), timing(true), elapsed(0.0)
        {
          policy_t::on_start(id);
          t.restart();
        }
        basic_profiler(const id_type& x)
          : underflow(false), overflow(false), timing(true), id(x), elapsed(0.0)
        {
          policy_t::on_start(id);
          t.restart();
        }
        ~basic_profiler() {
          if (timing) {
            stop();
          }
          policy_t::generate_report();
        }
        void stop() {
          duration_type tmp = t.elapsed();
          if (tmp <= t.elapsed_min()) {
            underflow = true;
          }
          if (tmp >= t.elapsed_max()) {
            overflow = true;
          }
          tmp += elapsed;
          timing = false;
          policy_t::on_stop(id, tmp, underflow, overflow);
        }
        void restart() {
          timing = true;
          elapsed = 0.0;
          policy_t::on_restart(id);
          t.restart();
        }
        void resume() {
          timing = true;
          policy_t::on_resume(id);
          t.restart();
        }
        void pause() {
          duration_type tmp = t.elapsed();
          if (tmp <= t.elapsed_min()) {
            underflow = true;
          }
          if (tmp >= t.elapsed_max()) {
            overflow = true;
          }
          elapsed += tmp;
          timing = false;
          policy_t::on_pause(id, tmp, underflow, overflow);
        }
      public:
        // extra functionality, may or may not stay in
        duration_type& recorded_elapsed() {
          return policy_t::stats[id].first;
        }
        duration_type& elapsed_so_far() {
          return elapsed;
        }
        int& recorded_count() {
          return policy_t::stats[id].second;
        }
      private:
        bool underflow;
        bool overflow;
        bool timing;
        id_type id;
        duration_type elapsed;
        timer_type t;
    };
  #else
    template<typeid logging_policy, typeid stats_policy, typeid timer_t>
    class basic_profiler {
      public:
        basic_profiler() { }
        void stop() { }
        void restart() { }
        void resume() { }
        void pause() { }
    };
  #endif

  typedef basic_profiler<profiler_manager> profiler;

} // namespace prof
} // namespace boost

#endif // #ifndef BOOST_PROFILER_HPP_INCLUDED
