/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 23 mai 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_TIMERTABLE_HPP
#define FEELPP_TIMERTABLE_HPP 1

#include <iomanip>
#include <iostream>

#include <boost/bind.hpp> 
#include <boost/ref.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/algorithm/string/trim_all.hpp>

#include <feel/feelcore/journalwatcher.hpp>

namespace Feel {

static uint16_type TIMERDATA_INSTANCE_NUMBER = 0;
static const std::string timerDefaultInstanceName() { return "timer-" + std::to_string(TIMERDATA_INSTANCE_NUMBER); }

//! TimerData is the value for the TimerTable map.
class TimerData : public std::vector<double>
{
  public:
    TimerData() = default;
    TimerData(TimerData const& ) = default;
    TimerData( std::string const& m, 
               std::string const& p = timerDefaultInstanceName() )
        : msg(m), instance_name(p)
    {
        TIMERDATA_INSTANCE_NUMBER++;
    }
    void add( std::pair<double,int> const& t ) { this->push_back( t.first ); level=t.second; }
    std::string msg;
    std::string instance_name;
    int level;
};

//! TimerTable is a map of timer.
class TimerTable : std::map<std::string, TimerData>,
    public JournalWatcher
{
    using super = std::map<std::string, TimerData>;
    using super2 = JournalWatcher;
public:

    //! Constructors
    //! @{

    //! Default constructor
    TimerTable() : super2( "TimerTable", false ) {}


    //! @}

    //! Destructor
    ~TimerTable() override { this->journalFinalize(); }

    //! Modifiers
    //! @{

    //! Add a timer for the given message
    void add( std::string const& msg,
              std::pair<double,int> const&  t,
              std::string const& uiname = "" )
        {
            if ( msg.empty() ) return;
            auto it = this->find( msg );
            if ( it != this->end() )
            {
                it->second.add( t );
            }
            else
            {
                
                TimerData T( msg );
                if( not uiname.empty() )
                    T.instance_name = uiname;
                T.add( t );
                this->insert( std::make_pair( msg, T ) );
                auto m = msg.size()+2*t.second;
                M_max_len = (m>M_max_len)?m:M_max_len;
            }
        }

    //! @}

    //! Methods
    //! @{

    //! Save the timers.
    //! If display is true, the timers are printed on the standard output.
    void save( bool display ) 
        {
            std::ostringstream os;
            using namespace boost::accumulators;
            using namespace std::string_literals;
            
            os << std::setw( M_max_len ) <<std::left << "Timer" << " " 
               << std::setw(7) << std::right << "Count" << " "
               << std::setw(9+2) << std::right << "Total(s)" << " "
               << std::setw(9+2) << std::right << "Max(s)" << " "
               << std::setw(9+2) << std::right << "Min(s)" << " "
               << std::setw(9+2) << std::right << "Mean(s)" << " "
               << std::setw(9+2) << std::right << "StdDev(s)" << "\n";
            
            std::map<double,TimerData, std::greater<double>> sortTotal;
            for( auto const& T: *this )
            {
                accumulator_set<double, stats<boost::accumulators::tag::count,
                                              boost::accumulators::tag::mean,
                                              boost::accumulators::tag::variance,
                                              boost::accumulators::tag::min,
                                              boost::accumulators::tag::max> > acc;
                for_each(T.second.begin(), T.second.end(), boost::bind<void>(boost::ref(acc), _1));
                double tot = boost::accumulators::sum(acc);
                sortTotal[tot] = T.second;
                os << std::setw( 2*T.second.level ) << " "
                   << std::setw( M_max_len -  2*T.second.level ) <<std::left << T.first << " " 
                   << std::setw(7) << std::right << boost::accumulators::count(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::sum(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::max(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::min(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::mean(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << std::sqrt(boost::accumulators::variance(acc)) << "\n";
            }
            os << "--------------------------------------------------------------------------------\n";
            os << std::setw( M_max_len ) <<std::left << "Timer" << " " 
               << std::setw(7) << std::right << "Count" << " "
               << std::setw(9+2) << std::right << "Total(s)" << " "
               << std::setw(9+2) << std::right << "Max(s)" << " "
               << std::setw(9+2) << std::right << "Min(s)" << " "
               << std::setw(9+2) << std::right << "Mean(s)" << " "
               << std::setw(9+2) << std::right << "StdDev(s)" << "\n";
            for( auto const& T: sortTotal )
            {
                accumulator_set<double, stats<boost::accumulators::tag::count,
                                              boost::accumulators::tag::mean,
                                              boost::accumulators::tag::variance,
                                              boost::accumulators::tag::min,
                                              boost::accumulators::tag::max> > acc;
                for_each(T.second.begin(), T.second.end(), boost::bind<void>(boost::ref(acc), _1));

                os << std::setw( 2*T.second.level ) << " "
                   << std::setw( M_max_len-2*T.second.level ) <<std::left << T.second.msg << " " 
                   << std::setw(7) << std::right << boost::accumulators::count(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::sum(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::max(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::min(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::mean(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << std::sqrt(boost::accumulators::variance(acc)) << "\n";

                std::string n = boost::trim_fill_copy( boost::trim_fill_copy_if( T.second.msg, "_", boost::is_any_of(":")), "_" );                
            }
            os << std::resetiosflags(os.flags());
            if ( display )
                if ( Environment::isMasterRank() )
                    std::cout << os.str() << std::endl;
        }

    //! Save timers in markdown format.
    //! \param os output stream
    //! \see save
    void saveMD( std::ostream & os ) 
        {
            using namespace boost::accumulators;
           
            os << std::setw( M_max_len ) << std::left  << "|Timer" << " " 
               << std::setw(7)           << std::right << "|Count" << " "
               << std::setw(9+2)         << std::right << "|Total(s)" << " "
               << std::setw(9+2)         << std::right << "|Max(s)" << " "
               << std::setw(9+2)         << std::right << "|Min(s)" << " "
               << std::setw(9+2)         << std::right << "|Mean(s)" << " "
               << std::setw(9+2)         << std::right << "|StdDev(s)|" << "\n";
           
            os << "|---|---|---|---|---|---|---|" << std::endl; 
            std::map<double,TimerData, std::greater<double>> sortTotal;
            for( auto const& T: *this )
            {
                accumulator_set<double, stats<boost::accumulators::tag::count,
                                              boost::accumulators::tag::mean,
                                              boost::accumulators::tag::variance,
                                              boost::accumulators::tag::min,
                                              boost::accumulators::tag::max> > acc;
                for_each(T.second.begin(), T.second.end(), boost::bind<void>(boost::ref(acc), _1));
                double tot = boost::accumulators::sum(acc);
                sortTotal[tot] = T.second;
                os << std::setw( M_max_len ) << std::left       << "|" << T.first << " | " 
                   << std::setw(7)           << std::right      << boost::accumulators::count(acc) << " | "
                   << std::setw(11)          << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::sum(acc) << " | "
                   << std::setw(11)          << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::max(acc) << " | "
                   << std::setw(11)          << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::min(acc) << " | "
                   << std::setw(11)          << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::mean(acc) << " | "
                   << std::setw(11)          << std::scientific << std::setprecision( 2 ) << std::right << std::sqrt(boost::accumulators::variance(acc)) << "|\n";
            }
            os << "\n";
            os << std::setw( M_max_len ) << std::left  << "|Timer" << " " 
               << std::setw(7)           << std::right << "|Count" << " "
               << std::setw(9+2)         << std::right << "|Total(s)" << " "
               << std::setw(9+2)         << std::right << "|Max(s)" << " "
               << std::setw(9+2)         << std::right << "|Min(s)" << " "
               << std::setw(9+2)         << std::right << "|Mean(s)" << " "
               << std::setw(9+2)         << std::right << "|StdDev(s)|" << "\n";
            os << "|---|---|---|---|---|---|---|" << std::endl; 
            for( auto const& T: sortTotal )
            {
                accumulator_set<double, stats<boost::accumulators::tag::count,
                                              boost::accumulators::tag::mean,
                                              boost::accumulators::tag::variance,
                                              boost::accumulators::tag::min,
                                              boost::accumulators::tag::max> > acc;
                for_each(T.second.begin(), T.second.end(), boost::bind<void>(boost::ref(acc), _1));

                os << std::setw( M_max_len ) <<std::left << "|" << T.second.msg << " | " 
                   << std::setw(7) << std::right << boost::accumulators::count(acc) << " | "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::sum(acc) << " | "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::max(acc) << " | "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::min(acc) << " | "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << boost::accumulators::mean(acc) << " | "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << std::sqrt(boost::accumulators::variance(acc)) << "|\n";
            }
            os << "\n";
            os << std::resetiosflags(os.flags());
        }

    void test()
    {}

private:

    //! Private Methods
    //! @{

    //! update information
    void updateInformationObject( nl::json & p ) const override
    {
        using namespace boost::accumulators;
        using namespace std::string_literals;

        p.clear();
        std::map<double,TimerData, std::greater<double>> sortTotal;

        // Get total time accumulating all timers.
        for( auto const& T: *this )
        {
            accumulator_set<double, stats<boost::accumulators::tag::count,
                                          boost::accumulators::tag::mean,
                                          boost::accumulators::tag::variance,
                                          boost::accumulators::tag::min,
                                          boost::accumulators::tag::max> > acc;
            for_each(T.second.begin(), T.second.end(), boost::bind<void>(boost::ref(acc), _1));
            double tot = boost::accumulators::sum(acc);
            sortTotal[tot] = T.second;

            const std::string& prefix = T.second.instance_name;
            auto & jData = p[prefix];
            jData[ "message" ] = T.second.msg;
            jData[ "count" ] = boost::accumulators::count(acc);
            jData[ "total" ] = boost::accumulators::sum(acc);
            jData[ "max" ] = boost::accumulators::max(acc);
            jData[ "min"] = boost::accumulators::min(acc);
            jData[ "mean" ] = boost::accumulators::mean(acc);
            jData[ "std_dev" ] = std::sqrt( boost::accumulators::variance(acc) );
        }

        for( auto const& T: sortTotal )
        {
            accumulator_set<double, stats<boost::accumulators::tag::count,
                                          boost::accumulators::tag::mean,
                                          boost::accumulators::tag::variance,
                                          boost::accumulators::tag::min,
                                          boost::accumulators::tag::max> > acc;
            for_each(T.second.begin(), T.second.end(), boost::bind<void>(boost::ref(acc), _1));

            const std::string& prefix = T.second.instance_name;
            auto & jData = p[prefix];
            jData[ "message" ] = T.second.msg;
            jData[ "count" ] = boost::accumulators::count(acc);
            jData[ "total" ] = boost::accumulators::sum(acc);
            jData[ "max" ] = boost::accumulators::max(acc);
            jData[ "min" ] = boost::accumulators::min(acc);
            jData[ "mean" ] = boost::accumulators::mean(acc);
            jData[ "std_dev" ] = std::sqrt( boost::accumulators::variance(acc) );
        }
    }

    //! @}

private:
    size_type M_max_len = 0;
};

} // Feel namespace.
#endif
