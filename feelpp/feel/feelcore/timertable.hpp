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
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/algorithm/string/trim_all.hpp>
namespace Feel {

class TimerData : public std::vector<double>
{
  public:
    TimerData() = default;
    TimerData(TimerData const& ) = default;
    TimerData( std::string const& n ) : name(n) {}
    void add( std::pair<double,int> const& t ) { this->push_back( t.first ); level=t.second; }
    std::string name;
    int level;
};
class TimerTable : std::map<std::string, TimerData>
{
public:
    TimerTable() = default;
    ~TimerTable() = default;
    
    void add( std::string const& msg, std::pair<double,int> const&  t )
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
                T.add( t );
                this->insert( std::make_pair( msg, T ) );
                auto m = msg.size()+2*t.second;
                M_max_len = (m>M_max_len)?m:M_max_len;
            }
        }
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
                double tot = sum(acc);
                sortTotal[tot] = T.second;
                os << std::setw( 2*T.second.level ) << " "
                   << std::setw( M_max_len -  2*T.second.level ) <<std::left << T.first << " " 
                   << std::setw(7) << std::right << boost::accumulators::count(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << sum(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << max(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << min(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << mean(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << sqrt(variance(acc)) << "\n";
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
                   << std::setw( M_max_len-2*T.second.level ) <<std::left << T.second.name << " " 
                   << std::setw(7) << std::right << boost::accumulators::count(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << sum(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << max(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << min(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << mean(acc) << " "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << sqrt(variance(acc)) << "\n";

                std::string n = boost::trim_fill_copy( boost::trim_fill_copy_if( T.second.name, "_", boost::is_any_of(":")), "_" );
                
                try {
                    Environment::summary().put( "application.timers."s+n+".name", T.second.name );
                    Environment::summary().put( "application.timers."s+n+".count", boost::accumulators::count(acc) );
                    Environment::summary().put( "application.timers."s+n+".total", sum(acc) );
                    Environment::summary().put( "application.timers."s+n+".max", max(acc) );
                    Environment::summary().put( "application.timers."s+n+".min", min(acc) );
                    Environment::summary().put( "application.timers."s+n+".mean", mean(acc) );
                    Environment::summary().put( "application.timers."s+n+".stddev", sqrt(variance(acc)) );
                }
                catch( pt::ptree_bad_data const& d )
                {
                    std::cout << "d: " << d.what() << std::endl;
                }
                
            }
            os << std::resetiosflags(os.flags());
            if ( display )
                if ( Environment::isMasterRank() )
                    std::cout << os.str() << std::endl;
        }
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
                double tot = sum(acc);
                sortTotal[tot] = T.second;
                os << std::setw( M_max_len ) << std::left       << "|" << T.first << " | " 
                   << std::setw(7)           << std::right      << boost::accumulators::count(acc) << " | "
                   << std::setw(11)          << std::scientific << std::setprecision( 2 ) << std::right << sum(acc) << " | "
                   << std::setw(11)          << std::scientific << std::setprecision( 2 ) << std::right << max(acc) << " | "
                   << std::setw(11)          << std::scientific << std::setprecision( 2 ) << std::right << min(acc) << " | "
                   << std::setw(11)          << std::scientific << std::setprecision( 2 ) << std::right << mean(acc) << " | "
                   << std::setw(11)          << std::scientific << std::setprecision( 2 ) << std::right << sqrt(variance(acc)) << "|\n";
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

                os << std::setw( M_max_len ) <<std::left << "|" << T.second.name << " | " 
                   << std::setw(7) << std::right << boost::accumulators::count(acc) << " | "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << sum(acc) << " | "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << max(acc) << " | "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << min(acc) << " | "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << mean(acc) << " | "
                   << std::setw(11) << std::scientific << std::setprecision( 2 ) << std::right << sqrt(variance(acc)) << "|\n";
            }
            os << "\n";
            os << std::resetiosflags(os.flags());
        }

private:
    size_type M_max_len = 0;
};



}
#endif
