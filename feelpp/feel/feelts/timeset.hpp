/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 04 Oct 2015

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
#ifndef FEELPP_FEELTS_TIMESET_HPP
#define FEELPP_FEELTS_TIMESET_HPP 1

#include <iostream>
#include <tuple>
#include <vector>
#include <boost/optional.hpp>

namespace Feel {
using optional_double_t = boost::optional<double>;

/**
 * class that stores time step information:
 *  - time step
 *  - time
 *  - optionally error if available
 */
struct TimeStepInfo : public std::tuple<double, double, optional_double_t>
{
    using super = std::tuple<double, double, optional_double_t>;

    TimeStepInfo( super && s ) : super( std::forward<super>(s) ) {}
    TimeStepInfo( std::tuple<double,double,double> const& s )
        : super( std::get<0>(s), std::get<1>(s), optional_double_t(std::get<2>(s)) ) {}

    double k() const { return std::get<0>( *this ); }
    double timeStep() const { return k(); }
    double& k()  { return std::get<0>( *this ); }
    double& timeStep() { return k(); }

    double t() const { return std::get<1>( *this ); }
    double time() const { return t(); }
    double& t()  { return std::get<1>( *this ); }
    double& time()  { return t(); }

    double e() const
        {
#if BOOST_VERSION > 105500
            return std::get<2>( *this ).value_or( 0. );
#else
            auto const& a = std::get<2>( *this );
            return a?*a:0;
#endif
        }
    double error() const { return e(); }
};

/**
 * @ingroup SpaceTime
 * @class TimeSet
 * \brief Handles a set of time steps
 * 
 */
class TimeSet : public std::vector<TimeStepInfo>
{
public:
    using super = std::vector<TimeStepInfo>;
    double  kprev(int n) const
    {
        CHECK(n>=0 && this->size()>=1) << "Invalid n " << n << " or size " << this->size();
        if ( index() == 0 ) return k();
        return std::prev(this->end(), n+1 )->k();
    }
    double& kprev(int n)
    {
        CHECK(n>=0 && this->size()>=1) << "Invalid n " << n << " or size " << this->size();
        if ( index() == 0 ) return k();
        return std::prev(this->end(), n+1 )->k();
    }

    double  tprev(int n) const
    {
        CHECK(n>=0 && this->size()>=1) << "Invalid n " << n << " or size " << this->size();
        if ( index() == 0 ) return t();
        return std::prev(this->end(), n+1 )->t();
    }
    double& tprev(int n)
    {
        CHECK(n>=0 && this->size()>=1) << "Invalid n " << n << " or size " << this->size();
        if ( index() == 0 ) return t();
        return std::prev(this->end(), n+1 )->t();
    }

    double  k() const { return this->back().k(); }
    double& k()       { return this->back().k(); }
    double  t() const { return this->back().t(); }
    double& t()       { return this->back().t(); }

    void push_back( double k1, double t1 )
    {
      super::push_back( std::make_tuple( k1, t1, optional_double_t() ) );
      this->print();
    }
    void push_back( double k1, double t1, double err )
    {
      super::push_back( std::make_tuple( k1, t1, err ) );
      this->print();
    }
    size_t index() const { return this->size()-1; }

    void print();

    //! save times and time steps in \p fname
    void save( std::string const& fname = "times.tsv" );

};




} // Feel

#endif
