/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-12-30

  Copyright (C) 2006 Université Joseph Fourier (Grenoble)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file bdf.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-12-30
 */
#include <stdexcept>
#include <sstream>

#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/utility.hpp>

#include <feel/feelalg/glas.hpp>
#include <feel/feeldiscr/bdf.hpp>

namespace Feel
{
Bdf::Bdf( const UInt n )
    :
    _M_order( n ),
    _M_size( 0 ),
    _M_alpha( n + 1 ),
    _M_beta( n )
{
    if ( n <= 0 || n > BDF_MAX_ORDER )
    {
        std::ostringstream __ex;
        __ex << "Error: wrong BDF order\n"
             << " you want to use BDF order " << n << "\n"
             << " we support BDF order from 1 to " << BDF_MAX_ORDER << "\n";
        throw std::invalid_argument( __ex.str() );
    }

    switch ( n )
    {
    case 1:
        _M_alpha[ 0 ] = 1.; // Backward Euler
        _M_alpha[ 1 ] = 1.;
        _M_beta[ 0 ] = 1.; // u^{n+1} \approx u^n
        break;

    case 2:
        _M_alpha[ 0 ] = 3. / 2.;
        _M_alpha[ 1 ] = 2.;
        _M_alpha[ 2 ] = -1. / 2.;
        _M_beta[ 0 ] = 2.;
        _M_beta[ 1 ] = -1.;
        break;

    case 3:
        _M_alpha[ 0 ] = 11. / 6.;
        _M_alpha[ 1 ] = 3.;
        _M_alpha[ 2 ] = -3. / 2.;
        _M_alpha[ 3 ] = 1. / 3.;
        _M_beta[ 0 ] = 3.;
        _M_beta[ 1 ] = -3.;
        _M_beta[ 2 ] = 1.;
        break;
    }

    _M_unknowns.resize( n );
}



Bdf::~Bdf()
{}


void
Bdf::initialize( Vector const& u0 )
{
    std::for_each( _M_unknowns.begin(), _M_unknowns.end(), boost::lambda::_1 = u0 );
}


void
Bdf::initialize( unknowns_type const& uv0 )
{
    // Check if uv0 has the right dimensions
    ASSERT( uv0.size() < _M_order, "Initial data set are not enough for the selected BDF" );
    ASSERT( uv0.size() > _M_order, "Initial data set is too large for the selected BDF" );

    std::copy( uv0.begin(), uv0.end(), _M_unknowns.begin() );
}
const
Bdf::unknowns_type&
Bdf::unknowns() const
{
    return _M_unknowns;
}


double
Bdf::derivateCoefficient( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i >= 0 && i < _M_order + 1,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_alpha[ i ];
}

double
Bdf::extrapolateCoefficient( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i >= 0 && i < _M_order,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_beta[ i ];
}

void
Bdf::showMe( std::ostream& __out ) const
{
    __out << "*** BDF Time discretization of order " << _M_order << " ***"
          << "  size : " << _M_unknowns[0].size() << "\n"
          << " alpha : " << _M_alpha << "\n"
          << "  beta : " << _M_beta << "\n";

}

void
Bdf::shiftRight( Vector const& __new_unk )
{
    using namespace boost::lambda;
    unknowns_type::reverse_iterator __it = boost::next( _M_unknowns.rbegin() );
    std::for_each( _M_unknowns.rbegin(), boost::prior( _M_unknowns.rend() ),
                   ( _1 = *var( __it ), ++var( __it ) ) );
    // u(t^{n}) coefficient is in _M_unknowns[0]
    _M_unknowns[0] = __new_unk;
}


Vector
Bdf::derivate( Real dt ) const
{
    return derivate()/dt;
}

Vector
Bdf::derivate() const
{
    Vector __t( _M_unknowns[0].size() );
    __t = ZeroVector( __t.size() );

    for ( UInt i = 0; i < _M_order; ++i )
        __t += _M_alpha[ i+1 ]*_M_unknowns[i];

    return __t;
}

Vector
Bdf::extrapolate() const
{
    Vector __t( _M_unknowns[0].size() );
    __t = ZeroVector( __t.size() );

    for ( UInt i = 0; i < _M_order; ++i )
        __t +=  _M_beta[ i ] * _M_unknowns[ i ];

    return __t;
}
}
