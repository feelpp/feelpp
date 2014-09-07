/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-23

  Copyright (C) 2005,2006 EPFL

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
   \file fmspoint.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-23
 */
#ifndef __Fms_Point_H
#define __Fms_Point_H 1

#include <boost/numeric/ublas/storage.hpp>

namespace Feel
{
namespace details
{

template<typename T, uint16_type Dim>
class FmsPoint
    : public boost::numeric::ublas::bounded_array<T, Dim>
{
public:
    typedef T value_type;
    static const uint16_type dim = Dim;
    typedef boost::numeric::ublas::bounded_array<T, Dim> super;
    FmsPoint()
        : super( dim, value_type(0.0) ) {}
    FmsPoint( FmsPoint const& init )
        : super( init )
    {}
    FmsPoint<T, Dim>& operator=( FmsPoint<T, Dim> r )
    {
        ((super*)this)->operator=( (super)r );
        return *this;
    }
    void operator+=( FmsPoint<T, Dim> const& s )
    {
        for (uint16_type i=0; i<dim; ++i)
            this->operator[]( i ) += s[i];
    }
    void operator-=( FmsPoint<T, Dim> const& s )
    {
        for (uint16_type i=0; i<dim; ++i)
            this->operator[]( i ) -= s[i];
    }
    void operator*=( value_type f )
    {
        for (uint16_type i=0; i<dim; ++i)
            this->operator[]( i ) *= f;
    }
}; // class FmsPoint

template<typename T, uint16_type Dim>
FmsPoint<T, Dim> operator-( FmsPoint<T, Dim> const& m,
                            FmsPoint<T, Dim> const& s )
{
    FmsPoint<T, Dim> d(m);
    d -= s;
    return d;
}

template<typename T, uint16_type Dim>
T dot( FmsPoint<T, Dim> const& a,
       FmsPoint<T, Dim> const& b )
{
    T retval(0.0);
    for (uint16_type i=0; i<Dim; ++i)
        retval += a[i]*b[i];
    return retval;
}

template<typename T, uint16_type Dim>
FmsPoint<T, Dim> operator*( FmsPoint<T, Dim> const& a, T f )
{
    FmsPoint<T, Dim> retval( a );
    retval *= f;
    return retval;
}

template<typename T, uint16_type Dim>
FmsPoint<T, Dim> operator*( T f, FmsPoint<T, Dim> const& a )
{
    FmsPoint<T, Dim> retval( a );
    retval *= f;
    return retval;
}

template<typename T, uint16_type Dim>
T norm( FmsPoint<T, Dim> const& a )
{
    T norm2( dot(a,a) );
    return std::sqrt( norm2 );
}

} // namespace details

} // namespace Feel


#endif /* __Fms_Point_H */

