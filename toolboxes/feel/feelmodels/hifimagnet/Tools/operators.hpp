/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:f
enc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Author(s): Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   Date: 2011-16-12

   Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)
   Copyright (C) CNRS

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
   \file operators.h
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2013-06-03
*/

#ifndef __HIFIMAGNET_OPERATORS_HPP
#define __HIFIMAGNET_OPERATORS_HPP 1

#include <feel/feelalg/vectorublas.hpp>

////////// unary element_wise ///////////
namespace Feel{

template <typename T, typename R, typename I>
VectorUblas<T>
unary_elementwise( VectorUblas<T> const& v1, R (*unary_func)(I) )
{
    typedef typename type_traits<T>::real_type real_type;

    if( !boost::is_convertible<R,real_type>::value || !boost::is_convertible<I,real_type>::value )
        throw std::logic_error( "Conflicting types : unary function can't be applied to vector elements" );

    VectorUblas<real_type> _t( v1.mapPtr() );
    size_type s = v1.localSize();
    size_type start = v1.firstLocalIndex();

    for ( size_type i = 0; i < s; ++i )
        _t.operator()( start+i ) = (*unary_func)(v1.operator()( start + i ));

    return _t;
}

template <typename T, typename R, typename I>
VectorUblas<T>
unary_elementwise( boost::shared_ptr<VectorUblas<T> > const& v1, R (*unary_func)(I) )
{
    return unary_elementwise( *v1, &unary_func );
}

/////////////////////////////
template <typename T, typename R, typename I1, typename I2>
VectorUblas<T>
binary_elementwise( VectorUblas<T> const& v1, VectorUblas<T> const& v2,
                    R (*binary_func)(I1, I2) )
{
    FEELPP_ASSERT( v1.localSize() == v2.localSize() &&
                   v1.size() == v2.size() )
        ( v1.localSize() )( v2.localSize() )
        ( v1.size() )( v2.size() ).error( "incompatible vector sizes" );

    typedef typename type_traits<T>::real_type real_type;

    if( !boost::is_convertible<R,real_type>::value || !boost::is_convertible<I1,real_type>::value || !boost::is_convertible<I2,real_type>::value)
        throw std::logic_error( "Conflicting types : binary function can't be applied to vector elements" );

    VectorUblas<real_type> _t( v1.mapPtr() );
    size_type s = v1.localSize();
    size_type start = v1.firstLocalIndex();


    for ( size_type i = 0; i < s; ++i )
        _t.operator()( start+i ) = (*binary_func)( v1.operator()( start + i ), v2.operator()( start + i ) );

    return _t;
}

template <typename T, typename R, typename I1, typename I2>
VectorUblas<T>
binary_elementwise( boost::shared_ptr<VectorUblas<T> > const& v1, boost::shared_ptr<VectorUblas<T> > const& v2,
                    R (*binary_func)(I1, I2) )
{
    return binary_elementwise( *v1, *v2, &binary_func );
}


////////// binary element_wise ///////////
}
#endif /* __HIFIMAGNET_OPERATORS_HPP 1 */
