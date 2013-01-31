/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2013 Feel++ Consortium

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
   \file adtraits.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef __ADTraits_H
#ifdef __GNUG__
#pragma interface
#endif
#define __ADTraits_H 1

#include <complex>

#if defined( FEELPP_USES_BOOST_INTERVAL )

#ifndef BOOST_UBLAS_USE_INTERVAL
# define BOOST_UBLAS_USE_INTERVAL
#endif

#include <boost/numeric/interval.hpp>
#endif

namespace Feel
{
#if defined( FEELPP_USES_BOOST_INTERVAL )
using namespace boost::numeric;
using namespace boost::numeric::interval_lib;

template<class T>
class my_checking_base
    :
public boost::numeric::interval_lib::checking_base<T>
{
public:
    static bool is_nan( const T& x )
    {
        //return std::numeric_limits<T>::has_quiet_NaN && (x != x);
        return false;
    }
    static bool is_empty( const T& l, const T& u )
    {
        //return !(l <= u); // safety for partial orders
        return false;
    }
};
typedef boost::numeric::interval_lib::policies<boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_exact<double> >, my_checking_base<double> > interval_policy_exact_type;
typedef boost::numeric::interval_lib::policies<boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double> >, my_checking_base<double> > interval_policy_opp_type;
typedef boost::numeric::interval_lib::policies<boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_std<double> >, my_checking_base<double> > interval_policy_std_type;

typedef boost::numeric::interval_lib::default_policies<double>::type interval_policy_default_type;

template<class T, class P>
struct interval
{
    typedef boost::numeric::interval<T,P> type;
};
#endif // FEELPP_USES_BOOST_INTERVAL

template <class A, class B>
class SNumericalTraits
{
public:
};

//Specialization
template <class T>
class SNumericalTraits<T,T>
{
public:
    typedef T promote;
};

#define NUMERICAL_TRAITS(type1,type2,type3)      \
template <> class SNumericalTraits<type1,type2> { \
public:                                          \
    typedef type3 promote;                       \
};                                               \
template <> class SNumericalTraits<type2,type1> { \
public:                                          \
    typedef type3 promote;                       \
};

#if defined( FEELPP_USES_BOOST_INTERVAL )
#define INTERVAL_TRAITS( TYPE )                                                           \
typedef interval<double, interval_policy_## TYPE ##_type>::type interval_## TYPE ##_type; \
NUMERICAL_TRAITS(interval_## TYPE ##_type,double,interval_## TYPE ##_type)                \
NUMERICAL_TRAITS(interval_## TYPE ##_type,float,interval_## TYPE ##_type)                 \
NUMERICAL_TRAITS(interval_## TYPE ##_type,long,interval_## TYPE ##_type)                  \
NUMERICAL_TRAITS(interval_## TYPE ##_type,int,interval_## TYPE ##_type)

INTERVAL_TRAITS ( default );
INTERVAL_TRAITS ( exact );
INTERVAL_TRAITS ( opp );
INTERVAL_TRAITS ( std );

#undef INTERVAL_TRAITS
#endif // FEELPP_USES_BOOST_INTERVAL

NUMERICAL_TRAITS( double,std::complex<float>,std::complex<double> )
NUMERICAL_TRAITS( double,float,double )
NUMERICAL_TRAITS( double,long,double )
NUMERICAL_TRAITS( double,int,double )
NUMERICAL_TRAITS( float,long,float )
NUMERICAL_TRAITS( float,int,float )


#if defined( FEELPP_USES_BOOST_INTERVAL )
#define FEELPP_INTERVAL_DEBUG( TYPE )                                                        \
    SDebugStream& operator<<( SDebugStream& o, Feel::interval_## TYPE ##_type  const& e ); \
    SNDebugStream& operator<<( SNDebugStream& o, Feel::interval_## TYPE ##_type const& e );

FEELPP_INTERVAL_DEBUG( default );
FEELPP_INTERVAL_DEBUG( std );
FEELPP_INTERVAL_DEBUG( opp );
FEELPP_INTERVAL_DEBUG( exact );
#undef FEELPP_INTERVAL_DEBUG

#endif // FEELPP_USES_BOOST_INTERVAL

}
namespace std
{
#if defined( FEELPP_USES_BOOST_INTERVAL )
template<typename T, typename P>
std::ostream& operator<< ( std::ostream& __os, boost::numeric::interval<T,P> const& __i )
{
    __os << '[' << __i.lower() << ',' << __i.upper() << ']';
    return __os;
}
#endif // FEELPP_USES_BOOST_INTERVAL
}
#endif /* __ADTraits_H */


