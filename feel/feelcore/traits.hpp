/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-28

  Copyright (C) 2007,2009 Université de Grenoble 1
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
   \file traits.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-28
 */
#ifndef FEELPP_TRAITS_HPP
#define FEELPP_TRAITS_HPP 1

#include <cmath>

#include <boost/rational.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/max_element.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/sizeof.hpp>

#include <feel/feelcore/feel.hpp>
#include <boost/hana/traits.hpp>

#include <feel/feelcore/typetraits.hpp>
#include <feel/feelcore/ublastraits.hpp>

namespace Feel
{
template <typename T1, typename T2, typename T3 = mpl::void_, typename T4 = mpl::void_, typename T5 = mpl::void_>
struct strongest_numeric_type
{
    typedef mpl::vector<T1, T2, T3, T4, T5> types;
    typedef typename mpl::max_element<mpl::transform_view< types,mpl::sizeof_<mpl::_1> > >::type iter;
    typedef typename mpl::deref<typename iter::base>::type type;
};

template <typename T1, typename T2>
struct strongest_numeric_type<T1,std::complex<T2> >
{
    typedef typename type_traits<T1>::real_type R1;
    typedef std::complex<typename strongest_numeric_type<R1,T2>::type > type;
};
template <typename T1, typename T2>
struct strongest_numeric_type<std::complex<T1>,T2 >
{
    typedef typename type_traits<T2>::real_type R2;
    typedef std::complex<typename strongest_numeric_type<T1,R2>::type > type;
};
template <typename T1, typename T2>
struct strongest_numeric_type<std::complex<T1>,std::complex<T2> >
{
    typedef std::complex<typename strongest_numeric_type<T1,T2>::type > type;
};


template <class T>
struct is_shared_ptr_t
        : mpl::false_
{
};

template <class T>
struct is_shared_ptr_t<boost::shared_ptr<T> >
        : mpl::true_
{
};
template <class T>
struct is_shared_ptr_t<std::shared_ptr<T> >
    : mpl::true_
{
};

template <class T>
using is_shared_ptr_type = is_shared_ptr_t<T>;

#if 0
template <class T>
struct remove_shared_ptr
{
    typedef T type;
};

template <class T>
struct remove_shared_ptr<boost::shared_ptr<T> >
{
    typedef T type;
};
#endif

template<class T>
constexpr bool is_shared_ptr_v = is_shared_ptr_t<T>::value;

template<typename T>
struct is_ptr_or_shared_ptr_t : mpl::or_<is_shared_ptr_t<T>, boost::is_pointer<T> >::type {};

template<class T>
constexpr bool is_ptr_or_shared_ptr_v = is_ptr_or_shared_ptr_t<T>::value;

template<typename T>
using remove_shared_ptr_type = typename mpl::if_<is_shared_ptr_t<std::remove_pointer_t<T>>, mpl::identity<typename std::remove_pointer_t<T>::element_type>, mpl::identity<T>>::type::type;

template<typename T>
using decay_type = typename std::decay<remove_shared_ptr_type<typename std::decay<T>::type>>::type;

/**
 * @return true if \p e is a shared ptr, false otherwise
 */
template<typename T>
inline constexpr bool 
is_shared_ptr( T&& e )
{
    return is_shared_ptr_v<T>;
}

template<typename T>
inline constexpr bool 
is_ptr( T&& e )
{
    return is_ptr_or_shared_ptr_v<T>;
}

/**
 * @return the pointee object is e is a shared pointer, e otherwise
 */
template<typename T>
inline constexpr decltype(auto) 
remove_shared_ptr( T&& e )
{
    return hana::if_( hana::bool_<is_shared_ptr_v<T>>{},
                     []( auto&& x ) { return *x; },
                     []( auto&& x ) { return x; } )( std::forward<T>(e) );

}

template<typename T>
inline constexpr decltype(auto) 
remove_ptr( T&& e )
{
    return hana::if_( hana::bool_<is_ptr_or_shared_ptr_v<std::remove_reference_t<T>>>{},
                      []( auto&& x ) { return *std::forward<decltype(x)>(x); },
                      []( auto&& x ) { return std::forward<decltype(x)>(x); } )( std::forward<T>(e) );

}

template<typename T>
inline constexpr T const& //decltype(auto) 
remove_ptr( T const& e )
{
    return hana::if_( hana::bool_<is_ptr_or_shared_ptr_v<T>>{},
                      []( decay_type<T> const& x ) { return *x; },
                      []( decay_type<T> const& x ) { return x; } )( e );

}

template<typename T>
inline constexpr T const& //decltype(auto) 
remove_ptr( T const* e )
{
    return *e;

}

} // namespace Feel
#endif
