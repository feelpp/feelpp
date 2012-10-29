/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
struct is_shared_ptr
        : mpl::false_
{
};

template <class T>
struct is_shared_ptr<boost::shared_ptr<T> >
        : mpl::true_
{
};

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

} // namespace Feel
#endif


