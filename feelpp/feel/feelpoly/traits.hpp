/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 19 Apr 2015
 
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
#ifndef __FEELPP_POLY_TRAITS_HPP
#define __FEELPP_POLY_TRAITS_HPP 1

#include <feel/feelmesh/traits.hpp>

namespace Feel
{
class ScalarBase {};
class VectorialBase {};
class Tensor2Base {};

template<typename T>
struct polynomial_order: std::integral_constant<int,T::order()> {};

template<typename T>
constexpr bool polynomial_order_v = polynomial_order<T>::value;


template<typename T>
struct is_linear_polynomial: std::integral_constant<bool,T::is_linear> {};
//struct is_linear_polynomial: std::integral_constant<bool,(polynomial_order_v<T> == 1)> {};

template<typename T>
constexpr bool is_linear_polynomial_v = is_linear_polynomial<T>::value;

//template<template<uint16_type> class PolySetType>
//struct is_scalar_field : std::is_base_of<ScalarBase,PolySetType<1> >::type {};

//template<template<uint16_type> class PolySetType>
//struct is_vector_field : std::is_base_of<ScalarBase,PolySetType<1> >::type {};
        
template<typename T>
struct is_scalar_polynomial : std::is_base_of<ScalarBase, T>::type {};

template<typename T>
struct is_vector_polynomial : std::is_base_of<VectorialBase, T>::type {};

template<typename T>
using is_vectorial_polynomial = is_vector_polynomial<T>;

template<typename T>
struct is_tensor2_polynomial : std::is_base_of<Tensor2Base, T>::type {};

template<typename T>
using is_matrix_polynomial = is_tensor2_polynomial<T>;

template<typename T>
using local_interpolant_t = typename T::local_interpolant_type;
template<typename T>
struct local_interpolant
{
    using type = local_interpolant_t<T>;
};
template<typename T>
using local_interpolants_t = typename T::local_interpolants_type;
template<typename T>
struct local_interpolants
{
    using type = local_interpolants_t<T>;
};

namespace detail {

inline constexpr uint16_type
symmetricIndex( uint16_type i, uint16_type j, uint16_type n)
{
    if ( j>=i )
        return  j + n*i - i*(i+1)/2;
    else
        return i + n*j - j*(j+1)/2;
}

}

} // Feel
#endif /* __Traits_H */
