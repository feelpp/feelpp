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
class ScalarBase
{
};
class VectorialBase
{
};
class Tensor2Base
{
};

//template<template<uint16_type> class PolySetType>
//struct is_scalar_field : std::is_base_of<ScalarBase,PolySetType<1> >::type {};

//template<template<uint16_type> class PolySetType>
//struct is_vector_field : std::is_base_of<ScalarBase,PolySetType<1> >::type {};

template <typename T>
struct is_scalar_polynomial : std::is_base_of<ScalarBase, T>::type
{
};

template <typename T>
struct is_vector_polynomial : std::is_base_of<VectorialBase, T>::type
{
};

template <typename T>
using is_vectorial_polynomial = is_vector_polynomial<T>;

template <typename T>
struct is_tensor2_polynomial : std::is_base_of<Tensor2Base, T>::type
{
};

template <typename T>
using is_matrix_polynomial = is_tensor2_polynomial<T>;

} // Feel
#endif /* __Traits_H */
