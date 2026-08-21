/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomm@feelpp.org>
 Date: 2026-03-26

 Copyright (C) 2026 Feel++ Consortium

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
#ifndef FEELPP_FEELALG_PRODUCTSPACECONCEPTS_HPP
#define FEELPP_FEELALG_PRODUCTSPACECONCEPTS_HPP 1

#include <concepts>
#include <vector>

#include <boost/hana/concept/foldable.hpp>
#include <boost/hana/integral_constant.hpp>

#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/condenser.hpp>
#include <feel/feelalg/products.hpp>
#include <feel/feeldiscr/traits.hpp>

namespace Feel
{

template <typename T>
concept StaticProductSpacesType =
    std::is_base_of_v<ProductSpacesBase, decay_type<T>> &&
    requires( decay_type<T>& ps )
    {
        typename decay_type<T>::tuple_spaces_type;
        { ps.numberOfSpaces() } -> std::convertible_to<int>;
        ps.tupleSpaces();
    };

template <typename T>
concept DynamicProductSpaceType =
    std::is_base_of_v<ProductSpaceBase, decay_type<T>> &&
    requires( decay_type<T>& ps, int i )
    {
        typename decay_type<T>::underlying_functionspace_type;
        { ps.numberOfSpaces() } -> std::convertible_to<int>;
        ps[i];
    };

template <typename T>
concept FoldableProductSpacesType =
    StaticProductSpacesType<T> &&
    requires { typename decay_type<T>::tuple_spaces_type; } &&
    boost::hana::Foldable<typename decay_type<T>::tuple_spaces_type>::value;

template <typename T>
concept CondenserTag = is_condenser_v<decay_type<T>>;

template <typename T>
concept StokesCondenserTag = std::same_as<decay_type<T>, condenser_stokes>;

template <typename T>
concept Sb9CondenserTag = std::same_as<decay_type<T>, condenser_sb9>;

template <typename T, int N>
concept ProductElementNSpaces =
    requires( decay_type<T>& e )
    {
        decay_type<T>::nspaces;
        e( boost::hana::int_c<0> );
    } &&
    ( decay_type<T>::nspaces == N );

template <typename T, int N>
concept ProductElementAtLeastNSpaces =
    requires( decay_type<T>& e )
    {
        decay_type<T>::nspaces;
        e( boost::hana::int_c<0> );
    } &&
    ( decay_type<T>::nspaces >= N );

template <typename T>
concept Tensor2SymmFieldType = is_tensor2symm_field_v<decay_type<T>>;

template <typename T>
concept LocalInterpolantFieldType =
    requires( decay_type<T>& e,
              typename decay_type<T>::index_type K,
              typename decay_type<T>::local_interpolant_type loc )
    {
        typename decay_type<T>::index_type;
        typename decay_type<T>::local_interpolant_type;
        e.dof();
        e.assignE( K, loc );
        e.element( std::vector<typename decay_type<T>::index_type>{ K }, loc );
    };

template <typename T>
concept Sb9CondensableProductElement =
    ProductElementNSpaces<T,2> &&
    requires( decay_type<T>& e )
    {
        requires LocalInterpolantFieldType<decltype( e( boost::hana::int_c<0> ) )>;
        requires LocalInterpolantFieldType<decltype( e( boost::hana::int_c<1> ) )>;
    };

} // namespace Feel

#endif
