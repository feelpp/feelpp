/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file concepts.hpp
 * @brief C++20 concepts for mathematical mesh-entity views.
 */
#if !defined( FEELPP_FEELMESH_CONCEPTS_HPP )
#define FEELPP_FEELMESH_CONCEPTS_HPP 1

#include <concepts>
#include <ranges>
#include <type_traits>

namespace Feel
{

/**
 * @brief Contract for a non-owning mathematical geometric-point view.
 * @tparam T candidate type
 */
template <typename T>
concept GeometricPointViewType =
    std::semiregular<std::remove_cvref_t<T>> &&
    requires( std::remove_cvref_t<T> const& view )
    {
        typename std::remove_cvref_t<T>::mesh_type;
        typename std::remove_cvref_t<T>::point_type;
        typename std::remove_cvref_t<T>::value_type;
        typename std::remove_cvref_t<T>::index_type;
        { view.id() } -> std::convertible_to<
            typename std::remove_cvref_t<T>::index_type>;
        { view.localIndex() } -> std::convertible_to<unsigned int>;
        view.coordinates();
        view.point();
    };

/**
 * @brief Contract for a non-owning edge, face, or cell view.
 * @tparam T candidate type
 */
template <typename T>
concept GeometricEntityViewType =
    std::semiregular<std::remove_cvref_t<T>> &&
    requires( std::remove_cvref_t<T> const& view )
    {
        typename std::remove_cvref_t<T>::mesh_type;
        typename std::remove_cvref_t<T>::entity_type;
        typename std::remove_cvref_t<T>::value_type;
        typename std::remove_cvref_t<T>::index_type;
        { std::remove_cvref_t<T>::nDim } ->
            std::convertible_to<unsigned int>;
        { std::remove_cvref_t<T>::nRealDim } ->
            std::convertible_to<unsigned int>;
        { view.id() } -> std::convertible_to<
            typename std::remove_cvref_t<T>::index_type>;
        { view.pointCount() } -> std::convertible_to<unsigned int>;
        { view.pointIds() } -> std::ranges::random_access_range;
        { view.pointLocalIndices() } -> std::ranges::random_access_range;
        { view.point( 0 ) } -> GeometricPointViewType;
        view.entity();
    };

/**
 * @brief Contract refining an entity view as a top-dimensional mesh cell.
 * @tparam T candidate type
 */
template <typename T>
concept GeometricCellViewType =
    GeometricEntityViewType<T> &&
    ( std::remove_cvref_t<T>::nDim ==
      std::remove_cvref_t<T>::mesh_type::nDim );

/**
 * @brief Contract refining an entity view as a codimension-one entity.
 * @tparam T candidate type
 */
template <typename T>
concept GeometricFaceViewType =
    GeometricEntityViewType<T> &&
    ( std::remove_cvref_t<T>::nDim + 1 ==
      std::remove_cvref_t<T>::mesh_type::nDim );

/**
 * @brief Contract refining an entity view as a one-dimensional entity.
 * @tparam T candidate type
 *
 * In a two-dimensional mesh an entity can satisfy both
 * `GeometricEdgeViewType` and `GeometricFaceViewType`, matching its
 * mathematical roles.
 */
template <typename T>
concept GeometricEdgeViewType =
    GeometricEntityViewType<T> &&
    ( std::remove_cvref_t<T>::nDim == 1 );

} // namespace Feel

#endif // FEELPP_FEELMESH_CONCEPTS_HPP
