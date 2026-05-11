/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2026-01-02

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
/**
 * @file concepts.hpp
 * @brief C++20 concepts for Feel++ view factor computations
 * @author Christophe Prud'homme
 * @date 2026-01-02
 */
#ifndef FEELPP_FEELVIEWFACTOR_CONCEPTS_HPP
#define FEELPP_FEELVIEWFACTOR_CONCEPTS_HPP 1

#include <concepts>
#include <memory>
#include <string>
#include <vector>
#include <feel/feeldiscr/traits.hpp>

namespace Feel
{

/**
 * @brief Concept for meshes compatible with view factor computation
 *
 * Requirements:
 * - Must have face_type for boundary iteration
 * - Must have nDim (2 or 3 only)
 * - Must support marker queries
 * - Must provide geometric mapping (gm())
 */
template<typename M>
concept ViewFactorMesh = requires(M m, std::string marker) {
    typename M::face_type;
    typename M::element_type;
    { M::nDim } -> std::convertible_to<int>;
    { m.hasMarker(marker) } -> std::convertible_to<bool>;
    { m.gm() };
    { m.numElements() } -> std::convertible_to<size_t>;
} && (M::nDim == 2 || M::nDim == 3);

/**
 * @brief Concept for view factor computation classes
 *
 * All view factor implementations must satisfy this concept.
 */
template<typename T>
concept ViewFactorComputable = requires(T t, bool elementwise) {
    typename T::value_type;
    typename T::mesh_t;
    { t.compute(elementwise) } -> std::same_as<void>;
    { t.viewFactors() };
    { t.areas() };
};

/**
 * @brief Concept for view factor base class compatibility
 */
template<typename T>
concept ViewFactorDerived = ViewFactorComputable<T> && requires(T t, unsigned int i, unsigned int j) {
    { t.devReciprocity(i, j) } -> std::convertible_to<typename T::value_type>;
    { t.maxDevReciprocity() } -> std::convertible_to<typename T::value_type>;
};

/**
 * @brief Concept for ray-traceable geometry
 *
 * Meshes that support BVH-based ray tracing.
 */
template<typename M>
concept RayTraceableMesh = ViewFactorMesh<M> && requires(M m) {
    typename trace_mesh_t<M>;
};

/**
 * @brief Concept for JSON configuration objects
 */
template<typename J>
concept ViewFactorConfig = requires(J j) {
    { j["viewfactor"]["markers"] };
    { j["viewfactor"]["quadrature_order"] };
};

} // namespace Feel

#endif /* FEELPP_FEELVIEWFACTOR_CONCEPTS_HPP */
