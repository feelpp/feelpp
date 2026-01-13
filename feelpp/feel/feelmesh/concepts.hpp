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
 * @brief C++20 concepts for Feel++ meshes and geometric entities
 * @author Christophe Prud'homme
 * @date 2026-01-02
 */
#ifndef FEELPP_FEELMESH_CONCEPTS_HPP
#define FEELPP_FEELMESH_CONCEPTS_HPP 1

#include <concepts>
#include <memory>
#include <feel/feelcore/concepts.hpp>

namespace Feel
{

//
// Mesh Concepts
//

/**
 * @brief A Feel++ mesh
 * 
 * @details Meshes represent the geometric discretization of domains.
 * They contain elements, faces, edges, and points.
 */
template <typename T>
concept Mesh = requires(T t) {
    typename T::shape_type;
    typename T::element_type;
    typename T::face_type;
    typename T::point_type;
    typename T::value_type;
    { t.numElements() } -> std::convertible_to<size_t>;
    { t.numFaces() } -> std::convertible_to<size_t>;
    { t.numPoints() } -> std::convertible_to<size_t>;
};

/**
 * @brief A pointer to a mesh
 */
template <typename T>
concept MeshPtr = requires {
    requires (std::is_pointer_v<T> && Mesh<std::remove_pointer_t<T>>)
          || (SharedPtr<T> && Mesh<typename T::element_type>);
};

//
// Dimension Concepts
//

/**
 * @brief A 1D mesh
 */
template <typename T>
concept Mesh1D = Mesh<T> && requires {
    requires T::nDim == 1;
};

/**
 * @brief A 2D mesh
 */
template <typename T>
concept Mesh2D = Mesh<T> && requires {
    requires T::nDim == 2;
};

/**
 * @brief A 3D mesh
 */
template <typename T>
concept Mesh3D = Mesh<T> && requires {
    requires T::nDim == 3;
};

//
// Entity Concepts
//

/**
 * @brief A mesh element (cell)
 */
template <typename T>
concept MeshElement = requires(T t) {
    typename T::mesh_type;
    typename T::point_type;
    { t.id() } -> std::convertible_to<size_t>;
    { t.marker() };
    requires Mesh<typename T::mesh_type>;
};

/**
 * @brief A mesh face (boundary of element)
 */
template <typename T>
concept MeshFace = requires(T t) {
    typename T::mesh_type;
    typename T::element_type;
    { t.id() } -> std::convertible_to<size_t>;
    { t.isOnBoundary() } -> std::convertible_to<bool>;
    requires Mesh<typename T::mesh_type>;
};

/**
 * @brief A mesh edge
 */
template <typename T>
concept MeshEdge = requires(T t) {
    typename T::mesh_type;
    { t.id() } -> std::convertible_to<size_t>;
    requires Mesh<typename T::mesh_type>;
};

/**
 * @brief A mesh point (vertex)
 */
template <typename T>
concept MeshPoint = requires(T t) {
    typename T::value_type;
    { t.id() } -> std::convertible_to<size_t>;
    { t.node() };
};

//
// Range Concepts
//

/**
 * @brief A range of mesh elements
 */
template <typename T>
concept ElementRange = Iterable<T> && requires(T t) {
    requires MeshElement<std::decay_t<decltype(*t.begin())>>;
};

/**
 * @brief A range of mesh faces
 */
template <typename T>
concept FaceRange = Iterable<T> && requires(T t) {
    requires MeshFace<std::decay_t<decltype(*t.begin())>>;
};

/**
 * @brief A range of mesh edges
 */
template <typename T>
concept EdgeRange = Iterable<T> && requires(T t) {
    requires MeshEdge<std::decay_t<decltype(*t.begin())>>;
};

/**
 * @brief A range of mesh points
 */
template <typename T>
concept PointRange = Iterable<T> && requires(T t) {
    requires MeshPoint<std::decay_t<decltype(*t.begin())>>;
};

//
// Shape Concepts
//

/**
 * @brief A simplex shape (triangle, tetrahedron)
 */
template <typename T>
concept Simplex = requires {
    requires T::is_simplex;
};

/**
 * @brief A hypercube shape (quadrilateral, hexahedron)
 */
template <typename T>
concept Hypercube = requires {
    requires T::is_hypercube;
};

/**
 * @brief A triangle mesh element
 */
template <typename T>
concept Triangle = Simplex<typename T::shape_type> && requires {
    requires T::shape_type::nDim == 2;
    requires T::shape_type::nVertices == 3;
};

/**
 * @brief A tetrahedron mesh element
 */
template <typename T>
concept Tetrahedron = Simplex<typename T::shape_type> && requires {
    requires T::shape_type::nDim == 3;
    requires T::shape_type::nVertices == 4;
};

/**
 * @brief A quadrilateral mesh element
 */
template <typename T>
concept Quadrilateral = Hypercube<typename T::shape_type> && requires {
    requires T::shape_type::nDim == 2;
    requires T::shape_type::nVertices == 4;
};

/**
 * @brief A hexahedron mesh element
 */
template <typename T>
concept Hexahedron = Hypercube<typename T::shape_type> && requires {
    requires T::shape_type::nDim == 3;
    requires T::shape_type::nVertices == 8;
};

//
// Marker Concepts
//

/**
 * @brief An entity with a marker (for boundary conditions, materials)
 */
template <typename T>
concept HasMarker = requires(T t) {
    { t.marker() };
};

/**
 * @brief A marked range (filtered by marker)
 */
template <typename T>
concept MarkedRange = Iterable<T> && requires(T t) {
    requires HasMarker<std::decay_t<decltype(*t.begin())>>;
};

//
// Submesh Concepts
//

/**
 * @brief A submesh (portion of a parent mesh)
 */
template <typename T>
concept SubMesh = Mesh<T> && requires(T t) {
    typename T::parent_mesh_type;
    requires Mesh<typename T::parent_mesh_type>;
};

//
// Backward Compatibility Bridges
//

#ifdef FEELPP_ENABLE_CONCEPT_COMPATIBILITY

template <typename T>
constexpr bool is_mesh_v = Mesh<T>;

template <typename T>
constexpr bool is_mesh_element_v = MeshElement<T>;

template <typename T>
constexpr bool is_mesh_face_v = MeshFace<T>;

template <typename T>
constexpr bool is_simplex_v = Simplex<T>;

template <typename T>
constexpr bool is_hypercube_v = Hypercube<T>;

#endif // FEELPP_ENABLE_CONCEPT_COMPATIBILITY

} // namespace Feel

#endif /* FEELPP_FEELMESH_CONCEPTS_HPP */
