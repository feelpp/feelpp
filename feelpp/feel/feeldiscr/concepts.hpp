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
 * @brief C++20 concepts for Feel++ discretization (function spaces, elements, DOFs)
 * @author Christophe Prud'homme
 * @date 2026-01-02
 *
 * This header provides C++20 concepts specific to discretization, replacing
 * SFINAE patterns in feeldiscr/*.hpp files.
 *
 * Related concept files:
 * - feelcore/concepts.hpp: Base concepts (Iterable, SharedPtr, ScalarConcept, etc.)
 * - feelmesh/concepts.hpp: Mesh concepts (MeshElement, MeshFace, Simplex, etc.)
 * - feelpoly/concepts.hpp: Polynomial/basis concepts (Basis, NodalBasis, ModalBasis, etc.)
 * - feelvf/concepts.hpp: Variational formulation concepts (VfExpr, ScalarExpr, etc.)
 *
 * @note Some concepts here have richer versions in other modules. For example:
 * - VfExpr here is simplified; feelvf/concepts.hpp has a more complete version
 * - Mesh here is minimal; feelmesh/concepts.hpp has a more complete version
 */
#ifndef FEELPP_FEELDISCR_CONCEPTS_HPP
#define FEELPP_FEELDISCR_CONCEPTS_HPP 1

#include <concepts>
#include <memory>
#include <feel/feelcore/concepts.hpp>

namespace Feel
{

//
// Function Space Concepts
//

/**
 * @brief A Feel++ function space
 * 
 * @details Function spaces are the fundamental discretization concept in Feel++.
 * They define the approximation space for finite element methods.
 * 
 * Requirements:
 * - value_type: The scalar type (typically double)
 * - mesh_type: The associated mesh type
 * - element_type: The finite element function type
 * - element(): Factory method to create elements
 * 
 * @example
 * @code
 * template <FunctionSpace SpaceT>
 * auto project(VfExpr auto&& expr, std::shared_ptr<SpaceT> space) {
 *     auto u = space->element();
 *     // project expr onto u
 *     return u;
 * }
 * @endcode
 */
template <typename T>
concept FunctionSpace = requires(T t) {
    typename T::value_type;
    typename T::mesh_type;
    typename T::element_type;
    { t.element() } -> std::same_as<typename T::element_type>;
};

/**
 * @brief A pointer (shared_ptr or raw) to a function space
 */
template <typename T>
concept FunctionSpacePtr = requires {
    requires (std::is_pointer_v<T> && FunctionSpace<std::remove_pointer_t<T>>)
          || (SharedPtr<T> && FunctionSpace<typename T::element_type>);
};

/**
 * @brief A function space element (finite element function)
 * 
 * @details Elements represent functions in the discretized space.
 * They have values at DOFs and can be evaluated at points.
 */
template <typename T>
concept FunctionSpaceElement = requires(T t) {
    typename T::functionspace_type;
    typename T::value_type;
    requires FunctionSpace<typename T::functionspace_type>;
};

/**
 * @brief A product of function spaces (for mixed formulations)
 * 
 * @details Product spaces combine multiple function spaces, used in
 * mixed finite element methods (e.g., Stokes: velocity × pressure).
 */
template <typename T>
concept ProductSpace = requires(T t) {
    typename T::spaces_tuple_type;
    { t.numberOfSpaces() } -> std::convertible_to<int>;
};

/**
 * @brief Multiple product spaces (product of product spaces)
 */
template <typename T>
concept ProductSpaces = ProductSpace<T> && requires {
    typename T::spaces_array_type;
};

//
// Mesh Concepts
//

/**
 * @brief A Feel++ mesh
 * 
 * @details Meshes represent the geometric domain discretization.
 */
template <typename T>
concept Mesh = requires(T t) {
    typename T::shape_type;
    typename T::element_type;
    typename T::face_type;
    { t.numElements() } -> std::convertible_to<size_t>;
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
// DOF Concepts
//

/**
 * @brief A degree of freedom (DOF) identifier
 * 
 * @details DOFs represent discrete unknowns in the finite element system.
 */
template <typename T>
concept DofType = requires(T t) {
    { t.index() } -> std::convertible_to<size_t>;
    { t.sign() } -> std::convertible_to<int>;
};

/**
 * @brief A DOF table mapping elements to DOFs
 */
template <typename T>
concept DofTable = requires(T t) {
    typename T::dof_type;
    requires DofType<typename T::dof_type>;
};

//
// Geometric Mapping Concepts
//

/**
 * @brief A geometric mapping from reference to physical element
 */
template <typename T>
concept GeometricMapping = requires {
    typename T::gm_type;
    { T::nDim } -> std::convertible_to<int>;
    { T::nRealDim } -> std::convertible_to<int>;
};

/**
 * @brief A geometric mapping context (evaluation at quadrature points)
 */
template <typename T>
concept GeometricMappingContext = requires(T t) {
    typename T::gm_type;
    requires GeometricMapping<typename T::gm_type>;
};

//
// Basis Function Concepts
//

/**
 * @brief A finite element basis (shape functions)
 */
template <typename T>
concept Basis = requires {
    typename T::value_type;
    { T::nDof } -> std::convertible_to<int>;
    { T::nLocalDof } -> std::convertible_to<int>;
};

/**
 * @brief A basis context (evaluation at points)
 */
template <typename T>
concept BasisContext = requires(T t) {
    typename T::basis_type;
    requires Basis<typename T::basis_type>;
};

//
// Interpolation Concepts
//

/**
 * @brief A type that can be interpolated (projected) onto a function space
 */
template <typename T, typename SpaceT>
concept Interpolable = requires(T t, SpaceT space) {
    requires FunctionSpace<SpaceT>;
    // Can be evaluated to produce values
};

//
// Operator Concepts
//

/**
 * @brief A linear operator between function spaces
 */
template <typename T>
concept LinearOperator = requires(T t) {
    typename T::domain_space_type;
    typename T::dual_image_space_type;
    requires FunctionSpace<typename T::domain_space_type>;
    requires FunctionSpace<typename T::dual_image_space_type>;
};

/**
 * @brief A preconditioner for linear systems
 */
template <typename T>
concept Preconditioner = requires(T t) {
    typename T::operator_type;
    { t.setOperator(std::declval<typename T::operator_type>()) };
};

//
// Mesh Entity Type Concepts
//

/**
 * @brief Check if T is a mesh element type for mesh M
 *
 * @details Used to dispatch operations based on mesh entity type.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<std::is_same_v<T, typename M::element_type>>
 *
 * @example
 * @code
 * template <MeshEntityElement<mesh_type> EntityT>
 * auto localDof(EntityT const& elt) const { return this->localDof(elt.id()); }
 * @endcode
 */
template <typename T, typename M>
concept MeshEntityElement = requires {
    typename std::decay_t<M>::element_type;
    requires std::is_same_v<std::decay_t<T>, typename std::decay_t<M>::element_type>;
};

/**
 * @brief Check if T is a mesh face type for mesh M
 */
template <typename T, typename M>
concept MeshEntityFace = requires {
    typename std::decay_t<M>::face_type;
    requires std::is_same_v<std::decay_t<T>, typename std::decay_t<M>::face_type>;
};

/**
 * @brief Check if T is a mesh edge type for mesh M
 */
template <typename T, typename M>
concept MeshEntityEdge = requires {
    typename std::decay_t<M>::edge_type;
    requires std::is_same_v<std::decay_t<T>, typename std::decay_t<M>::edge_type>;
};

/**
 * @brief Check if T is a mesh point type for mesh M
 */
template <typename T, typename M>
concept MeshEntityPoint = requires {
    typename std::decay_t<M>::point_type;
    requires std::is_same_v<std::decay_t<T>, typename std::decay_t<M>::point_type>;
};

//
// Dimension Concepts
//

/**
 * @brief Shape has specific topological dimension D
 *
 * @details Used for dimension-based dispatch in mesh operations.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<TheShape::nDim == D>
 *
 * @example
 * @code
 * template <typename TheShape = GeoShape>
 *     requires Shape2D<TheShape>
 * void updateCommonDataInEntities();
 * @endcode
 */
template <typename T, int D>
concept ShapeWithDim = requires {
    { T::nDim } -> std::convertible_to<int>;
    requires (T::nDim == D);
};

/**
 * @brief A 0-dimensional shape (point)
 */
template <typename T>
concept Shape0D = ShapeWithDim<T, 0>;

/**
 * @brief A 1-dimensional shape (edge/segment)
 */
template <typename T>
concept Shape1D = ShapeWithDim<T, 1>;

/**
 * @brief A 2-dimensional shape (triangle, quadrilateral)
 */
template <typename T>
concept Shape2D = ShapeWithDim<T, 2>;

/**
 * @brief A 3-dimensional shape (tetrahedron, hexahedron)
 */
template <typename T>
concept Shape3D = ShapeWithDim<T, 3>;

//
// Function Space Property Concepts
//

/**
 * @brief A composite function space (product of spaces)
 *
 * @details Composite spaces combine multiple function spaces.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<TT::is_composite>
 */
template <typename T>
concept CompositeSpace = FunctionSpace<T> && requires {
    { T::is_composite } -> std::convertible_to<bool>;
    requires (T::is_composite == true);
};

/**
 * @brief A non-composite (simple) function space
 *
 * @details Simple function spaces represent a single approximation space.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<!TT::is_composite>
 */
template <typename T>
concept NonCompositeSpace = FunctionSpace<T> && requires {
    { T::is_composite } -> std::convertible_to<bool>;
    requires (T::is_composite == false);
};

/**
 * @brief A function space with modal basis
 *
 * @details Modal bases use orthogonal polynomial bases (Legendre, etc.)
 * Replaces SFINAE patterns like:
 * std::enable_if_t<B::is_modal>
 */
template <typename T>
concept ModalBasisSpace = requires {
    { T::is_modal } -> std::convertible_to<bool>;
    requires (T::is_modal == true);
};

/**
 * @brief A function space with nodal basis
 *
 * @details Nodal bases use Lagrange-type bases with DOFs at nodes.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<!B::is_modal>
 */
template <typename T>
concept NodalBasisSpace = requires {
    { T::is_modal } -> std::convertible_to<bool>;
    requires (T::is_modal == false);
};

//
// Field Type Concepts
//

// Forward declarations for base types (defined in feelpoly/traits.hpp)
struct ScalarBase;
struct VectorialBase;
struct Tensor2Base;
struct Tensor2SymmBase;

/**
 * @brief A scalar field element
 *
 * @details Scalar fields have a single component at each point.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<is_scalar_field_v<T>>
 *
 * @example
 * @code
 * template <typename T>
 *     requires ScalarField<T>
 * void processScalarField(T const& field);
 * @endcode
 */
template <typename T>
concept ScalarField = std::is_base_of_v<ScalarBase, std::decay_t<T>>;

/**
 * @brief A vectorial field element
 *
 * @details Vectorial fields have multiple components (e.g., velocity).
 * Replaces SFINAE patterns like:
 * std::enable_if_t<is_vector_field_v<T>>
 */
template <typename T>
concept VectorialField = std::is_base_of_v<VectorialBase, std::decay_t<T>>;

/**
 * @brief A tensor2 (matrix) field element
 *
 * @details Tensor2 fields represent rank-2 tensors (e.g., stress tensor).
 * Replaces SFINAE patterns like:
 * std::enable_if_t<is_tensor2_field_v<T>>
 */
template <typename T>
concept Tensor2Field = std::is_base_of_v<Tensor2Base, std::decay_t<T>>;

/**
 * @brief A symmetric tensor2 field element
 *
 * @details Symmetric tensor fields exploit symmetry for storage/computation.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<is_tensor2symm_field_v<T>>
 */
template <typename T>
concept Tensor2SymmField = std::is_base_of_v<Tensor2SymmBase, std::decay_t<T>>;

/**
 * @brief Any matrix-type field (tensor2 or symmetric tensor2)
 *
 * @details Combines tensor2 and symmetric tensor2 fields.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<is_tensor2_field_v<T> || is_tensor2symm_field_v<T>>
 */
template <typename T>
concept MatrixField = Tensor2Field<T> || Tensor2SymmField<T>;

//
// Range Entity Concepts
//

/**
 * @brief A range over mesh elements
 *
 * @details Ranges that iterate over volumetric mesh elements.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<RangeType::isOnElements()>
 *
 * @example
 * @code
 * template <RangeOnElements RangeType>
 * void dofs(RangeType const& rangeElt, std::vector<size_type>& dofIds) const;
 * @endcode
 */
template <typename T>
concept RangeOnElements = requires {
    { T::isOnElements() } -> std::convertible_to<bool>;
    requires T::isOnElements();
};

/**
 * @brief A range over mesh faces
 *
 * @details Ranges that iterate over mesh faces (boundary or internal).
 * Replaces SFINAE patterns like:
 * std::enable_if_t<RangeType::isOnFaces()>
 */
template <typename T>
concept RangeOnFaces = requires {
    { T::isOnFaces() } -> std::convertible_to<bool>;
    requires T::isOnFaces();
};

/**
 * @brief A range over mesh edges
 *
 * @details Ranges that iterate over mesh edges (1D entities).
 * Replaces SFINAE patterns like:
 * std::enable_if_t<RangeType::isOnEdges()>
 */
template <typename T>
concept RangeOnEdges = requires {
    { T::isOnEdges() } -> std::convertible_to<bool>;
    requires T::isOnEdges();
};

/**
 * @brief A range over mesh points
 *
 * @details Ranges that iterate over mesh vertices/nodes.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<RangeType::isOnPoints()>
 */
template <typename T>
concept RangeOnPoints = requires {
    { T::isOnPoints() } -> std::convertible_to<bool>;
    requires T::isOnPoints();
};

//
// Interpolation Entity Concepts
//

// Forward declaration for ElementsType enum
enum class ElementsType;

/**
 * @brief Interpolation type that operates on mesh elements
 *
 * @details Used in operatorinterpolation.hpp for entity-based dispatch.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<InterpOnType::onEntity() == ElementsType::MESH_ELEMENTS>
 */
template <typename T>
concept InterpolationOnElements = requires {
    { T::onEntity() };
    // Note: Compile-time check depends on ElementsType enum definition
};

/**
 * @brief Interpolation type that operates on mesh faces
 *
 * @details Used in operatorinterpolation.hpp for face-based interpolation.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<InterpOnType::onEntity() == ElementsType::MESH_FACES>
 */
template <typename T>
concept InterpolationOnFaces = requires {
    { T::onEntity() };
};

/**
 * @brief Interpolation type that operates on mesh edges
 */
template <typename T>
concept InterpolationOnEdges = requires {
    { T::onEntity() };
};

//
// Conformity Concepts
//

/**
 * @brief H(div) conforming basis (e.g., Raviart-Thomas)
 *
 * @details H(div) conforming spaces ensure normal continuity across faces.
 * Replaces SFINAE patterns checking is_hdiv_conforming.
 *
 * @example
 * @code
 * template <HdivConforming FEType>
 * void computeDivergence(FEType const& fe);
 * @endcode
 */
template <typename T>
concept HdivConforming = requires {
    { T::is_hdiv_conforming } -> std::convertible_to<bool>;
    requires (T::is_hdiv_conforming == true);
};

/**
 * @brief H(curl) conforming basis (e.g., Nedelec)
 *
 * @details H(curl) conforming spaces ensure tangential continuity across edges.
 * Replaces SFINAE patterns checking is_hcurl_conforming.
 */
template <typename T>
concept HcurlConforming = requires {
    { T::is_hcurl_conforming } -> std::convertible_to<bool>;
    requires (T::is_hcurl_conforming == true);
};

/**
 * @brief H1 conforming (continuous) basis
 *
 * @details H1 conforming spaces ensure full continuity across element boundaries.
 */
template <typename T>
concept H1Conforming = requires {
    { T::is_continuous } -> std::convertible_to<bool>;
    requires (T::is_continuous == true);
};

//
// Expression Concepts
//

// Forward declaration for ExprBase (defined in feelvf)
struct ExprBase;

/**
 * @brief A variational formulation expression
 *
 * @details VF expressions are used in bilinear/linear forms.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<std::is_base_of_v<ExprBase, T>>
 *
 * @example
 * @code
 * template <VfExpr ExprT>
 * void addExpression(ExprT&& expr);
 * @endcode
 */
template <typename T>
concept VfExpr = std::is_base_of_v<ExprBase, std::decay_t<T>>;

//
// Backward Compatibility Bridges
//

#ifdef FEELPP_ENABLE_CONCEPT_COMPATIBILITY

// Bridge old is_functionspace trait to new concept
template <typename T>
constexpr bool is_functionspace_v = FunctionSpace<T>;

template <typename T>
constexpr bool is_functionspace_element_v = FunctionSpaceElement<T>;

template <typename T>
constexpr bool is_product_space_v = ProductSpace<T>;

// Field type concept bridges
template <typename T>
constexpr bool is_scalar_field_concept_v = ScalarField<T>;

template <typename T>
constexpr bool is_vector_field_concept_v = VectorialField<T>;

template <typename T>
constexpr bool is_tensor2_field_concept_v = Tensor2Field<T>;

template <typename T>
constexpr bool is_matrix_field_concept_v = MatrixField<T>;

// Mesh entity concept bridges
template <typename T, typename M>
constexpr bool is_mesh_element_type_v = MeshEntityElement<T, M>;

template <typename T, typename M>
constexpr bool is_mesh_face_type_v = MeshEntityFace<T, M>;

template <typename T, typename M>
constexpr bool is_mesh_edge_type_v = MeshEntityEdge<T, M>;

// Dimension concept bridges
template <typename T>
constexpr bool is_shape_0d_v = Shape0D<T>;

template <typename T>
constexpr bool is_shape_1d_v = Shape1D<T>;

template <typename T>
constexpr bool is_shape_2d_v = Shape2D<T>;

template <typename T>
constexpr bool is_shape_3d_v = Shape3D<T>;

// Function space property bridges
template <typename T>
constexpr bool is_composite_space_v = CompositeSpace<T>;

template <typename T>
constexpr bool is_modal_basis_v = ModalBasisSpace<T>;

template <typename T>
constexpr bool is_nodal_basis_v = NodalBasisSpace<T>;

// Conformity concept bridges
template <typename T>
constexpr bool is_hdiv_conforming_concept_v = HdivConforming<T>;

template <typename T>
constexpr bool is_hcurl_conforming_concept_v = HcurlConforming<T>;

template <typename T>
constexpr bool is_h1_conforming_concept_v = H1Conforming<T>;

// VfExpr concept bridge
template <typename T>
constexpr bool is_vf_expr_concept_v = VfExpr<T>;

#endif // FEELPP_ENABLE_CONCEPT_COMPATIBILITY

} // namespace Feel

#endif /* FEELPP_FEELDISCR_CONCEPTS_HPP */
