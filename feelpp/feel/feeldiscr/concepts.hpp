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
 * SFINAE patterns in feeldiscr headers.
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
#include <type_traits>
#include <feel/feelcore/concepts.hpp>
#include <feel/feeldiscr/doflayout.hpp>

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
concept FunctionSpaceConcept = requires(T t) {
    typename T::value_type;
    typename T::mesh_type;
    typename T::element_type;
    { t.element() } -> std::same_as<typename T::element_type>;
};

/**
 * @brief A pointer (shared_ptr or raw) to a function space
 */
template <typename T>
concept FunctionSpacePtrConcept = requires {
    requires (std::is_pointer_v<T> && FunctionSpaceConcept<std::remove_pointer_t<T>>)
          || (SharedPtr<T> && FunctionSpaceConcept<typename T::element_type>);
};

/**
 * @brief A function space element (finite element function)
 *
 * @details Elements represent functions in the discretized space.
 * They have values at DOFs and can be evaluated at points.
 */
template <typename T>
concept FunctionSpaceElementConcept = requires(T t) {
    typename T::functionspace_type;
    typename T::value_type;
    requires FunctionSpaceConcept<typename T::functionspace_type>;
};

/**
 * @brief A product of function spaces (for mixed formulations)
 *
 * @details Product spaces combine multiple function spaces, used in
 * mixed finite element methods (e.g., Stokes: velocity × pressure).
 */
template <typename T>
concept ProductSpaceConcept = requires(T t) {
    typename std::remove_cvref_t<T>::functionspace_type;
    typename std::remove_cvref_t<T>::value_type;
    typename std::remove_cvref_t<T>::element_type;
    { t.numberOfSpaces() } -> std::convertible_to<int>;
    { t.nDof() } -> std::convertible_to<std::size_t>;
    { t.nLocalDof() } -> std::convertible_to<std::size_t>;
    { t.nDofStart(0) } -> std::convertible_to<std::size_t>;
    { t.nLocalDofStart(0) } -> std::convertible_to<std::size_t>;
    { t.blockDofStart(0) } -> std::convertible_to<std::size_t>;
    { t.blockLocalDofStart(0) } -> std::convertible_to<std::size_t>;
    t.blockMapPtr(0);
    { t.element() } -> std::same_as<typename std::remove_cvref_t<T>::element_type>;
};

/**
 * @brief Multiple product spaces (product of product spaces)
 */
template <typename T>
concept ProductSpacesConcept = ProductSpaceConcept<T> && requires {
    typename std::remove_cvref_t<T>::tuple_spaces_type;
} && requires(T t) {
    t.tupleSpaces();
    t.template space<0>();
};

/**
 * @brief An element owned by a product-space object.
 */
template <typename T>
concept BlockElementConcept = requires(T t) {
    typename std::remove_cvref_t<T>::functionspace_type;
    typename std::remove_cvref_t<T>::value_type;
    requires ProductSpaceConcept<typename std::remove_cvref_t<T>::functionspace_type>;
    t.functionSpace();
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
concept MeshConcept = requires(T t) {
    typename T::shape_type;
    typename T::element_type;
    typename T::face_type;
    { t.numElements() } -> std::convertible_to<size_t>;
};

/**
 * @brief A pointer to a mesh
 */
template <typename T>
concept MeshPtrConcept = requires {
    requires (std::is_pointer_v<T> && MeshConcept<std::remove_pointer_t<T>>)
          || (SharedPtr<T> && MeshConcept<typename T::element_type>);
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
concept DofTypeConcept = requires(T t) {
    { t.index() } -> std::convertible_to<size_t>;
    { t.sign() } -> std::convertible_to<int>;
};

/**
 * @brief A DOF table mapping elements to DOFs
 */
template <typename T>
concept DofTableConcept = requires(T t) {
    typename T::dof_type;
    requires DofTypeConcept<typename T::dof_type>;
};

//
// Geometric Mapping Concepts
//

/**
 * @brief A geometric mapping from reference to physical element
 */
template <typename T>
concept GeoMapConcept = requires {
    typename T::gm_type;
    { T::nDim } -> std::convertible_to<int>;
    { T::nRealDim } -> std::convertible_to<int>;
};

/**
 * @brief A geometric mapping context (evaluation at quadrature points)
 */
template <typename T>
concept GeoMapContextConcept = requires(T t) {
    typename T::gm_type;
    requires GeoMapConcept<typename T::gm_type>;
};

//
// Basis Function Concepts
//

/**
 * @brief A finite element basis (shape functions)
 */
template <typename T>
concept FEBasisConcept = requires {
    typename T::value_type;
    { T::nDof } -> std::convertible_to<int>;
    { T::nLocalDof } -> std::convertible_to<int>;
};

/**
 * @brief A basis context (evaluation at points)
 */
template <typename T>
concept FEBasisContextConcept = requires(T t) {
    typename T::basis_type;
    requires FEBasisConcept<typename T::basis_type>;
};

//
// Interpolation Concepts
//

/**
 * @brief A type that can be interpolated (projected) onto a function space
 */
template <typename T, typename SpaceT>
concept InterpolableConcept = requires(T t, SpaceT space) {
    requires FunctionSpaceConcept<SpaceT>;
    // Can be evaluated to produce values
};

//
// Operator Concepts
//

/**
 * @brief A linear operator between function spaces
 */
template <typename T>
concept LinearOperatorConcept = requires(T t) {
    typename T::domain_space_type;
    typename T::dual_image_space_type;
    requires FunctionSpaceConcept<typename T::domain_space_type>;
    requires FunctionSpaceConcept<typename T::dual_image_space_type>;
};

/**
 * @brief A preconditioner for linear systems
 * @note See feelvf/concepts.hpp for PreconditionerConcept with backend_type
 */
template <typename T>
concept DiscrPreconditionerConcept = requires(T t) {
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
concept CompositeSpaceConcept = FunctionSpaceConcept<T> && requires {
    { T::is_composite } -> std::convertible_to<bool>;
    requires (T::is_composite == true);
};

/**
 * @brief A legacy composite implemented by multi-basis FunctionSpace internals.
 *
 * New mixed-space code should avoid this concept and use product-backed
 * composition instead.
 */
template <typename T>
concept LegacyCompositeFunctionSpaceConcept =
    CompositeSpaceConcept<std::remove_cvref_t<T>> && requires {
        { std::remove_cvref_t<T>::uses_internal_composite } -> std::convertible_to<bool>;
        { std::remove_cvref_t<T>::is_legacy_composite } -> std::convertible_to<bool>;
        requires (std::remove_cvref_t<T>::uses_internal_composite == true);
        requires (std::remove_cvref_t<T>::is_legacy_composite == true);
    };

/**
 * @brief A composite facade backed by product.hpp/ProductFunctionSpaces.
 *
 * This is the preferred compatibility shape for mixed spaces that still need
 * FunctionSpace-like accessors while storing composition outside FunctionSpace.
 */
template <typename T>
concept ProductBackedCompositeSpaceConcept =
    CompositeSpaceConcept<std::remove_cvref_t<T>> &&
    ProductSpaceConcept<std::remove_cvref_t<T>> &&
    requires {
        { std::remove_cvref_t<T>::uses_internal_composite } -> std::convertible_to<bool>;
        { std::remove_cvref_t<T>::is_product_backed_composite } -> std::convertible_to<bool>;
        requires (std::remove_cvref_t<T>::uses_internal_composite == false);
        requires (std::remove_cvref_t<T>::is_product_backed_composite == true);
    };

/**
 * @brief A non-composite (simple) function space
 *
 * @details Simple function spaces represent a single approximation space.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<!TT::is_composite>
 */
template <typename T>
concept NonCompositeSpaceConcept = FunctionSpaceConcept<T> && requires {
    { T::is_composite } -> std::convertible_to<bool>;
    requires (T::is_composite == false);
};

/**
 * @brief An explicit mortar function space.
 *
 * @details Mortar spaces are non-composite spaces with a mortar construction
 * policy.  New code should depend on this concept instead of detecting raw
 * FunctionSpace<..., mortars<Mortar>> instantiations.
 */
template <typename T>
concept MortarFunctionSpaceConcept = NonCompositeSpaceConcept<T> && requires {
    typename std::remove_cvref_t<T>::mortar_policy_type;
    typename std::remove_cvref_t<T>::mortar_0_type;
    { std::remove_cvref_t<T>::is_mortar } -> std::convertible_to<bool>;
    { std::remove_cvref_t<T>::is_explicit_mortar_space } -> std::convertible_to<bool>;
    requires (std::remove_cvref_t<T>::is_mortar == true);
    requires (std::remove_cvref_t<T>::is_explicit_mortar_space == true);
};

/**
 * @brief A function space with modal basis
 *
 * @details Modal bases use orthogonal polynomial bases (Legendre, etc.)
 * Replaces SFINAE patterns like:
 * std::enable_if_t<B::is_modal>
 */
template <typename T>
concept ModalBasisSpaceConcept = requires {
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
concept NodalBasisSpaceConcept = requires {
    { T::is_modal } -> std::convertible_to<bool>;
    requires (T::is_modal == false);
};

//
// Field Type Concepts
//

// Forward declarations for base types (defined in feelpoly/traits.hpp and policy.hpp)
class ScalarBase;
class VectorialBase;
class Tensor2Base;
struct Tensor2SymmBase; // struct in policy.hpp

/**
 * @brief A scalar field element
 *
 * @details Scalar fields have a single component at each point.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<is_scalar_field_v<T>>
 *
 * @note See feelpoly/concepts.hpp for ScalarFieldConcept
 *
 * @example
 * @code
 * template <typename T>
 *     requires DiscrScalarFieldConcept<T>
 * void processScalarField(T const& field);
 * @endcode
 */
template <typename T>
concept DiscrScalarFieldConcept = std::is_base_of_v<ScalarBase, std::decay_t<T>>;

/**
 * @brief A vectorial field element
 *
 * @details Vectorial fields have multiple components (e.g., velocity).
 * Replaces SFINAE patterns like:
 * std::enable_if_t<is_vector_field_v<T>>
 */
template <typename T>
concept VectorialFieldConcept = std::is_base_of_v<VectorialBase, std::decay_t<T>>;

/**
 * @brief A tensor2 (matrix) field element
 *
 * @details Tensor2 fields represent rank-2 tensors (e.g., stress tensor).
 * Replaces SFINAE patterns like:
 * std::enable_if_t<is_tensor2_field_v<T>>
 *
 * @note See feelpoly/concepts.hpp for Tensor2FieldConcept
 */
template <typename T>
concept DiscrTensor2FieldConcept = std::is_base_of_v<Tensor2Base, std::decay_t<T>>;

/**
 * @brief A symmetric tensor2 field element
 *
 * @details Symmetric tensor fields exploit symmetry for storage/computation.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<is_tensor2symm_field_v<T>>
 */
template <typename T>
concept Tensor2SymmFieldConcept = std::is_base_of_v<Tensor2SymmBase, std::decay_t<T>>;

/**
 * @brief Any matrix-type field (tensor2 or symmetric tensor2)
 *
 * @details Combines tensor2 and symmetric tensor2 fields.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<is_tensor2_field_v<T> || is_tensor2symm_field_v<T>>
 */
template <typename T>
concept MatrixFieldConcept = DiscrTensor2FieldConcept<T> || Tensor2SymmFieldConcept<T>;

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

// Note: ElementsType enum is defined in feelmesh/enums.hpp (unscoped enum)

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
concept HdivConforming =
    HDivFiniteElement<T> ||
    requires {
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
concept HcurlConforming =
    HCurlFiniteElement<T> ||
    requires {
        { T::is_hcurl_conforming } -> std::convertible_to<bool>;
        requires (T::is_hcurl_conforming == true);
    };

/**
 * @brief H1 conforming (continuous) basis
 *
 * @details H1 conforming spaces ensure full continuity across element boundaries.
 */
template <typename T>
concept H1Conforming =
    H1FiniteElement<T> ||
    requires {
        { T::is_continuous } -> std::convertible_to<bool>;
        requires (T::is_continuous == true);
    };

//
// Expression Concepts
//

// Forward declaration for ExprBase (defined in feelvf/exprbase.hpp)
class ExprBase;

/**
 * @brief A variational formulation expression (simplified version)
 *
 * @details VF expressions are used in bilinear/linear forms.
 * Replaces SFINAE patterns like:
 * std::enable_if_t<std::is_base_of_v<ExprBase, T>>
 *
 * @note See feelvf/concepts.hpp for VfExprConcept (primary definition)
 *
 * @example
 * @code
 * template <DiscrVfExprConcept ExprT>
 * void addExpression(ExprT&& expr);
 * @endcode
 */
template <typename T>
concept DiscrVfExprConcept = std::is_base_of_v<ExprBase, std::decay_t<T>>;

//
// Backward Compatibility Bridges
//

#ifdef FEELPP_ENABLE_CONCEPT_COMPATIBILITY

// Bridge old is_functionspace trait to new concept
template <typename T>
constexpr bool is_functionspace_v = FunctionSpaceConcept<T>;

template <typename T>
constexpr bool is_functionspace_element_v = FunctionSpaceElementConcept<T>;

template <typename T>
constexpr bool is_product_space_v = ProductSpaceConcept<T>;

// Field type concept bridges
template <typename T>
constexpr bool is_scalar_field_concept_v = DiscrScalarFieldConcept<T>;

template <typename T>
constexpr bool is_vector_field_concept_v = VectorialFieldConcept<T>;

template <typename T>
constexpr bool is_tensor2_field_concept_v = DiscrTensor2FieldConcept<T>;

template <typename T>
constexpr bool is_matrix_field_concept_v = MatrixFieldConcept<T>;

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
constexpr bool is_composite_space_v = CompositeSpaceConcept<T>;

template <typename T>
constexpr bool is_mortar_space_v = MortarFunctionSpaceConcept<T>;

template <typename T>
constexpr bool is_modal_basis_v = ModalBasisSpaceConcept<T>;

template <typename T>
constexpr bool is_nodal_basis_v = NodalBasisSpaceConcept<T>;

// Conformity concept bridges
template <typename T>
constexpr bool is_hdiv_conforming_concept_v = HdivConforming<T>;

template <typename T>
constexpr bool is_hcurl_conforming_concept_v = HcurlConforming<T>;

template <typename T>
constexpr bool is_h1_conforming_concept_v = H1Conforming<T>;

// VfExpr concept bridge
template <typename T>
constexpr bool is_vf_expr_concept_v = DiscrVfExprConcept<T>;

#endif // FEELPP_ENABLE_CONCEPT_COMPATIBILITY

} // namespace Feel

#endif /* FEELPP_FEELDISCR_CONCEPTS_HPP */
