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
#ifndef FEELPP_DISCRTRAITS_HPP
#define FEELPP_DISCRTRAITS_HPP 1

#include <feel/feelcore/unwrapptr.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelpoly/traits.hpp>
#include <feel/feelpoly/policy.hpp>
#include <feel/feelmesh/meshbase.hpp>
#include <feel/feeldiscr/functionspacebase.hpp>

namespace Feel {

/**
 * \ingroup Traits
 * @{
 */
/**
 * if \p T has base class \p MeshBase<> (hense if it is a function space)
 * then provides the member constant value equal to true, false otherwise
 */
template<typename MeshType>
using is_mesh = typename std::is_base_of<MeshBase<>,decay_type<MeshType>>::type;

/**
 * provides the mesh  type
 * if \p MeshType is a shared_ptr of a Mesh then provides the mesh type
 * \note it checks that the \p Mesh is indeed a mesh type and return void if it is not the case.
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
using mesh_t = decay_type<MeshType>;


/**
 * helper variable template for is_mesh
 */
template<typename MeshType>
constexpr bool is_mesh_v = is_mesh<MeshType>::value;

//!
//! @return the topogical dimension of the mesh \p m
//!
template <typename MeshType, typename = std::enable_if_t<is_mesh_v<MeshType>>>
inline constexpr int topodim( std::shared_ptr<MeshType> const& m )
{
    return MeshType::nDim;
}
//!
//! @return the real dimension in which the mesh is defined
//!
template <typename MeshType, typename = std::enable_if_t<is_mesh_v<MeshType>>>
inline constexpr int realdim( std::shared_ptr<MeshType> const& m )
{
    return MeshType::nRealDim;
}

template<typename MeshType>
constexpr bool is_simplex_mesh_v = is_simplex_v<typename MeshType::shape_type>;

/**
 * if \p T has base class \p ScalarBase then @return the member constant value equal
 * to true, false otherwise
 */
template<typename T>
using is_scalar_field = typename std::is_base_of<ScalarBase, T>::type;

/**
 * helper variable template for is_scalar_field
 */
template<typename T>
constexpr bool is_scalar_field_v = is_scalar_field<T>::value;

/**
 * if \p T has base class \p VectorialBase then @return the member constant value equal
 * to true, false otherwise
 */
template<typename T>
using is_vector_field = typename std::is_base_of<VectorialBase, T>::type;

/**
 * helper variable template for is_vector_field
 */
template<typename T>
constexpr bool is_vector_field_v = is_vector_field<T>::value;

/**
 * if \p T has base class \p Tensor2Base then @return the member constant value equal
 * to true, false otherwise
 */
template<typename T>
using is_tensor2_field  =typename std::is_base_of<Tensor2Base, T>::type;

/**
 * helper variable template for is_tensor2_field
 */
template<typename T>
constexpr bool is_tensor2_field_v = is_tensor2_field<T>::value;

/**
 * if \p T has base class \p Tensor2SymmBase then @return the member constant value equal
 * to true, false otherwise
 */
template<typename T>
using is_tensor2symm_field  =typename std::is_base_of<Tensor2SymmBase, T>::type;
using tensor2symm_true  = std::integral_constant<bool,true>;
using tensor2symm_false  = std::integral_constant<bool,false>;

/**
 * helper variable template for is_tensor2symm_field
 */
template<typename T>
constexpr bool is_tensor2symm_field_v = is_tensor2symm_field<T>::value;

/**
 * helper variable template, 
 * @return true of field is a matrix field field, false otherwise
 */
template<typename T>
constexpr bool is_matrix_field_v  = is_tensor2_field_v<T> || is_tensor2symm_field_v<T>;

/**
 * if \p T has base class \p FunctionSpaceBase (hense if it is a function space)
 * then provides the member constant value equal to true, false otherwise
 */
template<typename FuncSpaceType>
using is_functionspace = typename std::is_base_of<FunctionSpaceBase,decay_type<FuncSpaceType>>::type;

/**
 * helper variable template for is_functionspace
 */
template<typename FuncSpaceType>
constexpr bool is_functionspace_v = is_functionspace<FuncSpaceType>::value;

/**
 * provides the function space type
 * if \p FESpace is a shared_ptr of a function space then provides the function space type
 * \note it checks that FESpace is indeed a functionspace type and return void if it is not the case.
 */
template<typename FESpace>
using functionspace_type = typename mpl::if_<is_functionspace<decay_type<FESpace> >,
                                             mpl::identity<decay_type<FESpace>>,
                                             mpl::identity<void> >::type::type;


/**
 * if \p T has base class \p ElementBase (hense if it is an element of a function space)
 * then provides the member constant value equal to true, false otherwise
 */
template<typename ElementType>
using is_functionspace_element = typename std::is_base_of<FunctionSpaceBase::ElementBase,ElementType>::type;

/**
 * helper variable template for is_functionspace_element
 */
template<typename ElementType>
constexpr bool is_functionspace_element_v = is_functionspace_element<ElementType>::value;


/**
 * provides the function space type
 * if \p FESpace is a shared_ptr of a function space then provides the function space type
 * \note it checks that FESpace is indeed a functionspace type and return void if it is not the case.
 */
template<typename ElementT>
using functionspace_element_type = typename mpl::if_<is_functionspace_element<decay_type<ElementT> >,
                                                     mpl::identity<decay_type<ElementT>>,
                                                     mpl::identity<void> >::type::type;

/**
 * @brief Specialization for FunctionSpaces
 * 
 * @tparam SpaceType 
 */
template <typename SpaceType>
struct element_type_helper<SpaceType, std::enable_if_t<is_functionspace_v<SpaceType>>> {
    using type = typename decay_type<SpaceType>::element_type;
    using ptrtype = typename decay_type<SpaceType>::element_ptrtype;
};


// forward declaration of Mesh class
template <typename GeoShape, typename T, int Tag, typename IndexT>
class Mesh;

/**
 * @brief Specialization for Meshes
 * 
 * @tparam MeshType 
 */
template <typename MeshType>
struct element_type_helper<MeshType, std::enable_if_t<is_mesh_v<MeshType>>> 
{
    using type = typename decay_type<MeshType>::element_type;
    using ptrtype = typename decay_type<MeshType>::element_ptrtype;
};

/**
 * @brief get the mesh type of topological dimension \p d+1 if the real dimension of the mesh is \p d+1
 * 
 * @tparam MeshTyp type of the mesh
 * @tparam TheTag tag of the mesh
 */
template <typename MeshType, int TheTag = 0>
struct elements_mesh
{
    using type = MeshType;
    using ptrtype = std::shared_ptr<type>;
    using const_ptrtype = std::shared_ptr<const type>;
};
/**
 * @brief get the type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using elements_mesh_t = typename elements_mesh<MeshType,TheTag>::type;

/**
 * @brief get the shared_ptr type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using elements_mesh_ptr_t = typename elements_mesh<MeshType,TheTag>::ptrtype;

/**
 * @brief get the shared_ptr type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using elements_mesh_cptr_t = typename elements_mesh<MeshType,TheTag>::const_ptrtype;


/**
 * @brief get the mesh type of topological dimension \p d+1 if the real dimension of the mesh is \p d+1
 * 
 * @tparam MeshTyp type of the mesh
 * @tparam TheTag tag of the mesh
 */
template <typename MeshType, int TheTag = 0>
struct faces_mesh
{
    using index_t = typename MeshType::index_type;
    using value_type = typename MeshType::value_type;
    static constexpr int nOrder = MeshType::nOrder;
    using type = mp11::mp_if_c<
        is_simplex_mesh_v<MeshType>,
        Mesh<Simplex<2, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>,
        Mesh<Hypercube<2, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>
    >;
    using ptrtype = std::shared_ptr<type>;
    using const_ptrtype = std::shared_ptr<const type>;
};
/**
 * @brief get the type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using faces_mesh_t = typename faces_mesh<MeshType,TheTag>::type;

/**
 * @brief get the shared_ptr type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using faces_mesh_ptr_t = typename faces_mesh<MeshType,TheTag>::ptrtype;

/**
 * @brief get the shared_ptr type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using faces_mesh_cptr_t = typename faces_mesh<MeshType,TheTag>::const_ptrtype;


/**
 * @brief get the mesh type of topological dimension \p d+1 if the real dimension of the mesh is \p d+1
 * 
 * @tparam MeshTyp type of the mesh
 * @tparam TheTag tag of the mesh
 */
template <typename MeshType, int TheTag = 0>
struct facets_mesh
{
    using index_t = typename MeshType::index_type;
    using value_type = typename MeshType::value_type;
    static constexpr int nOrder = MeshType::nOrder;
    static constexpr int d = MeshType::nDim-1;
    using type = mp11::mp_if_c<
        is_simplex_mesh_v<MeshType>,
        Mesh<Simplex<d, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>,
        Mesh<Hypercube<d, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>
    >;
    using ptrtype = std::shared_ptr<type>;
    using const_ptrtype = std::shared_ptr<const type>;
};
/**
 * @brief get the type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using facets_mesh_t = typename facets_mesh<MeshType,TheTag>::type;

/**
 * @brief get the shared_ptr type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using facets_mesh_ptr_t = typename facets_mesh<MeshType,TheTag>::ptrtype;

/**
 * @brief get the shared_ptr type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using facets_mesh_cptr_t = typename facets_mesh<MeshType,TheTag>::const_ptrtype;

/**
 * @brief get the mesh type of topological dimension \p d+1 if the real dimension of the mesh is \p d+1
 * 
 * @tparam MeshTyp type of the mesh
 * @tparam TheTag tag of the mesh
 */
template <typename MeshType, int TheTag = 0>
struct edges_mesh
{
    using index_t = typename MeshType::index_type;
    using value_type = typename MeshType::value_type;
    static constexpr int nOrder = MeshType::nOrder;
    using type = mp11::mp_if_c<
        is_simplex_mesh_v<MeshType>,
        Mesh<Simplex<1, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>,
        Mesh<Hypercube<1, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>
    >;
    using ptrtype = std::shared_ptr<type>;
    using const_ptrtype = std::shared_ptr<const type>;
};
/**
 * @brief get the type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using edges_mesh_t = typename edges_mesh<MeshType,TheTag>::type;

/**
 * @brief get the shared_ptr type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using edges_mesh_ptr_t = typename edges_mesh<MeshType,TheTag>::ptrtype;

/**
 * @brief get the shared_ptr type of the mesh of topological dimension 2
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using edges_mesh_cptr_t = typename edges_mesh<MeshType,TheTag>::const_ptrtype;


/**
 * @brief get the mesh type of topological dimension \p d+1 if the real dimension of the mesh is \p d+1
 * 
 * @tparam MeshTyp type of the mesh
 * @tparam TheTag tag of the mesh
 */
template <typename MeshType, int TheTag = 0>
struct parent_mesh
{
    using index_t = typename MeshType::index_type;
    using value_type = typename MeshType::value_type;
    static constexpr int nOrder = MeshType::nOrder;
    static constexpr int d = (MeshType::nDim == MeshType::nRealDim - 1) ? MeshType::nDim + 1 : MeshType::nDim;
    using type = mp11::mp_if_c<
        is_simplex_mesh_v<MeshType>,
        Mesh<Simplex<d, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>,
        Mesh<Hypercube<d, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>
    >;
    using ptrtype = std::shared_ptr<type>;
    using const_ptrtype = std::shared_ptr<const type>;
};
/**
 * @brief get the type of the mesh of topological dimension \p d+1 if the real dimension of the mesh is \p d+1
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using parent_mesh_t = typename parent_mesh<MeshType,TheTag>::type;

/**
 * @brief get the shared_ptr type of the mesh of topological dimension \p d+1 if the real dimension of the mesh is \p d+1
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using parent_mesh_ptr_t = typename parent_mesh<MeshType,TheTag>::ptrtype;

/**
 * @brief get the shared_ptr type of the mesh of topological dimension \p d+1 if the real dimension of the mesh is \p d+1
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using parent_mesh_cptr_t = typename parent_mesh<MeshType,TheTag>::const_ptrtype;

/**
 * @brief helper class to get the type of the trace of a mesh
 *
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
struct trace_mesh
{
    using index_t = typename MeshType::index_type;
    using value_type = typename MeshType::value_type;
    static constexpr int nOrder = MeshType::nOrder;
    static constexpr int d = (MeshType::nDim == 0) ? 0 : MeshType::nDim - 1;
    using type = mp11::mp_if_c<
        is_simplex_mesh_v<MeshType>,
        Mesh<Simplex<d, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>,
        Mesh<Hypercube<d, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>
    >;
    using ptrtype = std::shared_ptr<type>;
    using const_ptrtype = std::shared_ptr<const type>;
};
/**
 * @brief get type of the trace of a mesh
 *
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType,int TheTag = 0>
using trace_mesh_t = typename trace_mesh<MeshType,TheTag>::type;

/**
 * @brief get the shared_ptr type of the trace of a mesh
 *
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType,int TheTag = 0>
using trace_mesh_ptr_t = typename trace_mesh<MeshType,TheTag>::ptrtype;

/**
 * @brief get the shared_ptr type of the trace of a const mesh
 *
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType,int TheTag = 0>
using trace_mesh_cptr_t = typename trace_mesh<MeshType,TheTag>::const_ptrtype;

/**
 * @brief helper class to get the type of the trace of the trace of a mesh
 *
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template <typename MeshType,int TheTag = 0>
struct trace_trace_mesh
{
    using index_t = typename MeshType::index_type;
    using value_type = typename MeshType::value_type;
    static constexpr int nOrder = MeshType::nOrder;
    static inline constexpr uint16_type nDim = (MeshType::nDim == 1) ? MeshType::nDim - 1 : MeshType::nDim - 2;
    using type = mp11::mp_if_c<
        is_simplex_mesh_v<MeshType>,
        Mesh<Simplex<nDim, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>,
        Mesh<Hypercube<nDim, nOrder, MeshType::nRealDim>, value_type, TheTag, index_t>
    >;
    using ptrtype = std::shared_ptr<type>;
    using const_ptrtype = std::shared_ptr<const type>;
};

/**
 * @brief get type of the trace of the trace of a mesh
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using trace_trace_mesh_t = typename trace_trace_mesh<MeshType,TheTag>::type;

/**
 * @brief get the shared_ptr type of the trace of the trace of a mesh
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using trace_trace_mesh_ptr_t = typename trace_trace_mesh<MeshType, TheTag>::ptrtype;

/**
 * @brief get the shared_ptr type of the trace of the trace of a const mesh
 * 
 * @tparam MeshType type of the mesh
 * @tparam TheTag tag of the mesh
 */
template<typename MeshType, int TheTag = 0>
using trace_trace_mesh_cptr_t = typename trace_trace_mesh<MeshType, TheTag>::const_ptrtype;
/**
 * @} // end Traits group
 */
}
#endif
