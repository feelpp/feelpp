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
using is_mesh = typename std::is_base_of<MeshBase<>,MeshType>::type;

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
using is_functionspace = typename std::is_base_of<FunctionSpaceBase,FuncSpaceType>::type;

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
 * get the element type of a functionspace
 */
template<typename SpaceT, typename = std::enable_if_t<is_functionspace_v<SpaceT>>>
using element_t = typename SpaceT::element_type;
/**
 * get the element shared ptr type of a functionspace
 */
template<typename SpaceT, typename = std::enable_if_t<is_functionspace_v<SpaceT>>>
using element_ptr_t = typename SpaceT::element_ptrtype;

/**
 * @} // end Traits group
 */
}
#endif
