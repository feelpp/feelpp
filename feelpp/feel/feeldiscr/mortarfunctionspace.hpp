/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

/*
  This file is part of the Feel library

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
   \file mortarfunctionspace.hpp
   \brief Explicit mortar-space facade for the compatibility migration.
 */
#if !defined(FEELPP_MORTARFUNCTIONSPACE_HPP)
#define FEELPP_MORTARFUNCTIONSPACE_HPP 1

#include <memory>
#include <type_traits>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/order.hpp>

namespace Feel
{

/**
 * \brief Explicit mortar function-space abstraction.
 *
 * This class keeps the existing FunctionSpace/DofTable mortar implementation
 * underneath while giving new code an explicit mortar type and construction
 * point.  The DofTable template still receives the mortar policy during this
 * compatibility phase.
 */
template<typename MeshTypes,
         typename BasisTypes = bases<Lagrange<1, Scalar>>,
         typename T = double,
         typename MortarPolicy = Mortar>
class MortarFunctionSpace
    : public FunctionSpace<MeshTypes, BasisTypes, T, mortars<MortarPolicy>>
{
    static_assert( MortarPolicy::is_mortar, "MortarFunctionSpace requires a mortar policy" );

public:
    using super_type = FunctionSpace<MeshTypes, BasisTypes, T, mortars<MortarPolicy>>;
    using self_type = MortarFunctionSpace<MeshTypes, BasisTypes, T, MortarPolicy>;
    using functionspace_type = self_type;
    using base_functionspace_type = super_type;
    using pointer_type = std::shared_ptr<self_type>;
    using mortar_policy_type = MortarPolicy;
    using typename super_type::mesh_ptrtype;
    using typename super_type::mesh_support_vector_type;

    static constexpr bool is_explicit_mortar_space = true;

    using super_type::super_type;

    /**
     * helper static function to create a std::shared_ptr<> out of
     * the explicit mortar function space.
     */
    template <typename ... Ts>
        requires ( sizeof...( Ts ) != 0 ) && ( NA::is_named_argument_v<Ts> && ... )
    static pointer_type
    New( Ts&&... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>( v )... );
        auto&& mesh = args.get( _mesh );
        worldscomm_ptr_t worldscomm =
            args.get_else_invocable( _worldscomm,
                                     [&mesh]() { return Feel::detail::createWorldsComm<self_type>( mesh ).worldsComm(); } );
        size_type components = args.get_else( _components, MESH_RENUMBER | MESH_CHECK );
        auto&& extended_doftable =
            args.get_else( _extended_doftable,
                           std::vector<DofTableExtendedType>( self_type::nSpaces, DofTableExtendedType::DEFAULT ) );
        auto&& range = args.get_else( _range, mesh_support_vector_type() );
        RuntimeOrder runtime_order = args.get_else( _runtime_order, RuntimeOrder{ 0 } );

        auto cms = Feel::detail::createMeshSupport<self_type>( mesh, range );
        std::vector<DofTableExtendedType> edt = Feel::detail::createInfoExtendedDofTable<self_type>( extended_doftable );
        return NewImpl( mesh, cms.M_meshSupportVector, worldscomm, components, edt, runtime_order );
    }

    static pointer_type
    New( mesh_ptrtype const& m )
    {
        return New( _mesh = m );
    }

    static pointer_type
    NewImpl( mesh_ptrtype const& m,
             mesh_support_vector_type const& meshSupport,
             worldscomm_ptr_t const& worldscomm = Environment::worldsComm( self_type::nSpaces ),
             size_type mesh_components = MESH_RENUMBER | MESH_CHECK,
             std::vector<DofTableExtendedType> extendedDofTable =
                 std::vector<DofTableExtendedType>( self_type::nSpaces, DofTableExtendedType::DEFAULT ),
             RuntimeOrder runtime_order = RuntimeOrder{ 0 } )
    {
        return pointer_type( new self_type( m, meshSupport, mesh_components, worldscomm, extendedDofTable, runtime_order ) );
    }
};

template<typename MeshType,
         int Order = Dynamic,
         template<class, int, class> class Pts = PointSetEquiSpaced,
         typename T = double>
using MortarLagrangeSpace =
    MortarFunctionSpace<MeshType, bases<Lagrange<Order, Scalar, Continuous, Pts>>, T, Mortar>;

template<typename MeshType,
         int Order = Dynamic,
         template<class, int, class> class Pts = PointSetEquiSpaced,
         typename T = double>
using MortarLagrangeSpace_ptrtype = std::shared_ptr<MortarLagrangeSpace<MeshType, Order, Pts, T>>;

template<int Order = Dynamic,
         template<class, int, class> class Pts = PointSetEquiSpaced,
         typename MeshType, typename T = double>
inline
MortarLagrangeSpace_ptrtype<MeshType, Order, Pts, T>
mortarFunctionSpace( std::shared_ptr<MeshType> const& mesh,
                     RuntimeOrder order,
                     DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return MortarLagrangeSpace<MeshType, Order, Pts, T>::New( _mesh = mesh,
                                                              _worldscomm = makeWorldsComm( 1, mesh->worldCommPtr() ),
                                                              _runtime_order = order,
                                                              _extended_doftable = dte );
}

template<int Order,
         template<class, int, class> class Pts = PointSetEquiSpaced,
         typename MeshType, typename T = double>
    requires ( Order >= 0 )
inline
MortarLagrangeSpace_ptrtype<MeshType, Order, Pts, T>
mortarFunctionSpace( std::shared_ptr<MeshType> const& mesh,
                     DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return mortarFunctionSpace<Order, Pts, MeshType, T>( mesh,
                                                         RuntimeOrder{ static_cast<uint16_type>( Order ) },
                                                         dte );
}

} // namespace Feel

#endif /* FEELPP_MORTARFUNCTIONSPACE_HPP */
