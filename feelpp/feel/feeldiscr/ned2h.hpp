/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file ned2h.hpp
 * @brief FunctionSpace aliases and factories for second-kind Nedelec elements.
 */
#ifndef FEELPP_FEELDISCR_NED2H_HPP
#define FEELPP_FEELDISCR_NED2H_HPP 1

#include <feel/feelpoly/nedelec2.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/order.hpp>

namespace Feel
{
/** @brief Second-kind Nedelec FunctionSpace type. */
template<typename MeshType, int Order = Dynamic, typename T = double>
using ned2h_type = FunctionSpace<MeshType,
                                 bases<NedelecSecondKind<Order>>,
                                 T,
                                 mortars<NoMortar>>;

/** @brief Shared pointer to a second-kind Nedelec FunctionSpace. */
template<typename MeshType, int Order = Dynamic, typename T = double>
using ned2h_ptrtype = std::shared_ptr<ned2h_type<MeshType, Order, T>>;

/**
 * @brief Build a runtime-order second-kind Nedelec FunctionSpace.
 */
template<int Order = Dynamic, typename MeshType, typename T = double>
[[nodiscard]] inline ned2h_ptrtype<MeshType, Order, T>
Ned2h( std::shared_ptr<MeshType> const& mesh,
       RuntimeOrder order,
       DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    using space_type = ned2h_type<MeshType, Order, T>;
    return space_type::New( _mesh = mesh,
                            _worldscomm = makeWorldsComm( 1, mesh->worldCommPtr() ),
                            _extended_doftable = dte,
                            _runtime_order = order );
}

/**
 * @brief Build a compile-time-order second-kind Nedelec FunctionSpace.
 */
template<int Order, typename MeshType, typename T = double>
    requires ( Order >= 0 )
[[nodiscard]] inline ned2h_ptrtype<MeshType, Order, T>
Ned2h( std::shared_ptr<MeshType> const& mesh,
       DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Ned2h<Order, MeshType, T>(
        mesh, RuntimeOrder{ static_cast<uint16_type>( Order ) }, dte );
}
} // namespace Feel

#endif // FEELPP_FEELDISCR_NED2H_HPP
