/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file bdmh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-02-18
 */
#ifndef FEELPP_BDMH_H
#define FEELPP_BDMH_H 1

#include <feel/feelpoly/brezzidouglasmarini.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/order.hpp>

namespace Feel {

template<int Order, typename MeshType, typename T = double>
using bdmh_type = FunctionSpace<MeshType, bases<BrezziDouglasMarini<Order>>, T>;

template<int Order, typename MeshType>
using bdmh_ptrtype = std::shared_ptr<bdmh_type<Order, MeshType>>;

template<typename MeshType, int Order = Dynamic, typename T = double>
using BDMh_type = FunctionSpace<MeshType, bases<BrezziDouglasMarini<Order>>, T>;

template<typename MeshType, int Order = Dynamic>
using BDMh_ptrtype = std::shared_ptr<bdmh_type<Order, MeshType>>;

template<int Order = Dynamic, typename MeshType>
inline
bdmh_ptrtype<Order, MeshType>
BDMh( std::shared_ptr<MeshType> const& mesh,
      RuntimeOrder order,
      DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return bdmh_type<Order, MeshType>::New( _mesh = mesh,
                                            _worldscomm = makeWorldsComm( 1, mesh->worldComm() ),
                                            _extended_doftable = dte,
                                            _runtime_order = order );
}

template<int Order, typename MeshType>
    requires ( Order >= 0 )
inline
bdmh_ptrtype<Order, MeshType>
BDMh( std::shared_ptr<MeshType> const& mesh,
      DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return BDMh<Order>( mesh, RuntimeOrder{ static_cast<uint16_type>( Order ) }, dte );
}

} // namespace Feel

#endif /* FEELPP_BDMH_H */
