/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2014-2016 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file dh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-01-30
 */
#ifndef FEELPP_DH_H
#define FEELPP_DH_H 1

#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/order.hpp>

namespace Feel {

template<int Order,typename MeshType,typename T = double>
using dh_type = FunctionSpace<MeshType,bases<RaviartThomas<Order>>,T>;

template<int Order,typename MeshType>
using dh_ptrtype = std::shared_ptr<dh_type<Order,MeshType>>;


template<typename MeshType, int Order = Dynamic,typename T = double>
using Dh_type = FunctionSpace<MeshType,bases<RaviartThomas<Order>>,T>;

template<typename MeshType, int Order = Dynamic>
using Dh_ptrtype = std::shared_ptr<dh_type<Order,MeshType>>;

/**
 * \fn Dh<k,MeshType>
 *
 * build a function space of continuous function which are piecewise polynomial
 * of degree (total or in each variable) less than k.
 */
template<int Order = Dynamic,typename MeshType>
inline
dh_ptrtype<Order,MeshType>
Dh( std::shared_ptr<MeshType> const& mesh,
    RuntimeOrder order,
    DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return dh_type<Order,MeshType>::New( _mesh=mesh,
                                         _worldscomm=makeWorldsComm( 1, mesh->worldComm() ),
                                         _extended_doftable=dte,
                                         _runtime_order=order );
}

template<int Order,typename MeshType>
    requires ( Order >= 0 )
inline
dh_ptrtype<Order,MeshType>
Dh( std::shared_ptr<MeshType> const& mesh,
    DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Dh<Order>( mesh, RuntimeOrder{ static_cast<uint16_type>( Order ) }, dte );
}

template<typename MeshType, int Order, typename T = double>
using RTh_type = Dh_type<MeshType, Order, T>;

template<typename MeshType, int Order, typename T = double>
using RTh_ptrtype = std::shared_ptr<RTh_type<MeshType, Order, T>>;

template<int Order = Dynamic, typename MeshType>
inline
RTh_ptrtype<MeshType, Order>
RTh( std::shared_ptr<MeshType> const& mesh,
     RuntimeOrder order,
     DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Dh<Order>( mesh, order, dte );
}

template<int Order, typename MeshType>
    requires ( Order >= 0 )
inline
RTh_ptrtype<MeshType, Order>
RTh( std::shared_ptr<MeshType> const& mesh,
     DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Dh<Order>( mesh, dte );
}

} // Feel

#endif /* FEELPP_DH_H */
