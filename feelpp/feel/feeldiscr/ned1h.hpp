/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2014-2016 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file ned1h.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-01-30
 */
#ifndef FEELPP_NED1H_H
#define FEELPP_NED1H_H 1

#include <feel/feelpoly/nedelec.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/fecapabilities.hpp>
#include <feel/feelpoly/order.hpp>

namespace Feel {

namespace meta {
template<typename MeshType,int Order = Dynamic,typename T = double>
struct Ned1h
{
    typedef FunctionSpace<MeshType,
                          bases<Nedelec<Order,NedelecKind::NED1>>,
                          T,
                          mortars<NoMortar>> type;
    typedef std::shared_ptr<type> ptrtype;
};

} // meta

template<typename MeshType,int Order = Dynamic,typename T = double>
using Ned1h_type = typename meta::Ned1h<MeshType,Order,T>::type;

template<typename MeshType,int Order = Dynamic, typename T = double>
using Ned1h_ptrtype = typename meta::Ned1h<MeshType,Order,T>::ptrtype;

/**
 * \fn Ned1h<k,MeshType>
 *
 */
template<int Order = Dynamic,typename MeshType,typename T=double>
inline
typename meta::Ned1h<MeshType,Order,T>::ptrtype
Ned1h( std::shared_ptr<MeshType> const& mesh,
       RuntimeOrder order,
       DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    fe::requireCapability<fe::Family::NedelecFirstKind>( *mesh, order );
    typedef typename meta::Ned1h<MeshType,Order,T>::type space_type;
    return space_type::New( _mesh=mesh,
                            _worldscomm=makeWorldsComm( 1, mesh->worldCommPtr() ),
                            _extended_doftable=dte,
                            _runtime_order=order );
}

template<int Order,typename MeshType,typename T=double>
    requires ( Order >= 0 )
inline
typename meta::Ned1h<MeshType,Order,T>::ptrtype
Ned1h( std::shared_ptr<MeshType> const& mesh,
       DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Ned1h<Order, MeshType, T>( mesh,
                                      RuntimeOrder{ static_cast<uint16_type>( Order ) },
                                      dte );
}

/**
 * \fn Ned1h<k,MeshType>
 *
 */
template<int Order,typename MeshType, typename RangeType, typename T = double>
    requires ( Order >= 0 )
inline
typename meta::Ned1h<MeshType,Order,T>::ptrtype
Ned1h( std::shared_ptr<MeshType> const& mesh, RangeType && rangeElt, DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    return Ned1h<Order, MeshType, RangeType, T>( mesh,
                                                 std::forward<RangeType>( rangeElt ),
                                                 RuntimeOrder{ static_cast<uint16_type>( Order ) },
                                                 dte );
}

/**
 * \fn Ned1h<k,MeshType>
 *
 */
template<int Order = Dynamic,typename MeshType, typename RangeType, typename T = double>
inline
typename meta::Ned1h<MeshType,Order,T>::ptrtype
Ned1h( std::shared_ptr<MeshType> const& mesh,
       RangeType && rangeElt,
       RuntimeOrder order,
       DofTableExtendedType dte = DofTableExtendedType::DEFAULT )
{
    typedef typename meta::Ned1h<MeshType,Order,T>::type space_type;
    return space_type::New( _mesh=mesh,
                            _range=std::forward<RangeType>(rangeElt),
                            _extended_doftable=dte,
                            _runtime_order=order );
}

} // Feel

#endif /* FEELPP_NED1H_H */
