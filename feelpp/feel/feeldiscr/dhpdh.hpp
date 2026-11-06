/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013-2016 Feel++ Consortium

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
#if !defined(FEELPP_DHPDH_HPP)
#define FEELPP_DHPDH_HPP 1

#include <feel/feeldiscr/dh.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/productfunctionspaces.hpp>

namespace Feel {

#if FEELPP_ENABLE_LEGACY_COMPOSITE_FUNCTIONSPACE
template<typename MeshType, int Order = Dynamic, typename T = double>
using DhPdh_type = FunctionSpace<MeshType,
                                 bases<RaviartThomas<Order>,Lagrange<Order,Scalar,Discontinuous>>,
                                 T>;

template<typename MeshType, int Order = Dynamic, typename T = double>
using DhPdh_ptrtype = std::shared_ptr<DhPdh_type<MeshType,Order,T>>;
#endif // FEELPP_ENABLE_LEGACY_COMPOSITE_FUNCTIONSPACE

template<typename MeshType, int Order = Dynamic, typename T = double>
using DhPdh_product_flux_space_type = RTh_type<MeshType,Order,T>;

template<typename MeshType, int Order = Dynamic, typename T = double>
using DhPdh_product_pressure_space_type = FunctionSpace<MeshType,
                                                        bases<Lagrange<Order,Scalar,Discontinuous>>,
                                                        T>;

template<typename MeshType, int Order = Dynamic, typename T = double>
using DhPdh_product_flux_space_ptrtype = std::shared_ptr<DhPdh_product_flux_space_type<MeshType,Order,T>>;

template<typename MeshType, int Order = Dynamic, typename T = double>
using DhPdh_product_pressure_space_ptrtype = std::shared_ptr<DhPdh_product_pressure_space_type<MeshType,Order,T>>;

template<typename MeshType, int Order = Dynamic, typename T = double>
using DhPdh_product_type = product_function_spaces_t<DhPdh_product_flux_space_ptrtype<MeshType,Order,T>,
                                                     DhPdh_product_pressure_space_ptrtype<MeshType,Order,T>>;

template<typename MeshType, int Order = Dynamic, typename T = double>
using DhPdh_product_ptrtype = std::shared_ptr<DhPdh_product_type<MeshType,Order,T>>;

template<typename MeshType, int Order = Dynamic, typename T = double>
using DhPdhProduct_type = DhPdh_product_type<MeshType,Order,T>;

template<typename MeshType, int Order = Dynamic, typename T = double>
using DhPdhProduct_ptrtype = DhPdh_product_ptrtype<MeshType,Order,T>;

#if FEELPP_ENABLE_LEGACY_COMPOSITE_FUNCTIONSPACE
/**

   \code
   #include <feel/feeldiscr/dhpdh.hpp>
   auto Xh = DhPdh<2>( mesh );
   \endcode
 */
template<int Order = Dynamic,typename MeshType, typename T = double>
inline
DhPdh_ptrtype<MeshType,Order,T>
DhPdh( std::shared_ptr<MeshType> mesh,
       RuntimeOrder order,
       std::vector<DofTableExtendedType> dte = std::vector<DofTableExtendedType>( 2, DofTableExtendedType::DEFAULT ) )
{
    CHECK( dte.size() == 2 ) << " vector activation for extended dof table must be equal to 2 but here " << dte.size();
    return DhPdh_type<MeshType,Order,T>::New( _mesh=mesh,
                                              _worldscomm=makeWorldsComm( 2,mesh->worldComm() ),
                                              _runtime_order=order,
                                              _extended_doftable=dte );
}
#endif // FEELPP_ENABLE_LEGACY_COMPOSITE_FUNCTIONSPACE

template<int Order = Dynamic,typename MeshType, typename T = double>
inline
DhPdh_product_ptrtype<MeshType,Order,T>
DhPdhProduct( std::shared_ptr<MeshType> mesh,
              RuntimeOrder order,
              std::vector<DofTableExtendedType> dte = std::vector<DofTableExtendedType>( 2, DofTableExtendedType::DEFAULT ) )
{
    CHECK( dte.size() == 2 ) << " vector activation for extended dof table must be equal to 2 but here " << dte.size();
    using pressure_space_type = DhPdh_product_pressure_space_type<MeshType,Order,T>;
    auto Dh = RTh<Order>( mesh, order, dte[0] );
    auto Pd = pressure_space_type::New( _mesh=mesh,
                                        _worldscomm=makeWorldsComm( 1,mesh->worldComm() ),
                                        _runtime_order=order,
                                        _extended_doftable=dte[1] );
    return productFunctionSpacesPtr( Dh, Pd );
}

#if FEELPP_ENABLE_LEGACY_COMPOSITE_FUNCTIONSPACE
template<int Order,typename MeshType, typename T = double>
requires ( Order >= 0 )
inline
DhPdh_ptrtype<MeshType,Order,T>
DhPdh( std::shared_ptr<MeshType> mesh,
       std::vector<DofTableExtendedType> dte = std::vector<DofTableExtendedType>( 2, DofTableExtendedType::DEFAULT ) )
{
    CHECK( dte.size() == 2 ) << " vector activation for extended dof table must be equal to 2 but here " << dte.size();
    return DhPdh_type<MeshType,Order,T>::New( _mesh=mesh,
                                              _worldscomm=makeWorldsComm( 2,mesh->worldComm() ),
                                              _runtime_order=RuntimeOrder{ static_cast<uint16_type>( Order ) },
                                              _extended_doftable=dte );
}
#endif // FEELPP_ENABLE_LEGACY_COMPOSITE_FUNCTIONSPACE

template<int Order,typename MeshType, typename T = double>
requires ( Order >= 0 )
inline
DhPdh_product_ptrtype<MeshType,Order,T>
DhPdhProduct( std::shared_ptr<MeshType> mesh,
              std::vector<DofTableExtendedType> dte = std::vector<DofTableExtendedType>( 2, DofTableExtendedType::DEFAULT ) )
{
    return DhPdhProduct<Order,MeshType,T>( mesh,
                                           RuntimeOrder{ static_cast<uint16_type>( Order ) },
                                           dte );
}


}
#endif /* FEELPP_DHPDH_HPP */
