/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 26 mai 2015
 
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
#ifndef FEELPP_NCH_HPP
#define FEELPP_NCH_HPP 1

#include <feel/feelpoly/crouzeixraviart.hpp>

#include <feel/feeldiscr/functionspace.hpp>

namespace Feel
{

namespace meta
{

template<typename MeshType,
         int Order,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced,
         typename T=double,
         int Tag = 0>
struct NChv
{
    typedef FunctionSpace<MeshType,
                          bases<CrouzeixRaviart<Order,Vectorial,Pts,Tag>>,
                          T,
                          Periodicity <NoPeriodicity>,
                          mortars<NoMortar>> type;
    typedef std::shared_ptr<type> ptrtype;
};

} // meta

template<typename MeshType,
         int Order,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced,
         typename T=double,
         int Tag = 0>
using NChv_type = typename meta::NChv<MeshType,Order,Pts,T,Tag>::type;
template<typename MeshType,
         int Order,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced,
         typename T=double,
         int Tag = 0>
using NChv_ptrtype = typename meta::NChv<MeshType,Order,Pts,T,Tag>::ptrtype;

/**
   Given a \p mesh, build a function space of vectorial continuous function
   which are piecewise polynomial of degree (total or in each variable) less
   than k using Lagrange basis functions
 */
template<int Order,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced,
         typename MeshType,
         typename T=double,
         int Tag = 0>
inline
NChv_ptrtype<MeshType,Order,Pts,T,Tag>
NChv( std::shared_ptr<MeshType> const& mesh, bool buildExtendedDofTable=false  )
{
    return NChv_type<MeshType,Order,Pts,T,Tag>::New( _mesh=mesh,
                                                   _worldscomm=makeWorldsComm( 1,mesh->worldComm() ),
                                                   _extended_doftable=std::vector<bool>( 1,buildExtendedDofTable ) );
}

}

#endif
