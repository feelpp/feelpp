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
/**
   \file pchv.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_PCHV_HPP)
#define FEELPP_PCHV_HPP 1

#include <boost/mp11/utility.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel
{

namespace meta
{

template<typename MeshType,
         int Order,         
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename T = double,
         int Tag = 0>
struct Pchv
{
    using type = FunctionSpace<MeshType, bases<Lagrange<Order, Vectorial, Continuous, Pts, Tag>>,T>;
    typedef std::shared_ptr<type> ptrtype;
};

} // meta

template<typename MeshType,
         int Order,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename T = double,
         int Tag = 0>
using Pchv_type = typename meta::Pchv<MeshType,Order,Pts,T,Tag>::type;
template<typename MeshType,
         int Order,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename T = double,
         int Tag = 0>
using Pchv_ptrtype = typename meta::Pchv<MeshType,Order,Pts,T,Tag>::ptrtype;

template<typename MeshType,int Order,template<class, uint16_type, class> class Pts = PointSetFekete,typename T = double,int Tag = 0>
using Pchv_element_t=typename Pchv_type<MeshType,Order,Pts,T,Tag>::element_type;

template<typename MeshType,int Order,template<class, uint16_type, class> class Pts = PointSetFekete,typename T = double,int Tag = 0>
using Pchv_element_type=Pchv_element_t<MeshType,Order,Pts,T,Tag>;


/**
   Given a \p mesh, build a function space of vectorial continuous function
   which are piecewise polynomial of degree (total or in each variable) less
   than k using Lagrange basis functions
 */
template<int Order,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename MeshType,
         typename T=double,
         int Tag = 0>
inline
Pchv_ptrtype<MeshType,Order,Pts,T,Tag>
Pchv( std::shared_ptr<MeshType> const& mesh, bool buildExtendedDofTable=false  )
{
    return Pchv_type<MeshType,Order,Pts,T,Tag>::New( _mesh=mesh,
                                                   _worldscomm=makeWorldsComm(1,mesh->worldComm() ),
                                                   _extended_doftable=buildExtendedDofTable );
}

/**
 Given a \p mesh, build a function space of vectorial continuous function
 which are piecewise polynomial of degree (total or in each variable) less
 than k using Lagrange basis functions
 */
template<int Order,
         template<class, uint16_type, class> class Pts = PointSetFekete,
         typename MeshType,typename RangeType,
         typename T=double,
         int Tag = 0>
inline
Pchv_ptrtype<MeshType,Order,Pts,T,Tag>
Pchv( std::shared_ptr<MeshType> const& mesh, RangeType && rangeElt, bool buildExtendedDofTable=false  )
{
    return Pchv_type<MeshType,Order,Pts,T,Tag>::New( _mesh=mesh,
                                                   _range=std::forward<RangeType>(rangeElt),
                                                   _worldscomm=makeWorldsComm( 1,mesh->worldComm() ),
                                                   _extended_doftable=buildExtendedDofTable );
}


}

#endif /* FEELPP_PCHV_HPP */
