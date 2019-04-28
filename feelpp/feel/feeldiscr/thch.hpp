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
   \file thch.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_THCH_HPP)
#define FEELPP_THCH_HPP 1

#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {

namespace meta {
template<int Order,typename MeshType>
struct THch
{
    typedef FunctionSpace<MeshType,
                          bases<Lagrange<Order+1,Vectorial>,Lagrange<Order,Scalar>>,
                          double,
                          Periodicity <NoPeriodicity,NoPeriodicity>,
                          mortars<NoMortar,NoMortar> > type;
    typedef std::shared_ptr<type> ptrtype;
};

} //meta

/**
 * Define the type for Taylor-Hood space 
 * \code
 * THch_type<1,Mesh<Simplex<2>>> // generates \f$P2P1\f$ over a mesh of triangles
 * \endcode
 */
template<int Order,typename MeshType>
using THch_type = FunctionSpace<MeshType,
                                bases<Lagrange<Order+1,Vectorial>,Lagrange<Order,Scalar>>,
                                double,
                                Periodicity <NoPeriodicity,NoPeriodicity>,
                                mortars<NoMortar,NoMortar> >;
/**
 * Define the shared_ptr type for Taylor-Hood space 
 * \code
 * THch_ptrtype<1,Mesh<Simplex<2>>> // defines the shared_ptr type of \f$P2P1\f$ over a mesh of triangles
 * \endcode
 */
template<int Order,typename MeshType>
using THch_ptrtype = std::shared_ptr<FunctionSpace<MeshType,
                                                     bases<Lagrange<Order+1,Vectorial>,Lagrange<Order,Scalar>>,
                                                     double,
                                                     Periodicity <NoPeriodicity,NoPeriodicity>,
                                                     mortars<NoMortar,NoMortar> >>;
template<int Order,typename MeshType>
using THch_velocity_space_t = typename THch_type<Order,MeshType>::template sub_functionspace_type<0>;
template<int Order,typename MeshType>
using THch_pressure_space_t = typename THch_type<Order,MeshType>::template sub_functionspace_type<1>;

template<int Order,typename MeshType>
using THch_element_t = typename THch_type<Order,MeshType>::element_type;
template<int Order,typename MeshType>
using THch_velocity_t = typename THch_type<Order,MeshType>::element_type::template sub_element_type<0>;
template<int Order,typename MeshType>
using THch_pressure_t = typename THch_type<Order,MeshType>::element_type::template sub_element_type<1>;

template<int Order,typename MeshType>
using THch_velocity_space_ptr_t = typename THch_type<Order,MeshType>::template sub_functionspace_ptrtype<0>;
template<int Order,typename MeshType>
using THch_pressure_space_ptr_t = typename THch_type<Order,MeshType>::template sub_functionspace_ptrtype<1>;
    
/**
   Given a \p mesh and polynomial order \f$k\f$(template argument), build a
   product function space of \f$[P_{k+1}]^d \times P_{k}]\f$ where $d$ is the
   dimension of the associated mesh. This kind of function space can be used for
   Stokes problems where the first space is associated to the velocity and the
   second one to the pressure.

   \code
   auto Xh = THch<2>( mesh );
   \endcode
 */
template<int Order,typename MeshType>
inline
THch_ptrtype<Order,MeshType>
THch( std::shared_ptr<MeshType> mesh,
      std::vector<bool> buildExtendedDofTable = std::vector<bool>( 2,false ) )
{
    CHECK( buildExtendedDofTable.size() == 2 ) << " vector activation for extended dof table must be equal to 2 but here " << buildExtendedDofTable.size() << "\n";
    return THch_type<Order,MeshType>::New( _mesh=mesh,
                                           _worldscomm=makeWorldsComm( 2,mesh->worldComm() ),
                                           _extended_doftable=buildExtendedDofTable );
}


}
#endif /* FEELPP_THCH_HPP */
