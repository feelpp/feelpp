//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 03 May 2017
//! @copyright 2017 Feel++ Consortium
//!

#if !defined(FEELPP_P2CH_HPP)
#define FEELPP_P2CH_HPP 1

#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {

namespace meta {
template<int Order1, int Order2,typename MeshType>
struct P2ch
{
    typedef FunctionSpace<MeshType,
                          bases<Lagrange<Order1,Vectorial>,Lagrange<Order2,Scalar>>,
                          double,
                          Periodicity <NoPeriodicity,NoPeriodicity>,
                          mortars<NoMortar,NoMortar> > type;
    typedef std::shared_ptr<type> ptrtype;
};

} //meta

/**
 * Define the type for P<Order1>P<Order2>
 * \code
 * P2ch_type<1,1,Mesh<Simplex<2>>> // generates \f$P1P1\f$ over a mesh of triangles
 * \endcode
 */
template<typename Base1,typename Base2,typename MeshType>
using P2ch_type = FunctionSpace<MeshType,
                                bases<Base1, Base2>,
                                double,
                                Periodicity <NoPeriodicity,NoPeriodicity>,
                                mortars<NoMortar,NoMortar> >;
/**
 * Define the shared_ptr type for Taylor-Hood space 
 * \code
 * P2ch_ptrtype<1,Mesh<Simplex<2>>> // defines the shared_ptr type of \f$P2P1\f$ over a mesh of triangles
 * \endcode
 */
template<typename Base1, typename Base2, typename MeshType>
using P2ch_ptrtype = std::shared_ptr<P2ch_type<Base1,Base2,MeshType>>;
    
/**
   Given a \p mesh and polynomial order \f$k\f$(template argument), build a
   product function space of \f$[P_{k}]^d \times P_{l}]\f$ where $d$ is the
   dimension of the associated mesh. This kind of function space can be used for
   Stokes problems where the first space is associated to the velocity and the
   second one to the pressure.

   \code
   auto Xh = P2ch<2,1>( mesh );
   \endcode
 */
template<typename Base1, typename Base2, typename MeshType>
inline
P2ch_ptrtype<Base1, Base2, MeshType>
P2ch( std::shared_ptr<MeshType> mesh,
      std::vector<bool> buildExtendedDofTable = std::vector<bool>( 2,false ) )
{
    CHECK( buildExtendedDofTable.size() == 2 ) << " vector activation for extended dof table must be equal to 2 but here " << buildExtendedDofTable.size() << "\n";
    return P2ch_type<Base1,Base2,MeshType>::New( _mesh=mesh,
                                                 _worldscomm=makeWorldsComm( 2,mesh->worldCommPtr() ),
                                                 _extended_doftable=buildExtendedDofTable );
}


}
#endif /* FEELPP_P2CH_HPP */
