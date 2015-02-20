/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-06-20

  Copyright (C) 2014 Feel++ Consortium

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
   \file test_prodspace.cpp
   \author Cecile Daversin <daversin@math.unistra.fr>
   \author Vincent Chabannes
   \date 2014-06-20
 */
#if 1
#define BOOST_TEST_MODULE test_project_prodspace
#include <testsuite/testsuite.hpp>
#endif

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

//int main( int argc, char **argv)
void run()
{
    /*mesh*/
    typedef Simplex<2,1> entity_type;
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef Lagrange<1, Scalar> basis_type;
    typedef bases<basis_type, basis_type> prod_basis_type;

    /*spaces*/
    typedef FunctionSpace<mesh_type, bases<basis_type>, double> space_type;
    typedef FunctionSpace<mesh_type, prod_basis_type, double > prod_space_type;

    auto mesh=unitCircle();
    auto pXh = prod_space_type::New( mesh );
    auto Xh1 = pXh->functionSpace<0>();
    auto Xh2 = pXh->functionSpace<1>();
    auto Xh = space_type::New( mesh );
    size_type nLocDofXh = Xh->nLocalDofWithGhost();
    size_type nDofXh = Xh->nDof();
    size_type nLocDofpXh1 = Xh1->nLocalDofWithGhost();
    size_type nLocDofpXh2 = Xh2->nLocalDofWithGhost();
    size_type nDofpXh1 = Xh1->nDof();
    size_type nDofpXh2 = Xh2->nDof();

    auto U1 = pXh->element("M_U1");
    auto U2 = Xh->element("M_U2");
    if ( Environment::isMasterRank() )
    {
        std::cout << "Xh nDof = " << Xh->nDof() << std::endl;
        std::cout << "pXh nDof = " << pXh->nDof() << std::endl;
    }

    for ( size_type k=0 ; k<nLocDofXh ; ++k )
        U2.set( k, 1.0 );
    for ( size_type k=0 ; k<nLocDofpXh1 ; ++k )
        U1.set( k, 2.0 );
    for ( size_type k=0 ; k<nLocDofpXh2 ; ++k )
        U1.set( nLocDofpXh1+k, 3.0 );


    auto myproj1 = vf::project(_space=Xh, _expr= idv(U2));
    size_type sumProj1 = myproj1.sum();
    BOOST_CHECK( sumProj1 == 1*nDofXh );
    auto myproj1b = Xh->element();
    myproj1b.on(_range=elements(mesh),_expr=idv(U2) );
    size_type sumProj1b = myproj1b.sum();
    BOOST_CHECK( sumProj1b == 1*nDofXh );
    if ( Environment::isMasterRank() )
        std::cout << "classical proj ok" << std::endl;

    auto myproj31 = vf::project(_space=Xh1, _expr= idv(U1.element<0>()));
    size_type sumProj31 = myproj31.sum();
    BOOST_CHECK( sumProj31 == 2*nDofpXh1 );
    auto myproj32 = vf::project(_space=Xh2, _expr= idv(U1.element<1>()));
    size_type sumProj32 = myproj32.sum();
    BOOST_CHECK( sumProj32 == 3*nDofpXh2 );

    auto myproj3b = pXh->element();
    myproj3b.element<0>().on(_range=elements(mesh),_expr=idv(U1.element<0>()) );
    myproj3b.element<1>().on(_range=elements(mesh),_expr=idv(U1.element<1>()) );
    size_type sumProj3b1 = myproj3b.element<0>().sum();
    BOOST_CHECK( sumProj3b1 == 2*nDofpXh1 );
    size_type sumProj3b2 = myproj3b.element<1>().sum();
    BOOST_CHECK( sumProj3b2 == 3*nDofpXh2 );
    size_type sumProj3b = myproj3b.sum();
    BOOST_CHECK( sumProj3b == ( 2*nDofpXh1+3*nDofpXh2) );
    if ( Environment::isMasterRank() )
        std::cout << "view proj (directly in expr) ok" << std::endl;

    auto myexpr_view = idv(U1.element<0>());
    auto myproj4 = vf::project(_space=Xh, _expr=myexpr_view );
    size_type sumProj4 = myproj4.sum();
    BOOST_CHECK( sumProj4 == 2*nDofpXh1 );
    if ( Environment::isMasterRank() )
        std::cout << "view proj (with copy) ok" << std::endl;

}


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( project_prodspace )

BOOST_AUTO_TEST_CASE( project_prodspace )
{
    run();
}

BOOST_AUTO_TEST_SUITE_END()
