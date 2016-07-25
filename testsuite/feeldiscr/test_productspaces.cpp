/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 24 Jul 2016

 Copyright (C) 2016 Feel++ Consortium

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
#define BOOST_TEST_MODULE test_product
#include <testsuite/testsuite.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test space product Options" );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_product" ,
                     "test_product" ,
                     "0.2",
                     "test product of spaces",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );

    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

 FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
 BOOST_AUTO_TEST_SUITE( productspace_suite )


 BOOST_AUTO_TEST_CASE( test1 )
 {
     using namespace Feel;
     using Feel::cout;

     auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
     auto Xh=Pch<1>(mesh);
     auto Wh=Pchv<2>(mesh);
     auto sp = product(Xh,Wh);
     BOOST_CHECK_EQUAL( sp.numberOfSpaces(), 2 );
     cout << tc::red << "number of spaces " << tc::reset << sp.numberOfSpaces() << std::endl;
     BOOST_CHECK_EQUAL( sp.nDof(), Xh->nDof()+Wh->nDof() );
     cout << tc::red << "total number of dof " << tc::reset << sp.nDof() << std::endl;
     BOOST_CHECK_EQUAL( sp.nLocalDof(), Xh->nLocalDof()+Wh->nLocalDof() );
     cout << tc::red << "local number of dof " << tc::reset << sp.nLocalDof() << std::endl;
     BOOST_CHECK_EQUAL( sp[0_c], Xh );
     BOOST_CHECK_EQUAL( sp[1_c], Wh );

     auto cp = hana::cartesian_product(hana::make_tuple(sp,sp));

     BOOST_CHECK_EQUAL( cp[0_c][0_c], Xh );
     BOOST_CHECK_EQUAL( cp[0_c][1_c], Xh );
     BOOST_CHECK_EQUAL( cp[1_c][0_c], Xh );
     BOOST_CHECK_EQUAL( cp[1_c][1_c], Wh );
     BOOST_CHECK_EQUAL( cp[2_c][0_c], Wh );
     BOOST_CHECK_EQUAL( cp[2_c][1_c], Xh );
     BOOST_CHECK_EQUAL( cp[3_c][0_c], Wh );
     BOOST_CHECK_EQUAL( cp[3_c][1_c], Wh );

     auto ps = product(Xh,Wh);
     auto u = Xh->element();
     auto v = Wh->element();
     auto bbf = blockform2( ps );
     bbf( 0_c, 0_c ) = integrate( _range=elements(mesh), _expr=id(u)*idt(u) );
     bbf( 0_c, 1_c ) += integrate( _range=elements(mesh), _expr=id(u)*(trans(idt(v))*one() ) ) ;
     bbf( 1_c, 0_c ) += integrate( _range=elements(mesh), _expr=1./2*(trans(id(v))*one())*idt(u) );
     bbf( 1_c, 1_c ) += integrate( _range=elements(mesh), _expr=trans(id(v))*idt(v) );
     bbf.close();

     auto blf = blockform1( ps );
     blf( 0_c) = integrate( _range=elements(mesh), _expr=id(u) );
     blf( 1_c) += integrate( _range=elements(mesh), _expr=trans(id(v))*one() );
     blf.close();

     auto U=ps.element();
     bbf.solve( _solution=U, _rhs=blf );
     auto ex = exporter(_mesh=mesh);
     ex->add("u",U[0_c]);
     ex->add("v",U[1_c]);
     ex->save();
 }

BOOST_AUTO_TEST_CASE( test2 )
{
    using namespace Feel;
    using Feel::cout;
    using namespace boost::hana::literals;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Xh=Pch<1>(mesh);
    auto Wh=Pchv<2>(mesh);
    auto Zh=Pdh<2>(mesh);
    auto sp = product(Xh,Wh,Zh);
    BOOST_CHECK_EQUAL( sp.numberOfSpaces(), 3 );
    cout << tc::red << "number of spaces " << tc::reset << sp.numberOfSpaces() << std::endl;
    BOOST_CHECK_EQUAL( sp.nDof(), Xh->nDof()+Wh->nDof()+Zh->nDof() );
    cout << tc::red << "total number of dof " << tc::reset << sp.nDof() << std::endl;
    BOOST_CHECK_EQUAL( sp.nLocalDof(), Xh->nLocalDof()+Wh->nLocalDof()+Zh->nLocalDof());
    cout << tc::red << "local number of dof " << tc::reset << sp.nLocalDof() << std::endl;
    BOOST_CHECK_EQUAL( sp[0_c], Xh );
    BOOST_CHECK_EQUAL( sp[1_c], Wh );
    BOOST_CHECK_EQUAL( sp[2_c], Zh );
}

BOOST_AUTO_TEST_SUITE_END()
