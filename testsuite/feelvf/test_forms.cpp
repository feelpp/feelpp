/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Jan 2016

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
#define BOOST_TEST_MODULE test_forms
#include <testsuite/testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Forms  Options" );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_forms" ,
                     "test_forms" ,
                     "0.2",
                     "nD(n=2,3) test forms",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );

    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}



FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( forms_suite )

using dim_t = boost::mpl::list<boost::mpl::int_<2>, boost::mpl::int_<3> >;
//using dim_t = boost::mpl::list<boost::mpl::int_<2>>;


BOOST_AUTO_TEST_CASE_TEMPLATE( test_form2_faces, T, dim_t )
{
    BOOST_MESSAGE( "test_form2_faces starts for dim=" << T::value);
    auto meshnd = unitHypercube<T::value>();
    auto mesh = createSubmesh( meshnd, internalfaces(meshnd), EXTRACTION_KEEP_MESH_RELATION, 0 );
    size_type nFaceInParallelMeshnd = nelements(internalfaces(meshnd),true) - nelements(interprocessfaces(meshnd),true)/2;
    size_type nInternalFaceInParallelMeshnd = nelements(internalfaces(meshnd),true) - nelements(interprocessfaces(meshnd),true)/2;
    BOOST_CHECK_EQUAL( nelements(elements(mesh),true), nFaceInParallelMeshnd  );
    auto Vh=Pdhv<1>(meshnd,true);
    auto u = Vh->element();
    auto Wh=Pdh<1>(meshnd,true);
    auto p = Wh->element();
    p.on(_range=elements(meshnd),_expr=cst(1.));
    auto Mh=Pdh<1>(mesh);
    auto l=Mh->element();
    l.on(_range=elements(mesh),_expr=cst(1.));

    // all the 4 integrals below should be the same
    // compute \int 1 over internalfaces
    auto I = integrate( _range=internalfaces(meshnd), _expr=cst(1.)).evaluate()(0,0);
    // compute \int 1 over d-1 mesh
    auto I1 = integrate( _range=elements(mesh), _expr=cst(1.)).evaluate()(0,0);
    // compute \int 1 over internalfaces using left element
    auto I2 = integrate( _range=internalfaces(meshnd), _expr=leftfacev(cst(1.))).evaluate()(0,0);
    // compute \int 1 over internalfaces using right element
    auto I3 = integrate( _range=internalfaces(meshnd), _expr=rightfacev(cst(1.))).evaluate()(0,0);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "I=" << I << " I1=" << I1 << " I2=" << I2 << " I3=" << I3 );
    BOOST_CHECK_CLOSE( I, I1, 1e-10 );
    BOOST_CHECK_CLOSE( I, I2, 1e-10 );
    BOOST_CHECK_CLOSE( I, I3, 1e-10 );

    auto a = form2(_test=Mh, _trial=Wh );
    a = integrate( _range=internalfaces(meshnd), _expr=(id(l)/2)*(leftfacet(idt(p))));//+rightfacet(idt(p))));
    a.close();

    if ( Environment::isMasterRank() )
        std::cout << "a(1,1) = " << a(l,p) << std::endl;

    auto a1 = form2(_test=Mh, _trial=Mh );
    a1 = integrate( _range=internalfaces(meshnd), _expr=id(l)*idt(l)/4);
    a1.close();
    auto a1en = a1(l,l);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a1(1,1)=" << a1en );

    auto a11 = form2(_test=Mh, _trial=Mh );
    a11 = integrate( _range=elements(mesh), _expr=id(l)*idt(l));
    a11.close();
    auto a11en = a11(l,l);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a11(1,1)=" << a11en << " int =" << I1 );

    auto a111 = form2(_test=Mh, _trial=Mh );
    a111 = integrate( _range=internalfaces(meshnd), _expr=id(l)*idt(l));
    a111.close();
    auto a111en = a111(l,l);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a111(1,1)=" << a111en << " int =" << I1 );

    // - Tests A22 = <ph, w> with ph, w \in Wh
    // - >>> FAILS <<<
    //
    auto a2 = form2(_test=Wh, _trial=Wh, _pattern=size_type(Pattern::EXTENDED)  );
    a2 = integrate( _range=internalfaces(meshnd), _expr=leftface(id(p))*leftfacet(idt(p)));
    a2.close();
    auto a2en = a2(p,p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a2(1,1)=" << a2en );

    // - Same as a2, but with rightface instead of leftface.
    // - >>> FAILS <<<
    //
    auto a3 = form2(_test=Wh, _trial=Wh, _pattern=size_type(Pattern::EXTENDED)  );
    a3 = integrate(_range=internalfaces(meshnd), _expr=rightface(id(p))*rightfacet(idt(p)));
    a3.close();
    auto a3en = a3(p,p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a3(1,1)=" << a3en );

    auto a23 = form2(_test=Wh, _trial=Wh, _pattern=size_type(Pattern::EXTENDED)  );
    a23 = integrate(_range=internalfaces(meshnd), _expr=rightface(id(p))*leftfacet(idt(p)));
    a23.close();
    auto a23en = a23(p,p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a23(1,1)=" << a23en );

    // - Tests A23 = <phat, w>, with phat \in Mh, w \in Wh
    // - >>> FAILS BUT IT RETURNS A TODO MESSAGE <<<
    //
    auto a4 = form2(_test=Wh, _trial=Mh, _pattern=size_type(Pattern::EXTENDED)  );
    a4 = integrate( _range=internalfaces(meshnd), _expr=0.5*leftface(id(p))*idt(l) );
    a4.close();
    auto a4en = a4(p, l);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a4(1,1)=" << a4en );

    // - Tests A32 = <ph, \mu>, with ph \in Wh, \mu \in Mh
    auto a5 = form2(_test=Mh, _trial=Wh, _pattern=size_type(Pattern::EXTENDED)  );
    a5 = integrate(_range=internalfaces(meshnd), _expr=0.5*leftfacet(idt(p))*id(l) );
    a5.close();
    auto a5en = a5(l, p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a5(1,1)=" << a5en );

    // - Tests (p, p)_Omega. Should give the measure of the domain
    auto a6 = form2(_test=Wh, _trial=Wh);
    a6 = integrate(_range=elements(meshnd), _expr=idt(p)*id(p));
    a6.close();
    auto a6en = a6(p, p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a6(1,1)=" << a6en );

    // The five following tests should all provide the (n-1)-Lebesgue
    // measure of the boundary of the domain.
    // The goal is to check whether we need to multiply by 0.5 or 0.25 or nothing in
    // boundary integrals. Also
    auto a7 = form2( _test=Wh, _trial=Wh);
    a7 = integrate(_range=boundaryfaces(meshnd), _expr=idt(p)*id(p));
    a7.close();
    auto a7eval = a7(p,p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a7(1,1)=" << a7eval );

    // FAILS:
    //
    // Returned message:
    // unknown location(0): fatal error: in "forms_suite/test_form2_faces<N4mpl_4int_ILi2EEE>": memory access violation at address: 0x00000000: no mapping at fault address
    // /data/atlas_dprada/codes/feelpp/testsuite/feelvf/test_forms.cpp(97): last checkpoint
    //
    // auto a8 = form2( _test=Mh, _trial=Mh);
    // a8 = integrate(_range=boundaryfaces(meshnd), _expr=idt(l)*id(l));
    // a8.close();
    // auto a8eval = a8(l,l);
    // if ( Environment::isMasterRank() )
    //     BOOST_TEST_MESSAGE( "a8(1,1)=" << a8eval );

    // FAILS:
    //
    // Just like a8
    // auto a9 = form2(_test=Mh, _trial=Wh);
    // a9 = integrate(_range=boundaryfaces(meshnd), _expr=idt(p)*id(l));
    // a9.close();
    // auto a9eval = a9(l,p);
    // if ( Environment::isMasterRank() )
    //     BOOST_TEST_MESSAGE( "a9(1,1)=" << a9eval );

    // FAILS:
    //
    // Just like a8
    // auto a10 = form1(_test=Mh);
    // a10 = integrate(_range=boundaryfaces(meshnd), _expr=id(l));
    // a10.close();
    // auto a10eval = a10(l);
    // if ( Environment::isMasterRank() )
    //     BOOST_TEST_MESSAGE( "a10(1)=" << a10eval );

    auto a12 = form1(_test=Wh);
    a12 = integrate(_range=boundaryfaces(meshnd), _expr=id(p));
    a12.close();
    auto a12eval = a12(p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a12(1)=" << a12eval );



    // ------------------------------------------------------------------------------------

    // The following tests help to understand what we have to do when mixing
    // integrals that require the extended pattern with others that do not

    // Works fine
    auto a7b = form2( _test=Wh, _trial=Wh, _pattern=size_type(Pattern::EXTENDED));
    a7b = integrate(_range=boundaryfaces(meshnd), _expr=idt(p)*id(p));
    a7b.close();
    auto a7beval = a7b(p,p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a7b(1,1)=" << a7beval );

    // FAILS
    //
    // Just like a8
    // auto a8b = form2( _test=Mh, _trial=Mh, _pattern=size_type(Pattern::EXTENDED));
    // a8b = integrate(_range=boundaryfaces(meshnd), _expr=idt(l)*id(l));
    // a8b.close();
    // auto a8beval = a8b(l,l);
    // if ( Environment::isMasterRank() )
    //     BOOST_TEST_MESSAGE( "a8b(1,1)=" << a8beval );

    // FAILS
    //
    // Just like a8
    // auto a9b = form2(_test=Mh, _trial=Wh, _pattern=size_type(Pattern::EXTENDED));
    // a9b = integrate(_range=boundaryfaces(meshnd), _expr=idt(p)*id(l));
    // a9b.close();
    // auto a9beval = a9b(l,p);
    // if ( Environment::isMasterRank() )
    //     BOOST_TEST_MESSAGE( "a9b(1,1)=" << a9beval );

    // FAILS
    //
    // Just like a8
    // auto a10b = form1(_test=Mh, _pattern=size_type(Pattern::EXTENDED));
    // a10b = integrate(_range=boundaryfaces(meshnd), _expr=id(l));
    // a10b.close();
    // auto a10beval = a10b(l);
    // if ( Environment::isMasterRank() )
    //     BOOST_TEST_MESSAGE( "a10b(1)=" << a10beval );

    auto a12b = form1(_test=Wh, _pattern=size_type(Pattern::EXTENDED));
    a12b = integrate(_range=boundaryfaces(meshnd), _expr=id(p));
    a12b.close();
    auto a12beval = a12b(p);
    if ( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( "a12b(1)=" << a12beval );

    BOOST_MESSAGE( "test_form2_faces ends for dim=" << T::value);
}


BOOST_AUTO_TEST_SUITE_END()
