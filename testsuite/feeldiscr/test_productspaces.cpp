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
#include <feel/feelcore/testsuite.hpp>

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

     auto cp = hana::cartesian_product(hana::make_tuple(sp.tupleSpaces(),sp.tupleSpaces()));

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
     ex->add("u",U(0_c));
     ex->add("v",U(1_c));
     ex->save();
 }

BOOST_AUTO_TEST_CASE( test3 )
{
    using namespace Feel;
    using Feel::cout;
    using namespace boost::hana::literals;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);

    backend(_rebuild=true);
    int n = int(doption("parameters.n"));
    if ( n <= 0 || n >= 10 ) return;
    ProductSpace<Pch_ptrtype<Mesh<Simplex<2>>,1>, true> ps( n, mesh );
    BOOST_CHECK_EQUAL( ps.numberOfSpaces(), n );
    cout << tc::red << "number of spaces " << tc::reset << ps.numberOfSpaces() << std::endl;


    auto U = ps.element();
    auto u = U[0];
    auto b = blockform2( ps );
    auto l = blockform1( ps );
    std::vector<std::string> alphabet { "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega" };

    for( int i = 0; i < n; ++i )
    {
        b(i,i) += integrate( _range=elements(mesh), _expr=idt(u)*id(u));
        l(i) += integrate( _range=elements(mesh), _expr=expr(soption("functions."+alphabet[i]))*id(u));
    }
    b.solve( _rhs=l, _solution=U );
    auto ex = exporter(_mesh=mesh);
    for(int i = 0; i < n; ++i )
    {
        ex->add(alphabet[i],U[i]);
        U[i].printMatlab(alphabet[i]+".m");
    }
    ex->save();
}

BOOST_AUTO_TEST_CASE( test4 )
{
    using namespace Feel;
    using Feel::cout;
    using namespace boost::hana::literals;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);

    backend(_rebuild=true);
    int n = int(doption("parameters.n"));
    auto Xh = Pch<1>(mesh);
    auto Zh = Pch<1>(mesh);
    auto Yh = Pchv<3>(mesh);
    auto ps = std::make_shared<ProductSpace<decltype(Pch<2>(mesh)), true>>( n, mesh );
    auto p = product2( ps, Xh, Yh, Zh );

    BOOST_CHECK_EQUAL( p.numberOfSpaces(), n+3 );
    cout << tc::red << "number of spaces " << tc::reset << p.numberOfSpaces() << std::endl;

    auto U = p.element();

    std::vector<std::string> alphabet { "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega" };

    auto u = U(0_c);
    auto w = U(1_c);
    auto z = U(2_c);
    auto l = blockform1( p );
    auto b = blockform2( p );


    l(0_c) = integrate( _range=elements(mesh), _expr=expr(soption("functions."+alphabet[0]))*id(u));
    b(0_c,0_c) += integrate( _range=elements(mesh), _expr=idt(u)*id(u));
    l(1_c) = integrate( _range=elements(mesh), _expr=expr(soption("functions."+alphabet[1]))*trans(id(w))*one());
    b(1_c,1_c) += integrate( _range=elements(mesh), _expr=trans(idt(w))*id(w));
    l(2_c) = integrate( _range=elements(mesh), _expr=expr(soption("functions."+alphabet[0]))*id(z));
    b(2_c,2_c) += integrate( _range=elements(mesh), _expr=idt(z)*id(z));
    for( int i = 0; i < n; ++i )
    {
        auto v = U(3_c,i);

        b(3_c,3_c,i,i) += integrate( _range=elements(mesh), _expr=idt(v)*id(v));
        l(3_c,i) += integrate( _range=elements(mesh), _expr=expr(soption("functions."+alphabet[2+i]))*id(v));
    }
    b.solve( _rhs=l, _solution=U );
    auto ex = exporter(_mesh=mesh);
    ex->add(alphabet[0],U(0_c));
    ex->add(alphabet[1],U(1_c));
    ex->add(alphabet[2],U(2_c));
    for(int i = 0; i < ps->numberOfSpaces(); ++i )
    {
        ex->add(alphabet[i+2],U(3_c,i));
        //U[i].printMatlab(alphabet[i]+".m");
    }
    ex->save();

}

BOOST_AUTO_TEST_SUITE_END()
