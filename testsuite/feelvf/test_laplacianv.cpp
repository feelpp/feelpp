/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-14

  Copyright (C) 2014-2016 Feel++ Consortium

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define BOOST_TEST_MODULE test_laplacianv
#include <testsuite/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/norm2.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/matvec.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Laplacian Options" );
    options.add_options()
        ( "f1",po::value<std::string>()->default_value( "" ),"test function 1D" )
        ( "f2",po::value<std::string>()->default_value( "" ),"test function 2D" )
        ( "f3",po::value<std::string>()->default_value( "" ),"test function 3D" )
        ;
    options.add( feel_options() );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_laplacianv" ,
                     "test_laplacianv" ,
                     "0.2",
                     "nD(n=2,3) test laplacianv",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

template<int Dim, int O = 1, template<uint16_type,uint16_type,uint16_type> class C = Simplex>
class Test:
    public Simget
{
public :

    void run()
    {
        auto mesh = createGMSHMesh( _mesh=new Mesh<C<Dim,O,Dim>>,
                                    _desc=domain( _name=( boost::format( "%1%-%2%-%3%-%4%" ) % soption(_name="gmsh.domain.shape") % Dim % O  % C<Dim,O,Dim>::type() ).str() ,
                                                  _dim=Dim, _order=O, _convex=C<Dim,O,Dim>::type()  ) );

        auto Xh = Pch<2>( mesh );
        auto u = Xh->element();
        auto v = Xh->element();

        auto g=soption(_name="f1");
        auto lapg = laplacian(expr(g));
        u.on(_range=elements(mesh), _expr=expr( g ));

        boost::mpi::timer ti;
        if ( Environment::isMasterRank()  )
            BOOST_TEST_MESSAGE( "Check integral with mesh dim " << Dim << " order " << O << " convex " << (C<Dim,O,Dim>::type()) );

        ti.restart();
        auto a = integrate( _range= elements( mesh ), _expr=laplacianv(u) ).evaluate();

        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "  . [time int laplacian(u)=" << ti.elapsed() << "s] a=" << a );

        ti.restart();
        auto b = integrate( _range= elements( mesh ), _expr= lapg ).evaluate();

        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "  . [time int " << lapg <<" =" << ti.elapsed() << "s] b=" <<  b );
        BOOST_CHECK_SMALL( (a-b).norm(), 1e-10 );
    }
};


template<int Dim, int O = 1, template<uint16_type,uint16_type,uint16_type> class C = Simplex>
class TestV:
    public Simget
{
public :

    void run()
    {
        auto mesh = createGMSHMesh( _mesh=new Mesh<C<Dim,O,Dim>>,
                                    _desc=domain( _name=( boost::format( "%1%-%2%-%3%-%4%" ) % soption(_name="gmsh.domain.shape") % Dim % O % C<Dim,O,Dim>::type() ).str() ,
                                                  _dim=Dim, _order=O, _convex=(C<Dim,O,Dim>::type() )) );

        auto Xh = Pchv<2>( mesh );
        auto u = Xh->element();
        auto v = Xh->element();

        auto g=soption( _name=( boost::format("f%1%") %Dim ).str() );
        auto lapg = laplacian(expr<Dim,1>(g));
        u.on(_range=elements(mesh), _expr=expr<Dim,1>( g ));

        boost::mpi::timer ti;
        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "Check integral mesh dim " << Dim << " order " << O << " convex " << (C<Dim,O,Dim>::type())  );

        ti.restart();
        auto a = integrate( _range= elements( mesh ), _expr=laplacianv(u) ).evaluate();

        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "  . [time int laplacian(" << g << ")=" << ti.elapsed() << "s] a=" << a );

        ti.restart();
        auto b = integrate( _range= elements( mesh ), _expr= lapg ).evaluate();

        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "  . [time int " << lapg <<" =" << ti.elapsed() << "s] b=" <<  b );
        BOOST_CHECK_SMALL( (a-b).norm(), 1e-10 );
    }
};

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( inner_suite )

#if 1
BOOST_AUTO_TEST_CASE( test_21 )
{
    Test<2> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_21_h_q1 )
{
    Test<2,1,Hypercube> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_21_g2 )
{
    Test<2,2> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_22 )
{
    TestV<2> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_22_g2 )
{
    TestV<2,2> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_22_g2_q2 )
{
    TestV<2,2,Hypercube> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_31 )
{
    Test<3> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_31_q1 )
{
    Test<3,1,Hypercube> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_33 )
{
    TestV<3> test;
    test.run();
}
BOOST_AUTO_TEST_CASE( test_33_q1 )
{
    TestV<3,1,Hypercube> test;
    test.run();
}
#else
BOOST_AUTO_TEST_CASE( test_21_g2 )
{
    Test<2,1,Hypercube> test;
    test.run();
}

#endif
BOOST_AUTO_TEST_SUITE_END()
