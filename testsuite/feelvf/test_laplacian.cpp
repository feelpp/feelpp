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
#define USE_BOOST_TEST 1
#if USE_BOOST_TEST
#define BOOST_TEST_MODULE test_laplacian
#include <feel/feelcore/testsuite.hpp>
#endif

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
#include <feel/feelvf/one.hpp>
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
    AboutData about( "test_laplacian" ,
                     "test_laplacian" ,
                     "0.2",
                     "nD(n=2,3) test laplacian",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2014-2016 Feel++ Consortium" );

    about.addAuthor( "JB Wahl", "developer", "wahl.jb@gmail.com", "" );
    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@cemosis.fr", "" );
    return about;
}

template<int Dim, typename ConvexT = Simplex<Dim>>
class Test:
    public Simget
{
public :

    void run()
        {
            auto mesh = createGMSHMesh( _mesh=new Mesh<ConvexT>,
                                        _desc=domain( _name=( boost::format( "%1%-%2%-%3%-%4%" ) % soption(_name="gmsh.domain.shape") % Dim % ConvexT::nOrder % ConvexT::type() ).str() ,
                                                      _dim=Dim,
                                                      _order=ConvexT::nOrder,
                                                      _convex=ConvexT::type() ) );

            auto Xh = Pch<2>( mesh );
            auto u = Xh->element();
            auto v = Xh->element();

            auto f=soption(_name="f1");
            auto lapf = laplacian(expr(f));
            u.on(_range=elements(mesh), _expr=expr( f ));

            boost::mpi::timer ti;

            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "Scalar Laplacian "<< ConvexT::name() );
            
            ti.restart();
            auto f1 = form1( _test=Xh );
            f1 = integrate( _range=elements(mesh), _expr=laplacian(u) );

            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "  . [time assemble laplacian(u)=" << ti.elapsed() << "s]" );


            ti.restart();
            auto a = f1(u);

            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "  . [time inner product  =" << ti.elapsed() << "s] a=" <<  a );


            ti.restart();
            auto b = integrate( _range= elements( mesh ), _expr= lapf ).evaluate()(0,0);

            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "  . [time int " << lapf <<" =" << ti.elapsed() << "s] b=" <<  b );


            BOOST_CHECK_SMALL( a-b, 1e-9 );
        }
};

template<int Dim, typename ConvexT=Simplex<Dim>>
class TestV:
    public Simget
{
public :

    void run()
        {
            auto mesh = createGMSHMesh( _mesh=new Mesh<ConvexT>,
                                        _desc=domain( _name=( boost::format( "%1%-%2%-%3%-%4%" ) % soption(_name="gmsh.domain.shape") % Dim % ConvexT::nOrder % ConvexT::type()  ).str() ,
                                                      _dim=Dim,
                                                      _order=ConvexT::nOrder,
                                                      _convex=ConvexT::type()) );

            auto Xh = Pchv<2>( mesh );
            auto u = Xh->element();
            auto v = Xh->element();

            auto f=soption( _name=( boost::format("f%1%") %Dim ).str()   );
            auto lapf = laplacian(expr<Dim,1>(f));
            u.on(_range=elements(mesh), _expr=expr<Dim,1>( f ));

            boost::mpi::timer ti;

            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "Vectorial Laplacian "<< ConvexT::name() );

            ti.restart();
            auto f1 = form1( _test=Xh );
            f1 = integrate( _range=elements(mesh), _expr=inner(laplacian(u),one()) );

            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "  . [time assemble laplacian(u)=" << ti.elapsed() << "s]" );

            ti.restart();
            auto a = f1(u);

            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "  . [time inner product  =" << ti.elapsed() << "s] a=" <<  a );

            ti.restart();
            auto b = integrate( _range= elements( mesh ), _expr= inner(lapf,one()) ).evaluate()(0,0);

            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "  . [time int " << lapf <<" =" << ti.elapsed() << "s] b=" <<  b );


            BOOST_CHECK_SMALL( a-b, 1e-9 );

        }
};


 FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
 BOOST_AUTO_TEST_SUITE( inner_suite )


 BOOST_AUTO_TEST_CASE( test_21 )
 {
     Test<2> test;
     test.run();
 }
BOOST_AUTO_TEST_CASE( test_22 )
{
    Test<2,Simplex<2,2>> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( testH_21 )
{
    Test<2,Hypercube<2>> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_31 )
{
    Test<3> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_32 )
{
    Test<3,Simplex<3,2>> test;
    test.run();
}
BOOST_AUTO_TEST_CASE( testH_31 )
{
    Test<3,Hypercube<3,1>> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( testv_21 )
{
    TestV<2> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( testv_22 )
{
    TestV<2,Simplex<2,2>> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( testHv_21 )
{
    TestV<2,Hypercube<2>> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( testv_31 )
{
    TestV<3> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( testv_32 )
{
    TestV<3,Simplex<3,2>> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( testHv_31 )
{
    TestV<3,Hypercube<3,1>> test;
    test.run();
}

BOOST_AUTO_TEST_SUITE_END()
