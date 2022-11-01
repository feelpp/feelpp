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
#define BOOST_TEST_MODULE test_laplacian
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/thch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/loadmesh.hpp>
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

template<int Dim, int O = 1, template<uint16_type,int,uint16_type> class C = Simplex>
class Test:
    public Simget
{
public :

    void run()
    {
        auto mesh = createGMSHMesh( _mesh=new Mesh<C<Dim,O,Dim>>,
                                    _desc=domain( _name=( boost::format( "%1%-%2%-%3%-%4%" ) % soption(_name="gmsh.domain.shape") % Dim % O  % C<Dim,O,Dim>::type() ).str() ,
                                                  _dim=Dim, _order=O, _convex=C<Dim,O,Dim>::type()  ) );

        auto Xh = Pch<3>( mesh );
        auto u = Xh->element();
        auto v = Xh->element();

        auto g=soption(_name="f1");
        auto lapg = laplacian(expr(g));
        u.on(_range=elements(mesh), _expr=expr( g ));

        boost::mpi::timer ti;
        if ( Environment::isMasterRank()  )
            BOOST_TEST_MESSAGE( "Check integral with mesh dim " << Dim << " order " << O << " convex " << (C<Dim,O,Dim>::type()) );

        ti.restart();
        auto c = normL2( _range=elements(mesh), _expr=laplacianv(u)-lapg, _quad=_Q<4>() );
        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "  . [time normL2 laplacianv(u)-" << lapg <<" =" << ti.elapsed() << "s] c=" <<  c );

        BOOST_CHECK_SMALL( c, 1e-9 );
    }
};


template<int Dim, int O = 1, template<uint16_type,int,uint16_type> class C = Simplex>
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
        auto c = normL2( _range=elements(mesh), _expr=laplacianv(u)-lapg );

        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "  . [time normL2 laplacianv(u)-" << lapg <<" =" << ti.elapsed() << "s] c=" <<  c );

        BOOST_CHECK_SMALL( c, 1e-9 );
    }
};


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( inner_suite )

#if 1
BOOST_AUTO_TEST_CASE( test_21_s )
{
    Test<2,1> test;
    test.run();
    BOOST_TEST_MESSAGE( "test_21" << tc::green << " OK" << tc::reset);
}

BOOST_AUTO_TEST_CASE( test_22_s )
{
    Test<2,2> test;
    test.run();
    BOOST_TEST_MESSAGE( "test_22_s" << tc::green << " OK" << tc::reset);
}

BOOST_AUTO_TEST_CASE( test_31_s )
{
    Test<3,1> test;
    test.run();
    BOOST_TEST_MESSAGE( "test_31_s" << tc::green << " OK" << tc::reset);
}

BOOST_AUTO_TEST_CASE( test_21_h )
{
    Test<2,1,Hypercube> test;
    test.run();
    BOOST_TEST_MESSAGE( "test_21_h" << tc::green << " OK" << tc::reset);
}

BOOST_AUTO_TEST_CASE( test_22_h )
{
    Test<2,2,Hypercube> test;
    test.run();
    BOOST_TEST_MESSAGE( "test_22_h" << tc::green << " OK" << tc::reset);
}

BOOST_AUTO_TEST_CASE( test_31_h )
{
    Test<3,1,Hypercube> test;
    test.run();
    BOOST_TEST_MESSAGE( "test_31_h" << tc::green << " OK" << tc::reset);
}


// Vectorial

BOOST_AUTO_TEST_CASE( testv_21_s )
{
    TestV<2,1> testv;
    testv.run();
    BOOST_TEST_MESSAGE( "testv_21" << tc::green << " OK" << tc::reset);
}

BOOST_AUTO_TEST_CASE( testv_22_s )
{
    TestV<2,2> testv;
    testv.run();
    BOOST_TEST_MESSAGE( "testv_22_s" << tc::green << " OK" << tc::reset);
}

BOOST_AUTO_TEST_CASE( testv_31_s )
{
    TestV<3,1> testv;
    testv.run();
    BOOST_TEST_MESSAGE( "testv_31_s" << tc::green << " OK" << tc::reset);
}

BOOST_AUTO_TEST_CASE( testv_21_h )
{
    TestV<2,1,Hypercube> testv;
    testv.run();
    BOOST_TEST_MESSAGE( "testv_21_h" << tc::green << " OK" << tc::reset);
}

BOOST_AUTO_TEST_CASE( testv_22_h )
{
    TestV<2,2,Hypercube> testv;
    testv.run();
    BOOST_TEST_MESSAGE( "testv_22_h" << tc::green << " OK" << tc::reset);
}

BOOST_AUTO_TEST_CASE( testv_31_h )
{
    TestV<3,1,Hypercube> testv;
    testv.run();
    BOOST_TEST_MESSAGE( "testv_31_h" << tc::green << " OK" << tc::reset);
}


#else
BOOST_AUTO_TEST_CASE( test_21_g2 )
{
    Test<2,1,Hypercube> test;
    test.run();
}

#endif
BOOST_AUTO_TEST_SUITE_END()
