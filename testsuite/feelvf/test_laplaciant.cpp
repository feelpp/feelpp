/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2014-01-14

 Copyright (C) 2014-2015 Feel++ Consortium

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
#define USE_BOOST_TEST 0
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_laplacian
#include <testsuite/testsuite.hpp>
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
#include <feel/feelvf/trans.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace vf;

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
                     "Copyright (c) 2014 Feel++ Consortium" );

    about.addAuthor( "JB Wahl", "developer", "wahl.jb@gmail.com", "" );
    return about;
}

template<int Dim>
class Test:
    public Simget
{
public :

    void run()
        {
            auto mesh = createGMSHMesh( _mesh=new Mesh<Simplex<Dim,1>>,
                                        _desc=domain( _name=( boost::format( "%1%-%2%" ) % soption(_name="gmsh.domain.shape") % Dim ).str() ,
                                                      _dim=Dim ) );

            auto Xh = Pch<2>( mesh );
            auto u = Xh->element();
            auto v = Xh->element();

            auto f=soption(_name="f1");
            auto lapf = laplacian(expr(f));
            u.on(_range=elements(mesh), _expr=expr( f ));
            v.on(_range=elements(mesh), _expr=cst(1.) );

            boost::mpi::timer ti;
            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "Scalar Laplacian "<<Dim<<"D" );

            ti.restart();
            auto f2 = form2( _test=Xh, _trial=Xh );
            f2 = integrate( _range=elements(mesh), _expr=laplaciant(u)*id(u) );
            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "[time assemble laplaciant(u)*id(u)=" << ti.elapsed() << "s]" );

            ti.restart();
            auto a = f2.matrixPtr()->energy(v,u);
            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "[time compute energy  =" << ti.elapsed() << "s] a=" <<  a );

            ti.restart();
            auto b = integrate( _range= elements( mesh ), _expr= lapf ).evaluate()(0,0);
            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "[time int " << lapf <<" =" << ti.elapsed() << "s] b=" <<  b );

            BOOST_CHECK_SMALL( a-b, 1e-9 );
        }
};

template<int Dim>
class TestV:
    public Simget
{
public :

    void run()
        {
            auto mesh = createGMSHMesh( _mesh=new Mesh<Simplex<Dim,1>>,
                                        _desc=domain( _name=( boost::format( "%1%-%2%" ) % soption(_name="gmsh.domain.shape") % Dim ).str() ,
                                                      _dim=Dim ) );

            auto Xh = Pchv<2>( mesh );
            auto u = Xh->element();
            auto v = Xh->element();

            auto f=soption( _name=( boost::format("f%1%") %Dim ).str()   );
            auto lapf = laplacian(expr<Dim,1>(f));
            u.on(_range=elements(mesh), _expr=expr<Dim,1>( f ) );
            v.on( _range=elements(mesh), _expr=one() );

            boost::mpi::timer ti;
            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "Vectorial Laplacian "<<Dim<<"D" );

            ti.restart();
            auto f2 = form2( _test=Xh, _trial=Xh );
            f2 = integrate( _range=elements(mesh), _expr=trans(laplaciant(u))*id(u) );
            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "[time assemble laplaciant(u)*id(u)=" << ti.elapsed() << "s]" );

            ti.restart();
            auto a = f2.matrixPtr()->energy(v,u);
            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "[time compute  energy =" << ti.elapsed() << "s] a=" <<  a );

            ti.restart();
            auto b = integrate( _range= elements( mesh ), _expr= inner(lapf,one()) ).evaluate()(0,0);
            if ( Environment::rank() == 0 )
                BOOST_TEST_MESSAGE( "[time int " << lapf <<" =" << ti.elapsed() << "s] b=" <<  b );

            BOOST_CHECK_SMALL( a-b, 1e-9 );
        }
};


#if defined(USE_BOOST_TEST)
 FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
 BOOST_AUTO_TEST_SUITE( inner_suite )


 BOOST_AUTO_TEST_CASE( test_21 )
 {
     Test<2> test;
     test.run();
 }

BOOST_AUTO_TEST_CASE( test_31 )
{
    Test<3> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_22 )
{
    TestV<2> test;
    test.run();
}

BOOST_AUTO_TEST_CASE( test_33 )
{
    TestV<3> test;
    test.run();
}

BOOST_AUTO_TEST_SUITE_END()
#else

int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=makeOptions() );

    Test<2>  test2s;
    test2s.run();

    Test<3>  test3s;
    test3s.run();

    TestV<2>  test2v;
    test2v.run();

    TestV<3>  test3v;
    test3v.run();

}
#endif
