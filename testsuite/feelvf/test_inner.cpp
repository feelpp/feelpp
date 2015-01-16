/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-14

  Copyright (C) 2014 Feel++ Consortium

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
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_inner
#include <testsuite/testsuite.hpp>
#endif

#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
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
AboutData
makeAbout()
{
    AboutData about( "test_inner" ,
                     "test_inner" ,
                     "0.2",
                     "nD(n=2,3) test inner",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
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

        auto g=soption(_name="functions.g");
        u.on(_range=elements(mesh), _expr=expr( g ));

        boost::mpi::timer ti;
        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "Check integral" );

        ti.restart();
        double a = integrate( _range= elements( mesh ), _expr=norm2(gradv(u)) ).evaluate()( 0,0 );

        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "time norm2(u)=" << ti.elapsed() << "s" );

        ti.restart();
        double b = integrate( _range= elements( mesh ), _expr= sqrt(inner(gradv(u),gradv(u))) ).evaluate()( 0,0 );

        if ( Environment::rank() == 0 )
            BOOST_TEST_MESSAGE( "time sqrt(inner(u,u))=" << ti.elapsed() << "s" );
        BOOST_CHECK_CLOSE( a, b, 1e-4 );
    }
};


#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() );
BOOST_AUTO_TEST_SUITE( inner_suite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    auto mesh = unitSquare();
    double a = integrate( _range= elements( mesh ), _expr= norm2(vec(cst(1.),cst(1.))) ).evaluate()( 0,0 );
    double b = integrate( _range= elements( mesh ), _expr= sqrt(inner(vec(cst(1.),cst(1.)),vec(cst(1.),cst(1.)) ) ) ).evaluate()( 0,0 );
    BOOST_CHECK_CLOSE( a, sqrt(2.), 1e-12);
    BOOST_CHECK_CLOSE( b, sqrt(2.), 1e-12);
    BOOST_CHECK_CLOSE( a, b, 1e-13 );

}

BOOST_AUTO_TEST_CASE( test_1 )
{
    Test<3> test;
    test.run();
}

BOOST_AUTO_TEST_SUITE_END()
#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=feel_options() );
    Test<2>  test;
    test.run();
}
#endif





