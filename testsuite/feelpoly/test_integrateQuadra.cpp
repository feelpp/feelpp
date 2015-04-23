/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)

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
/// [all]
#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_integrateQuadra
#include <testsuite/testsuite.hpp>
#endif


#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelpoly/QuadraRect.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelfilters/exporter.hpp>
using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "test_integrateQuadra" ,
                     "test_integrateQuadra" ,
                     "0.2",
                     "nD(n=2,3) test integrate Quadra",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}




class Test
{
 public :
       
        void run(std::string s)

        {
    // Initialize Feel++ Environment
    /*Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="myintegrals" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );
    */
    /// [mesh]
    // create the mesh (specify the dimension of geometric entity)

    auto mesh = createGMSHMesh( _mesh=new Mesh<Hypercube<2>>,  
                                 _desc=domain(_name="polymere",
                                              _xmax=1,
                                              _ymax=1));

    /// [mesh]

    /// [expression]
    // our function to integrate
    auto g = expr( s  );
    /// [expression]

    /// [integrals]
    // compute integral of g (global contribution): \(\int_{\partial \Omega} g\)
    auto intf_1 = integrate( _range = elements( mesh ),
                                 _expr = g,
                                 _quad=_Q<1,QuadraRect>() ).evaluate();
    auto intf_12 = integrate( _range = elements( mesh ),
                                 _expr = g).evaluate();
    // compute integral g on boundary: \( \int_{\partial \Omega} g \)
    auto intf_2 = integrate( _range = boundaryfaces( mesh ),
                             _expr = g,
                             _quad=_Q<1,QuadraRect>()  ).evaluate();
    auto intf_22 = integrate( _range = boundaryfaces( mesh ),
                             _expr = g).evaluate();
    // compute integral of grad f (global contribution): \( \int_{\Omega} \nabla g \)
    auto grad_g = grad<2>(g);
    auto intgrad_f = integrate( _range = elements( mesh ),
                                _expr = grad_g,
                                _quad=_Q<1,QuadraRect>()  ).evaluate();
    auto intgrad_f2 = integrate( _range = elements( mesh ),
                                _expr = grad_g).evaluate();

    // only the process with rank 0 prints to the screen to avoid clutter
    //if ( Environment::isMasterRank() )
        std::cout << "int_Omega " << g << " = " << intf_1  << std::endl
                  << "int_{boundary of Omega} " << g << " = " << intf_2 << std::endl
                  << "int_Omega grad " << g << " = "
                  << "int_Omega  " << grad_g << " = "
                  << intgrad_f  << std::endl;
    /// [integrals]

        BOOST_CHECK_CLOSE( intf_1(0,0), intf_12(0,0), 1e-2 );
        BOOST_CHECK_CLOSE( intf_2(0,0), intf_22(0,0), 1e-2 );
        //BOOST_CHECK_CLOSE( intgrad_f, intgrad_f2, 1e-4 );
    }

};

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() );
BOOST_AUTO_TEST_SUITE( integrQuadra_suite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    Test t0 ;
    t0.run("x:x:y");
}

BOOST_AUTO_TEST_CASE( test_1 ) 
{
    Test t1 ;
    t1.run("x*x:x:y");
}

BOOST_AUTO_TEST_CASE( test_2 )
{
    Test t2 ;
    t2.run("x*x*x*x*x:x:y");
}

BOOST_AUTO_TEST_CASE( test_3 )
{
    Test t3 ;
    t3.run("x+y:x:y");
}
BOOST_AUTO_TEST_SUITE_END()
#else
#endif

