/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Thomas Lantz
       Date: 2015-04-27

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
#include <feel/feelcore/testsuite.hpp>
#endif


#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelpoly/multiscalequadrature.hpp>
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
                     "test integrate Quadra",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Thomas Lantz", "student", "", "" );
    return about;
}

class Test
{
 public :
       
    void run(std::string s)

        {
            /// [mesh] 
            auto mesh = createGMSHMesh( _mesh=new Mesh<Hypercube<2>>,  
                                        _desc=domain(_name="polymere",
                                                     _xmax=1,
                                                     _ymax=1));

  
            /// [expression]
            // our function to integrate
            auto g = expr( s  );

            /// [integrals]
            // compute on \Omega
            auto intf_1 = integrate( _range = elements( mesh ),
                                     _expr = g,
                                     _quad=_Q<1,MultiScaleQuadrature>() ).evaluate();
            auto intf_12 = integrate( _range = elements( mesh ),
                                      _expr = g).evaluate();
            // compute on boundary
            auto intf_2 = integrate( _range = boundaryfaces( mesh ),
                                     _expr = g,
                                     _quad=_Q<1,MultiScaleQuadrature>()  ).evaluate();
            auto intf_22 = integrate( _range = boundaryfaces( mesh ),
                                      _expr = g).evaluate();
    
            // compute integral of grad f
            auto grad_g = grad<2>(g);
            auto intgrad_f = integrate( _range = elements( mesh ),
                                        _expr = grad_g,
                                        _quad=_Q<1,MultiScaleQuadrature>()  ).evaluate();
            auto intgrad_f2 = integrate( _range = elements( mesh ),
                                         _expr = grad_g).evaluate();

            // values view    
            std::cout << "int_Omega " << g << " = " << intf_1  << std::endl
                      << "int_{boundary of Omega} " << g << " = " << intf_2 << std::endl
                      << "int_Omega grad " << g << " = "
                      << "int_Omega  " << grad_g << " = "
                      << intgrad_f  << std::endl;

            BOOST_CHECK_CLOSE( intf_1(0,0), intf_12(0,0), 5 );
            BOOST_CHECK_CLOSE( intf_2(0,0), intf_22(0,0), 5 );

            std::cout <<"" << std::endl;

        }
    
    void resol(std::string s)

        {
            /// [mesh] 
            auto mesh = createGMSHMesh( _mesh=new Mesh<Hypercube<2>>,  
                                        _desc=domain(_name="polymere",
                                                     _xmax=1,
                                                     _ymax=1));

            auto Vh = Pch<1>( mesh );
            auto u=Vh->element();
            auto v=Vh->element();
         
            /// [expression]
            // our function to integrate
            auto g = expr( s );
            auto gProj = vf::project( _space=Vh, _range=elements( mesh ), _expr=g);

            auto a = form2( _trial=Vh, _test=Vh );
            a=integrate( _range=elements( mesh ),
                         _expr=idt(u)*id(v),
                         _quad=_Q<1,MultiScaleQuadrature>() );

            auto l = form1( _test=Vh );
            l= integrate( _range=elements( mesh ),
                          _expr=g*id(v),
                          _quad=_Q<1,MultiScaleQuadrature>() ); 
       
            a.solve( _rhs=l, _solution=u );

            std::cout << "|| u-f ||^2 =" << normL2 ( _range=elements( mesh ), _expr=idv(u)-idv(gProj)) << std::endl;

            std::cout <<"" << std::endl;


            //        BOOST_CHECK_CLOSE( idv(u), idv(gProj), 1e-2 );
        }


};

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() )
BOOST_AUTO_TEST_SUITE( integrQuadra_suite )

BOOST_AUTO_TEST_CASE( test_run0 )
{
    Test t0 ;
    t0.run("x:x:y");
}


BOOST_AUTO_TEST_CASE( test_run1 ) 
{
    Test t1 ;
    t1.run("x*x:x:y");
}


BOOST_AUTO_TEST_CASE( test_run2 )
{
    Test t2 ;
    t2.run("x*x*x*x*x:x:y");
}


BOOST_AUTO_TEST_CASE( test_run3 )
{
    Test t3 ;
    t3.run("x+y:x:y");
}


BOOST_AUTO_TEST_CASE( test_run4 )
{
    Test t4 ;
    t4.run("cos(x)*sin(y):x:y");
}


BOOST_AUTO_TEST_CASE( test_run5 )
{
    Test t5 ;
    t5.run("cos(x*y):x:y");
}


BOOST_AUTO_TEST_CASE( test_run6 )
{
    Test t6 ;
    t6.run("tan(x*x*x):x:y");
}


BOOST_AUTO_TEST_CASE( test_run7 )
{
    Test t7 ;
    t7.run("y*exp(x):x:y");
}


BOOST_AUTO_TEST_CASE( test_resol0 )
{
    Test t0 ;
    t0.resol("sin(x):x:y");
}


BOOST_AUTO_TEST_CASE( test_resol1 )
{
    Test t1 ;
    t1.resol("x*x*x*x:x:y");
}


BOOST_AUTO_TEST_CASE( test_resol2 )
{
    Test t2 ;
    t2.resol("x+y:x:y");
}

BOOST_AUTO_TEST_CASE( test_resol3 )
{
    Test t3 ;
    t3.resol("x*y:x:y");
}


BOOST_AUTO_TEST_CASE( test_resol4 )
{
    Test t4 ;
    t4.resol("cos(x)*sin(y):x:y");
}

BOOST_AUTO_TEST_CASE( test_resol5 )
{
    Test t5 ;
    t5.resol("cos(x*y):x:y");
}

BOOST_AUTO_TEST_CASE( test_resol6 )
{
    Test t6 ;
    t6.resol("tan(x*x*x):x:y");
}


BOOST_AUTO_TEST_CASE( test_resol7 )
{
    Test t7 ;
    t7.resol("y*exp(x):x:y");
}


BOOST_AUTO_TEST_SUITE_END()
#else
std::cout << "USE_BOOST_TEST non define" << std::endl;
#endif

