/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2013-03-14

  Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)

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
/**
   \file test_fspace_context.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2012-03-14
*/

#define BOOST_TEST_MODULE test_fspace_context
#include <testsuite/testsuite.hpp>

#include <fstream>

#include <feel/feel.hpp>

/** use Feel namespace */
using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description testfspacecontext( "FESpace context options" );
    testfspacecontext.add_options()
        ( "hsize", po::value<double>()->default_value( 0.2 ), "mesh size" )
        ( "shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)" )
        ;
    return testfspacecontext.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_fspace_context" ,
                     "test_fspace_context" ,
                     "0.2",
                     "nD(n=1,2,3) test context of functionspace",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}


template<int Dim>
void
testFspaceContext()
{
    //auto mesh=unitHypercube<Dim>();
    auto mesh=unitSquare();
    auto Xh = Pch<1>( mesh );
    auto ctx = Xh->context();
    BOOST_TEST_MESSAGE( "functionspace defined\n" );
    //expression we want to evaluate
    auto x = Px();
    auto y = Py();
    auto theta = 2*atan(Py()/(Px()+sqrt(Px()*Px()+Py()*Py())));
    auto r= sqrt(Px()*Px()+Py()*Py());

    //projection on the mesh
    auto px = vf::project( Xh , elements(mesh), x );
    auto py = vf::project( Xh , elements(mesh), y );
    auto ptheta = vf::project( Xh , elements(mesh), theta );
    auto pr = vf::project( Xh , elements(mesh), r );

    BOOST_TEST_MESSAGE( "expression defined done\n" );
    //nodes where we want evaluate expressions x,y,r,theta
#if 0
#if DIM==1
    node_type t1; /*x*/ t1(0)=0.1;
    ctx.add( t1 );
    node_type t2; /*x*/ t2(0)=0.2;
    ctx.add( t2 );
    node_type t3; /*x*/ t3(0)=1;
    ctx.add( t3 );
#elif DIM==2
    node_type t1; /*x*/ t1(0)=0.1; /*y*/ t1(1)=0.2;
    ctx.add( t1 );
    node_type t2; /*x*/ t2(0)=0.1; /*y*/ t2(1)=0.8;
    ctx.add( t2 );
    node_type t3; /*x*/ t3(0)=1; /*y*/ t3(1)=1;
    ctx.add( t3 );
#elif DIM==3
    node_type t1; /*x*/ t1(0)=0.1; /*y*/ t1(1)=0.2; /*z*/ t1(2)=1;
    ctx.add( t1 );
    node_type t2; /*x*/ t2(0)=0.1; /*y*/ t2(1)=0.8; /*z*/ t2(2)=0.4;
    ctx.add( t2 );
    node_type t3; /*x*/ t3(0)=1; /*y*/ t3(1)=1; /*z*/ t2(2)=0.6;
    ctx.add( t3 );
#else
    throw std::logic_error("ERROR with the dimension ( dim > 3 ) " );
#endif
#else

    node_type t1(Dim), t2(Dim), t3(Dim);
    if ( Dim >= 2 ) { /*x*/ t1(0)=0.1; /*y*/ t1(1)=0.2;}
    if ( Dim >= 2 ) { /*x*/ t2(0)=0.1; /*y*/ t2(1)=0.8; }
    if ( Dim >= 2 ) { /*x*/ t3(0)=1; /*y*/ t3(1)=1; }
    ctx.add( t1 );
    ctx.add( t2 );
    ctx.add( t3 );

#endif

    BOOST_TEST_MESSAGE( "define pts done\n" );

    //evaluation on all nodes via functionspace and ctx
    auto px_evaluate = px.evaluate( ctx );
    auto py_evaluate = py.evaluate( ctx );
    auto ptheta_evaluate = ptheta.evaluate( ctx );
    auto pr_evaluate = pr.evaluate( ctx );

    LOG(INFO) << "px_evaluate=" << px_evaluate << "\n";

    BOOST_TEST_MESSAGE( "evaluate expressions at pts done\n" );

    //true expressions (for verification)
    double x_t1=t1(0);
    double x_t2=t2(0);
    double x_t3=t3(0);

    double y_t1=t1(1);
    double y_t2=t2(1);
    double y_t3=t3(1);
    double r_t1 = sqrt( x_t1*x_t1 + y_t1*y_t1 );
    double r_t2 = sqrt( x_t2*x_t2 + y_t2*y_t2 );
    double r_t3 = sqrt( x_t3*x_t3 + y_t3*y_t3 );
    double theta_t1 = 2*atan(y_t1 / (x_t1 + r_t1 ) );
    double theta_t2 = 2*atan(y_t2 / (x_t2 + r_t2 ) );
    double theta_t3 = 2*atan(y_t3 / (x_t3 + r_t3 ) );


    //store true expressions in a vectore
    std::vector<double> solution_x, solution_y, solution_theta, solution_r;
    solution_x.resize(ctx.nPoints());
    solution_y.resize(ctx.nPoints());
    solution_theta.resize(ctx.nPoints());
    solution_r.resize(ctx.nPoints());

    //fill vectors solution
    solution_x[0] = x_t1;
    solution_x[1] = x_t2;
    solution_x[2] = x_t3;

    // dim >=2
    solution_y[0] = y_t1;
    solution_y[1] = y_t2;
    solution_y[2] = y_t3;
    solution_theta[0] = theta_t1;
    solution_theta[1] = theta_t2;
    solution_theta[2] = theta_t3;
    solution_r[0] = r_t1;
    solution_r[1] = r_t2;
    solution_r[2] = r_t3;

    boost::timer t;
    auto v1 = evaluateFromContext( _context=ctx, _expr=Px() );
    std::cout << "v1 = " << v1 << " time: " << t.elapsed() << "s\n";t.restart();
    auto v2 = evaluateFromContext( _context=ctx, _expr=idv(px) );
    std::cout << "v2 = " << v2 << " time: " << t.elapsed() << "s\n";
    BOOST_CHECK_SMALL( (v1-v2).norm(), 1e-13 );

    BOOST_TEST_MESSAGE( "start check\n" );
    //verification step
    for( int i=0; i<ctx.nPoints(); i++)
     {
         //check for expression x
         double evaluation_x = px_evaluate( i );
         double evaluation_x_node = px.evaluate(ctx , i);
         BOOST_CHECK_CLOSE( evaluation_x, evaluation_x_node, 1e-13 );
         BOOST_CHECK_CLOSE( evaluation_x, solution_x[i], 1e-13 );
         BOOST_CHECK_CLOSE( evaluation_x_node, solution_x[i], 1e-13 );


#if 1 //DIM >= 2
         //check for expression y
         double evaluation_y = py_evaluate( i );
         double evaluation_y_node = py.evaluate(ctx , i);
         BOOST_CHECK_CLOSE( evaluation_y, evaluation_y_node, 1e-13 );
         BOOST_CHECK_CLOSE( evaluation_y, solution_y[i], 1e-13 );

         //check for expression theta
         double evaluation_theta = ptheta_evaluate( i );
         double evaluation_theta_node = ptheta.evaluate(ctx , i);
         BOOST_CHECK_CLOSE( evaluation_theta, evaluation_theta_node, 1e-13 );
         BOOST_CHECK_CLOSE( evaluation_theta, solution_theta[i], 1e-13 );

         //check for expression r
         double evaluation_r = pr_evaluate( i );
         double evaluation_r_node = pr.evaluate(ctx , i);
         BOOST_CHECK_CLOSE( evaluation_r, evaluation_r_node, 1e-13 );
         BOOST_CHECK_CLOSE( evaluation_r, solution_r[i], 1e-13 );
#endif
     }


} // TestFspaceContext ::run


/**
 * main code
 */

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_SUITE( fspace_context )

BOOST_AUTO_TEST_CASE( test_1 )
{
    testFspaceContext<2>();

}

BOOST_AUTO_TEST_SUITE_END()


