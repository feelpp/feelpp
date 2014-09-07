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
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

/** use Feel namespace */
using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description testfspacecontext( "FESpace context options" );
    testfspacecontext.add_options()
        //( "mesh2d.hsize", po::value<double>()->default_value( 0.2 ), "mesh size" )
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


Eigen::VectorXd
r( std::vector<node_type> vec )
{
    int size = vec.size();
    Eigen::VectorXd res( vec.size() );
    for(int i=0; i<size; i++)
    {
        double x = vec[i](0);
        double y = vec[i](1);
        res(i)=sqrt(x*x+y*y);
    }
    return res;
}

Eigen::VectorXd
theta( std::vector<node_type> vec )
{
    int size = vec.size();
    Eigen::VectorXd res( vec.size() );
    for(int i=0; i<size; i++)
    {
        double x = vec[i](0);
        double y = vec[i](1);
        res(i)=2*atan(y/(x+sqrt(x*x+y*y)));
    }
    return res;
}

Eigen::VectorXd
x( std::vector<node_type> vec )
{
    int size = vec.size();
    Eigen::VectorXd res( vec.size() );
    for(int i=0; i<size; i++)
        res(i)=vec[i](0);
    return res;
}

Eigen::VectorXd
y( std::vector<node_type> vec )
{
    int size = vec.size();
    Eigen::VectorXd res( vec.size() );
    for(int i=0; i<size; i++)
        res(i)=vec[i](1);
    return res;
}


Eigen::VectorXd
x2( std::vector<node_type> vec )
{
    int size = vec.size();
    Eigen::VectorXd res( vec.size() );
    for(int i=0; i<size; i++)
        res(i)=vec[i](0)*vec[i](0);
    return res;
}

Eigen::VectorXd
x3( std::vector<node_type> vec )
{
    int size = vec.size();
    Eigen::VectorXd res( vec.size() );
    for(int i=0; i<size; i++)
        res(i)=vec[i](0)*vec[i](0)*vec[i](0);
    return res;
}

Eigen::VectorXd
x4( std::vector<node_type> vec )
{
    int size = vec.size();
    Eigen::VectorXd res( vec.size() );
    for(int i=0; i<size; i++)
        res(i)=vec[i](0)*vec[i](0)*vec[i](0)*vec[i](0);
    return res;
}

Eigen::VectorXd
xy( std::vector<node_type> vec )
{
    int size = vec.size();
    Eigen::VectorXd res( vec.size() );
    for(int i=0; i<size; i++)
        res(i)=vec[i](0)*vec[i](1);
    return res;
}

Eigen::VectorXd
x2y2( std::vector<node_type> vec )
{
    int size = vec.size();
    Eigen::VectorXd res( vec.size() );
    for(int i=0; i<size; i++)
        res(i)=vec[i](0)*vec[i](0)*vec[i](1)*vec[i](1);
    return res;
}

Eigen::VectorXd
xy3( std::vector<node_type> vec )
{
    int size = vec.size();
    Eigen::VectorXd res( vec.size() );
    for(int i=0; i<size; i++)
        res(i)=vec[i](0)*vec[i](1)*vec[i](1)*vec[i](1);
    return res;
}

Eigen::VectorXd
sin2pix( std::vector<node_type> vec )
{
    int size = vec.size();
    Eigen::VectorXd res( vec.size() );
    for(int i=0; i<size; i++)
    {
        double x = vec[i](0);
        res(i) = sin( 2*pi*x );
    }
    return res;
}



template<int Dim, int Order>
void
testFspaceContext()
{
    //auto mesh=unitHypercube<Dim>();
    auto mesh=unitSquare();
    auto Xh = Pch<Order>( mesh );
    LOG(INFO)<<"nDof : "<<Xh->nDof();
    auto ctx = Xh->context();
    BOOST_TEST_MESSAGE( "functionspace defined\n" );

    //expression we want to evaluate
    auto exprX = Px();
    auto exprX2 = Px()*Px();
    auto exprX3 = Px()*Px()*Px();
    auto exprX4 = Px()*Px()*Px()*Px();
    auto exprXY = Px()*Py();
    auto exprX2Y2 = Px()*Px()*Py()*Py();
    auto exprXY3 = Px()*Py()*Py()*Py();
    auto exprY = Py();
    auto exprTheta = 2*atan(Py()/(Px()+sqrt(Px()*Px()+Py()*Py())));
    auto exprR= sqrt(Px()*Px()+Py()*Py());
    auto exprSin2PiX= sin( 2 * pi * Px() );


    //projection on the mesh
    auto px = vf::project( Xh , elements(mesh), exprX );
    auto px2 = vf::project( Xh , elements(mesh), exprX2 );
    auto px3 = vf::project( Xh , elements(mesh), exprX3 );
    auto px4 = vf::project( Xh , elements(mesh), exprX4 );
    auto pxy = vf::project( Xh , elements(mesh), exprXY );
    auto px2y2 = vf::project( Xh , elements(mesh), exprX2Y2 );
    auto pxy3 = vf::project( Xh , elements(mesh), exprXY3 );
    auto py = vf::project( Xh , elements(mesh), exprY );
    auto ptheta = vf::project( Xh , elements(mesh), exprTheta );
    auto pr = vf::project( Xh , elements(mesh), exprR );
    auto psin2pix = vf::project( Xh , elements(mesh), exprSin2PiX );

    BOOST_TEST_MESSAGE( "expression defined done\n" );


    std::vector< node_type > vec_t;
    node_type t1(Dim), t2(Dim), t3(Dim);
    /*x*/ t1(0)=0.1; /*y*/ t1(1)=0.2;
    /*x*/ t2(0)=0.1; /*y*/ t2(1)=0.8;
    /*x*/ t3(0)=1; /*y*/ t3(1)=1;
    vec_t.push_back(t1);
    vec_t.push_back(t2);
    vec_t.push_back(t3);
    BOOST_TEST_MESSAGE( "define pts done\n" );

    ctx.add( t1 );
    ctx.add( t2 );
    ctx.add( t3 );
    BOOST_TEST_MESSAGE( "pts added to ctx\n" );

    //evaluation on all nodes via functionspace and ctx
    auto px_evaluate = px.evaluate( ctx );
    auto py_evaluate = py.evaluate( ctx );
    auto ptheta_evaluate = ptheta.evaluate( ctx );
    auto pr_evaluate = pr.evaluate( ctx );


    //boost::timer t;
    auto evaluateX = evaluateFromContext( _context=ctx, _expr= exprX );
    auto evaluateX2 = evaluateFromContext( _context=ctx, _expr= exprX2 );

    auto evaluateX3 = evaluateFromContext( _context=ctx, _expr= exprX3 );
    auto evaluateX4 = evaluateFromContext( _context=ctx, _expr= exprX4 );
    auto evaluateXY = evaluateFromContext( _context=ctx, _expr= exprXY );
    auto evaluateX2Y2 = evaluateFromContext( _context=ctx, _expr= exprX2Y2 );
    auto evaluateXY3 = evaluateFromContext( _context=ctx, _expr= exprXY3 );
    auto evaluateY = evaluateFromContext( _context=ctx, _expr= exprY );
    auto evaluateTheta = evaluateFromContext( _context=ctx, _expr= exprTheta );
    auto evaluateR = evaluateFromContext( _context=ctx, _expr= exprR );
    auto evaluateSin2PiX = evaluateFromContext( _context=ctx, _expr= exprSin2PiX );


    auto evaluateProjX = evaluateFromContext( _context=ctx, _expr=idv(px) );
    auto evaluateProjX2 = evaluateFromContext( _context=ctx, _expr=idv(px2) );
    auto evaluateProjX3 = evaluateFromContext( _context=ctx, _expr=idv(px3) );
    auto evaluateProjX4 = evaluateFromContext( _context=ctx, _expr=idv(px4) );
    auto evaluateProjXY = evaluateFromContext( _context=ctx, _expr=idv(pxy) );
    auto evaluateProjX2Y2 = evaluateFromContext( _context=ctx, _expr=idv(px2y2) );
    auto evaluateProjXY3 = evaluateFromContext( _context=ctx, _expr=idv(pxy3) );
    auto evaluateProjY = evaluateFromContext( _context=ctx, _expr=idv(py) );
    auto evaluateProjTheta = evaluateFromContext( _context=ctx, _expr=idv(ptheta) );
    auto evaluateProjR = evaluateFromContext( _context=ctx, _expr=idv(pr) );
    auto evaluateProjSin2PiX = evaluateFromContext( _context=ctx, _expr=idv(psin2pix) );


    //true expressions evaluated at points vec_t
    Eigen::VectorXd true_x ( vec_t.size() );       true_x = x ( vec_t );
    Eigen::VectorXd true_x2 ( vec_t.size() );      true_x2 = x2 ( vec_t );
    Eigen::VectorXd true_x3 ( vec_t.size() );      true_x3 = x3 ( vec_t );
    Eigen::VectorXd true_x4 ( vec_t.size() );      true_x4 = x4 ( vec_t );
    Eigen::VectorXd true_xy ( vec_t.size() );      true_xy = xy ( vec_t );
    Eigen::VectorXd true_x2y2 ( vec_t.size() );    true_x2y2 = x2y2 ( vec_t );
    Eigen::VectorXd true_xy3 ( vec_t.size() );     true_xy3 = xy3 ( vec_t );
    Eigen::VectorXd true_y ( vec_t.size() );       true_y = y ( vec_t );
    Eigen::VectorXd true_theta ( vec_t.size() );   true_theta = theta( vec_t );
    Eigen::VectorXd true_r ( vec_t.size() );       true_r = r( vec_t );
    Eigen::VectorXd true_sin2pix ( vec_t.size() ); true_sin2pix = sin2pix( vec_t );

    //verification
    double eval1 = ptheta_evaluate(0);
    double eval2 = ptheta.evaluate( ctx , 0);
    BOOST_CHECK_SMALL( (eval1-eval2), 1e-13 );

    eval1 = ptheta_evaluate(2);
    eval2 = ptheta.evaluate( ctx , 2);
    BOOST_CHECK_SMALL( (eval1-eval2), 1e-13 );

    BOOST_CHECK_SMALL( (evaluateX-true_x).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateY-true_y).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateX2-true_x2).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateX3-true_x3).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateX4-true_x4).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateXY-true_xy).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateX2Y2-true_x2y2).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateXY3-true_xy3).norm(), 1e-13 );

    BOOST_CHECK_SMALL( (evaluateTheta-true_theta).norm(), 1e-7 );
    BOOST_CHECK_SMALL( (evaluateR-true_r).norm(), 1e-7 );
    BOOST_CHECK_SMALL( (evaluateSin2PiX-true_sin2pix).norm(), 5e-6 );

    BOOST_CHECK_SMALL( (evaluateX-evaluateProjX).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateY-evaluateProjY).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateX2-evaluateProjX2).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateX3-evaluateProjX3).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateX4-evaluateProjX4).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateXY-evaluateProjXY).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateX2Y2-evaluateProjX2Y2).norm(), 1e-13 );
    BOOST_CHECK_SMALL( (evaluateXY3-evaluateProjXY3).norm(), 1e-13 );

    BOOST_CHECK_SMALL( (evaluateTheta-evaluateProjTheta).norm(), 1e-5 );
    BOOST_CHECK_SMALL( (evaluateR-evaluateProjR).norm(), 1e-7 );
    BOOST_CHECK_SMALL( (evaluateSin2PiX-evaluateProjSin2PiX).norm(), 5e-6 );


    auto evaluateProjSin2PiX_ = evaluateFromContext( _context=ctx, _expr=exprSin2PiX , _projection=true);
    BOOST_CHECK_SMALL( (evaluateProjSin2PiX_-evaluateProjSin2PiX).norm(), 1e-13 );

    //now vector field

    auto Xhv = Pchv<Order>( mesh );
    auto ctxv=Xhv->context();
    ctxv.add( t1 );
    ctxv.add( t2 );
    ctxv.add( t3 );

    auto vector1 = vf::project( Xhv , elements(mesh), vec( sin(Px()) , cos(Py()) ) );
    auto EvaluateProjVector1 = evaluateFromContext( _context=ctxv, _expr=idv(vector1) );
    auto EvaluateVector1 = evaluateFromContext( _context=ctxv, _expr=vec( sin(Px()) , cos(Py()) ) , _projection=true);
    BOOST_CHECK_SMALL( (EvaluateProjVector1-EvaluateVector1).norm() , 1e-13 );

} // TestFspaceContext ::run


/**
 * main code
 */

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_SUITE( fspace_context )

BOOST_AUTO_TEST_CASE( test_1 )
{
    testFspaceContext<2,5>();
}

BOOST_AUTO_TEST_SUITE_END()


