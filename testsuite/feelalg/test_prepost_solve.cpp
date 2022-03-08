/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2016-01-08
 
 Copyright (C) 2016 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3.0 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE prepost solve  testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <feel/feelcore/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/exporter.hpp>

const double DEFAULT_MESH_SIZE=0.5;

using namespace Feel;



inline
Feel::po::options_description
makeOptions()
{
  Feel::po::options_description options( "Options" );
  options.add_options()
      ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "h value" )
    ;
  return options;
}

inline
Feel::AboutData
makeAbout()
{
  Feel::AboutData about( "test_prepost_solve" ,
                         "test_prepost_solve" ,
                         "0.1",
                         "pre/post solve tests",
                         Feel::AboutData::License_GPL,
                         "Copyright (C) 2016 Feel++ Consortium" );

  about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
  return about;

}


    
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( prepost_solve_suite )

BOOST_AUTO_TEST_CASE( test_prepostsolve )
{
    std::string f_str = soption("functions.f");
    std::string f_name = "f";
    BOOST_MESSAGE("test_prepostsolve, checking with " << f_str);
    
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = Xh->element();
    auto v = Xh->element();

    auto f = expr(f_str,f_name);
    auto gradf = grad<2>(f);
    auto lapf = laplacian(f);

    auto a = form2( _test=Xh, _trial = Xh);
    a = integrate(_range=elements(mesh), _expr=gradt(u)*trans(grad(u)));
    auto l = form1( _test=Xh );
    l= integrate(_range=elements(mesh), _expr=-lapf*id(v));
    l+=integrate(_range=boundaryfaces(mesh), _expr=gradf*N()*id(v));
    auto zeromean =  [&Xh,&u]( vector_ptrtype rhs, vector_ptrtype sol )
        {
            Feel::cout << " . Call postsolve: remove mean value to input vector\n";
            double m = mean( _range=elements(Xh->mesh()), _expr=idv(u))(0,0);
            Feel::cout << "[zeromean] MEAN1:" << m << std::endl;
            u.add( -m );
            m = mean( _range=elements(Xh->mesh()), _expr=idv(u))(0,0);
            Feel::cout << "[zeromean] MEAN2:" << m << std::endl;
        };
    a.solve( _solution=u, _rhs=l, _pre=zeromean, _post=zeromean);
    double m = mean( _range=elements(Xh->mesh()), _expr=idv(u))(0,0);
    BOOST_CHECK_SMALL( m, 1e-10 );
    BOOST_MESSAGE( "MEAN:" << m );
    double mex = mean( _range=elements(Xh->mesh()), _expr=f) (0,0);
    BOOST_MESSAGE( "MEAN exact:" << mex );
    double l2error = normL2( _range=elements(mesh), _expr=(idv(u)-m)-(f-mex) );

    Feel::cout << "||u-(f-mex)||=" << l2error << std::endl;

    v.on(_range=elements(mesh), _expr=f-mex);

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->add( "v", v );
    e->save();
}

BOOST_AUTO_TEST_CASE( test_prepost_nlsolve )
{
    BOOST_MESSAGE("test_prepost_nlsolve, checking with " << soption("functions.f"));
    auto f = expr(soption("functions.f"),"f");
    auto gradf = grad<2>(f);
    auto lapf = laplacian(f);
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = Xh->element();
    auto v = Xh->element();
    auto a = form2( _test=Xh, _trial = Xh);
    a = integrate(_range=elements(mesh), _expr=gradt(u)*trans(grad(u)));
    auto l = form1( _test=Xh );
    a.close();
    l.close();
    auto zeromean =  [&Xh,&v]( vector_ptrtype rhs, vector_ptrtype sol )
        {
            if ( Environment::isMasterRank() )
                Feel::cout << " . Call postsolve(nonlinear): remove mean value to input vector\n";
            v = *sol;
            v.close();
            double m = mean( _range=elements(Xh->mesh()), _expr=idv(v))(0,0);
            v.add( -m );
            v.close();
            *sol = v;
        };

    auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            if (!J) J = backend()->newMatrix( _test=Xh, _trial=Xh );
            auto j = form2( _test=Xh, _trial=Xh, _matrix=J );
            j = a ;
        };
    auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
        {
            auto u = Xh->element(X);
            auto r = form1( _test=Xh, _vector=R );
            r= integrate(_range=elements(mesh), _expr=lapf*id(v)+gradv(u)*trans(grad(v)));
            r+=integrate(_range=boundaryfaces(mesh), _expr=-gradf*N()*id(v));
        };
    u.zero();
    backend()->nlSolver()->residual = Residual;
    backend()->nlSolver()->jacobian = Jacobian;
    backend()->nlSolve( _solution=u, _pre=zeromean, _post=zeromean );

    double m = mean( _range=elements(Xh->mesh()), _expr=idv(u))(0,0);
    BOOST_CHECK_SMALL( m, 1e-10 );
    BOOST_MESSAGE( "MEAN(u) must be zero:" << m );
    double mex = mean( _range=elements(Xh->mesh()), _expr=f) (0,0);
    BOOST_MESSAGE( "MEAN exact:" << mex );
    double l2error = normL2( _range=elements(mesh), _expr=(idv(u)-m)-(f-mex) );

    Feel::cout << "||u-(f-mex)||=" << l2error << Feel::endl;

    v.on(_range=elements(mesh), _expr=f-mex);

    auto e = exporter( _mesh=mesh, _name="nlsolve" );
    e->add( "u", u );
    e->add( "v", v );
    e->save();
}


BOOST_AUTO_TEST_SUITE_END()

