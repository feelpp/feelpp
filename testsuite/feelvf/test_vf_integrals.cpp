/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 SPDX-FileCopyrightText: 2015-2026 University of Strasbourg
 SPDX-License-Identifier: LGPL-3.0-or-later
 */
#define BOOST_TEST_MODULE test_integrals
#include <feel/feelcore/testsuite.hpp>

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
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/norm2.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/one.hpp>
#include <feel/feelvf/matvec.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/cross.hpp>
#include <feel/feelvf/print.hpp>
#include <list>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Laplacian Options" );
    options.add_options()
        ( "N",po::value<int>()->default_value( 10 ),"number of integrals to compute" )
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
    AboutData about( "test_integrals" ,
                     "test_integrals" ,
                     "0.2",
                     "nD(n=2,3) test integrals",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );

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

            auto Xh = Pchv<2>( mesh );
            auto u = Xh->element();
            u.on(_range=elements(mesh),_expr=vec(cst(1.),cst(1.),cst(1.)));
            
            Eigen::Matrix<double,3,1> n;
            n << 1,2,3;
            
            std::vector<Eigen::Matrix<double,3,1>> x( ioption("N"),  n);


            tic();
            auto v = integrate(_range=elements(mesh),_expr=_e1v,_quad=_Q<1>()).evaluate( x );
            if ( Environment::numberOfProcessors() == 1 )
                std::for_each( v.begin(), v.end(), 
                               [&n]( Eigen::Matrix<double,3,1> const& nn ) 
                               { 
                                   BOOST_CHECK_SMALL( (n-nn).norm(), 1e-11); 
                               });
            toc("lambda integral _e1=vec(x,y,z)",Environment::logVerbosityLevel()>0);
            
            int i=0;
            for( auto e : x )
            {
                tic();
                auto vv = integrate(_range=elements(mesh),_expr=vec(cst(e(0)),cst(e(1)),cst(e(2))),_quad=_Q<1>()).evaluate();
                toc("integral vec(x,y,z)",Environment::logVerbosityLevel()>0);
                BOOST_CHECK_SMALL( (vv-v[i]).norm(), 1e-11);
                ++i;
            }

            tic();
            v = integrate(_range=elements(mesh),_expr=cross(idv(u),_e1v-P()),_quad=_Q<3>()).evaluate( x );
            toc("lambda integral u x (_1-X)", Environment::logVerbosityLevel()>0);
            i = 0;
            for( auto e : x )
            {
                tic();
                auto vv = integrate(_range=elements(mesh),_expr=cross(idv(u),vec(cst(e(0)),cst(e(1)),cst(e(2)))-P()),_quad=_Q<3>()).evaluate();
                toc("integral u x ((x,y,z)-X) ",Environment::logVerbosityLevel()>0);
                BOOST_CHECK_SMALL( (v[i]-vv).norm(), 1e-10 );
                ++i;
            }
            tic();
            auto mid = mat<3,3>( cst(1.), cst(0.), cst(0.),
                                 cst(0.), cst(1.), cst(0.),
                                 cst(0.), cst(0.), cst(1.) );
            v = integrate(_range=elements(mesh),_expr=cross(trans(trans(idv(u))*mid),_e1v-P()),_quad=_Q<3>()).evaluate( x );
            toc("lambda integral u x (_1-X)", Environment::logVerbosityLevel()>0);
            i = 0;
            for( auto e : x )
            {
                tic();
                auto vv = integrate(_range=elements(mesh),_expr=cross(idv(u),vec(cst(e(0)),cst(e(1)),cst(e(2)))-P()),_quad=_Q<3>()).evaluate();
                toc("integral u x ((x,y,z)-X) ",Environment::logVerbosityLevel()>0);
                BOOST_CHECK_SMALL( (v[i]-vv).norm(), 1e-10 );
                ++i;
            }
        }
};


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( inner_suite )


BOOST_AUTO_TEST_CASE( test_3 )
{
    Test<3> test;
    test.run();
}


BOOST_AUTO_TEST_CASE( runtime_geometry_quadrature )
{
    // Integrate a known degree-six moment with automatic runtime-geometry
    // quadrature, including face, range-list and empty-range paths.
    for ( int order : { 1, 2 } )
    {
        auto mesh = createGMSHMesh(
            _mesh = new Mesh<Simplex<2, Dynamic>>( RuntimeOrder{order} ),
            _desc = domain( _name="quadrature_square_g" + std::to_string( order ),
                            _shape="hypercube", _dim=2, _h=0.4,
                            _xmin=0., _xmax=1., _ymin=0., _ymax=1. ),
            _force_rebuild=true, _straighten=false );
        auto range = elements( mesh );
        auto moment = Px()*Px()*Px()*Px()*Px()*Px();
        BOOST_CHECK_SMALL( integrate( _range=range, _expr=moment ).evaluate()(0,0) - 1./7., 1e-12 );
        BOOST_CHECK_SMALL( integrate( _range=boundaryfaces( mesh ), _expr=moment ).evaluate()(0,0)
                              - 9./7., 1e-12 );
        std::list<decltype(range)> ranges{range};
        BOOST_CHECK_SMALL( integrate( _range=ranges, _expr=moment ).evaluate()(0,0) - 1./7., 1e-12 );
        BOOST_CHECK_SMALL( integrate( _range=range, _expr=moment, _quad=_Q<12>() ).evaluate()(0,0)
                              - 1./7., 1e-12 );
        auto emptyRange = markedelements( mesh, "absent-marker" );
        BOOST_CHECK_SMALL( integrate( _range=emptyRange, _expr=moment ).evaluate()(0,0), 1e-12 );
    }
}

BOOST_AUTO_TEST_CASE( tensor_product_derivative_quadrature )
{
    auto mesh = createGMSHMesh(
        _mesh = new Mesh<Hypercube<2, Dynamic>>( RuntimeOrder{1} ),
        _desc = domain( _name="quadrature_tensor_square", _convex="Hypercube",
                        _shape="hypercube", _dim=2, _h=0.4,
                        _xmin=0., _xmax=1., _ymin=0., _ymax=1. ),
        _force_rebuild=true );
    auto space = Pch<Dynamic>( mesh, RuntimeOrder{2} );
    auto u = space->element();
    u.on( _range=elements(mesh), _expr=Px()*Px()*Py()*Py() );
    // grad(x^2 y^2) has degree two in the undifferentiated coordinate.
    auto energy = integrate( _range=elements(mesh), _expr=inner(gradv(u)) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( energy - 8./15., 1e-12 );
}

BOOST_AUTO_TEST_SUITE_END()
