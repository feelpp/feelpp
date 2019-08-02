//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 31 May 2019
//! @copyright 2019 Feel++ Consortium
//!


#define BOOST_TEST_MODULE test_integration_on_manifolds
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsegment.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelpython/pyexpr.hpp>
#include <feel/feelfilters/geotool.hpp>

/** use Feel namespace */
using namespace Feel;

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options("test_expr options");
    options.add_options()
        ("pyexpr.filename", Feel::po::value<std::string>()->default_value( "tetra.py" ), "filename to handle sympy")
        ;
    return options.add( Feel::feel_options() );
}

inline AboutData
makeAbout()
{
    AboutData about( "test_integration_on_manifolds",
                     "test_integration_on_manifolds",
                     "0.1",
                     "test integration on manifolds",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2019 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( int_manifold )

template<int Dim = 3,int OrderGeo=1>
class TestIntegrate
{
public:
    using mesh_type=Mesh<Simplex<Dim,OrderGeo>>;
    using manifold_mesh_type=Mesh<Simplex<Dim-1,OrderGeo,Dim>>;
    TestIntegrate( std::string shape = "", double meshSize = 0.1 )
        {
            if constexpr( Dim == 3 )
            {
                if ( shape.empty() )
                {
                GeoTool::Node x1( 1, 1, 1);
                GeoTool::Node x2( 1, -1, -1 );
                GeoTool::Node x3( -1, 1, -1 );
                GeoTool::Node x4( -1, -1, 1 );
                GeoTool::Tetrahedron C( meshSize, "OMEGA", x1, x2, x3, x4 );
                mesh = C.createMesh( _mesh = new mesh_type,
                                     _name = "omega_" + mesh_type::shape_type::name(),
                                     _partitions = Environment::worldComm().localSize() );
                }
                if ( shape == "cylinder" )
                {
                    GeoTool::Node Centre(0,0,0);
                    GeoTool::Node Rayon( 0.5);
                    GeoTool::Node Dir(1,0,0);
                    GeoTool::Node Lg(1,0,0);
                    GeoTool::Cylindre C( meshSize,"Cyl",Centre,Dir,Rayon,Lg);
                    mesh = C.createMesh( _mesh = new mesh_type,
                                         _name = "omega_" + mesh_type::shape_type::name(),
                                         _partitions = Environment::worldComm().localSize() );
                }
                
                
            }
            else if constexpr( Dim == 2 )
            {
                if ( shape == "circle" )
                {
                    GeoTool::Node x1( 0, 0);
                    GeoTool::Node x2( 0, 1 );
                    GeoTool::Circle C( meshSize, "OMEGA", x1, x2 );
                    mesh = C.createMesh( _mesh = new mesh_type,
                                         _name = "omega_" + mesh_type::shape_type::name(),
                                         _partitions = Environment::worldComm().localSize() );
                }
            }
            else if constexpr ( Dim ==1 )
                mesh = unitSegment();
            smesh = createSubmesh( _mesh=mesh, _range=boundaryfaces(mesh ), _update=0 );            
        }
    template<typename ExprT>
    double evaluateOnManifold( ExprT const& e )
        {
            return integrate( _range=elements(smesh), _expr=e ).evaluate()( 0, 0 );
        }
    template<typename ExprT>
    double evaluateInDomain( ExprT const& e )
        {
            return integrate( _range=elements(mesh), _expr=e ).evaluate()( 0, 0 );
        }
    template<typename ExprT>
    double evaluateOnBoundary( ExprT const& e )
        {
            return integrate( _range=boundaryfaces(mesh), _expr=e ).evaluate()( 0, 0 );
        }
    
    
    double evaluateWithSympy( std::string const& e )
        {
            std::map<std::string,std::string> locals{{"dim",std::to_string(Dim)},{"e",e},{"expected_value","0."}};
            Feel::pyexprFromFile( Environment::expand(soption("pyexpr.filename")), locals );
            return std::stod( locals.at("expected_value") );
        }
    
    template<typename ExprT>
    void check( ExprT const& e, double expected_value, double tolerance = 1e-10 )
        {
            double v = integrate( _range=elements(mesh), _expr=e ).evaluate()( 0, 0 );
            BOOST_CHECK_CLOSE( v, expected_value, tolerance );
        }
    template<typename ExprT>
    void checkOnBoundary( ExprT const& e, double expected_value, double tolerance = 1e-10 )
        {
            double v = integrate( _range=boundaryfaces(mesh), _expr=e ).evaluate()( 0, 0 );
            BOOST_CHECK_CLOSE( v, expected_value, tolerance );
        }
    template<typename ExprT>
    void checkSmallOnBoundary( ExprT const& e, double tolerance = 1e-10 )
        {
            double v = integrate( _range=boundaryfaces(mesh), _expr=e ).evaluate()( 0, 0 );
            BOOST_CHECK_SMALL( v, tolerance );
        }
    std::shared_ptr<mesh_type> mesh;
    std::shared_ptr<manifold_mesh_type> smesh;
};

    
BOOST_AUTO_TEST_CASE( t0 )
{
    using namespace Feel::vf;
    TestIntegrate<3> t;

    double v1 = t.evaluateOnManifold( cst(1.) );
    double v2 = t.evaluateOnBoundary( cst(1.) );
    BOOST_CHECK_CLOSE( v1, math::sqrt(3.)*8, 1e-10 ); // area = sqrt(3)*(2*sqrt(2))^2
    BOOST_CHECK_CLOSE( v1, v2, 1e-10 );
}
BOOST_AUTO_TEST_CASE( t1 )
{
    using namespace Feel::vf;
    TestIntegrate<3,1> t("cylinder",5e-2);

    double v1 = t.evaluateOnManifold( cst(1.) );
    double v2 = t.evaluateOnBoundary( cst(1.) );
    BOOST_CHECK_CLOSE( v1, 2*M_PI*0.5*0.5+2*M_PI*0.5*1, 1e-1 ); // area = 2*pi*r^2 + 2*pi*r*h
    BOOST_CHECK_CLOSE( v1, v2, 1e-10 );
}
BOOST_AUTO_TEST_CASE( t2 )
{
    using namespace Feel::vf;
    TestIntegrate<2,1> t("circle",5e-2);

    double v1 = t.evaluateOnManifold( cst(1.) );
    double v2 = t.evaluateOnBoundary( cst(1.) );
    BOOST_CHECK_CLOSE( v1, 2*M_PI*1, 1e-1 ); // length = 2*pi*r
    BOOST_CHECK_CLOSE( v1, v2, 1e-10 );
}

#if 0
BOOST_AUTO_TEST_CASE( t1 )
{
    using namespace Feel::vf;
    TestIntegrate<3> t;
    auto e00 = cst( 1.0 ) * expr( "1." ) - expr( "2." );
    BOOST_CHECK_EQUAL( hasDynamicContext( e00 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e00 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e00 ) ), false );
    BOOST_CHECK_EQUAL( hasValue( staticContext( e00 ), vm::POINT ), false );
    double expected_v = t.evaluateWithSympy( "1*1-2" );
    t.check( e00, expected_v );
    
    auto e01 = Px() + expr( "1." );
    BOOST_CHECK_EQUAL( hasDynamicContext( e01 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e01 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e01 ) ), true );
    BOOST_CHECK_EQUAL( hasValue( staticContext( e01 ), vm::POINT ), true );
    expected_v = t.evaluateWithSympy( "x+1" );
    t.check( e01, expected_v );
    
    auto e1 = sin( cst( 1. ) + expr( "z:x:y:z" ) );
    BOOST_CHECK_EQUAL( hasDynamicContext( e1 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e1 ) ), true );
    BOOST_CHECK_EQUAL( hasValue( staticContext( e1 ), vm::POINT ), false );
    expected_v = t.evaluateWithSympy( "sin(1+z)" );
    t.check( e1, expected_v, 1e-2 );
}

BOOST_AUTO_TEST_CASE( t2 )
{
    using namespace Feel::vf;

    TestIntegrate<3> t;
    double expected_v = 0.;
    for ( std::string s : {"nx:nx", "ny:ny", "nz:nz", "nx+ny:nx:ny", "nx+ny+nz:nx:ny:nz", "2*x*nx+y*ny-2*z*nz:x:y:z:nx:ny:nz"} )
    {
        BOOST_TEST_MESSAGE( "testing expression " << s );
        auto e3 = expr( s );
        BOOST_CHECK_EQUAL( hasDynamicContext( e3 ), true );
        BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e3 ) ), false );
        BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e3 ) ), true );
        BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e3 ) ), true );

        if ( s == "2*x*nx+y*ny-2*z*nz:x:y:z:nx:ny:nz" )
        {
            BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e3 ) ), true );
            expected_v = t.evaluateOnBoundary( 2*Px()*Nx()+Py()*Ny()-2*Pz()*Nz() );
            double expected_v2 = t.evaluateWithSympy( "1" );
            BOOST_CHECK_CLOSE( expected_v, expected_v2, 1e-10 );
            t.checkOnBoundary( e3, expected_v, 1e-10 );
        }
        else // integrals are equals to 0
            t.checkSmallOnBoundary( e3, 1e-10 );
    }
}

BOOST_AUTO_TEST_CASE( t3 )
{
    using namespace Feel::vf;
    TestIntegrate<3> t;
    auto mesh = t.mesh;
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    
    u.on(_range=elements(mesh), _expr=expr("x:x"));
    auto e0 = expr( "u:u", "u", idv( u ) );
    
    BOOST_CHECK_EQUAL( hasDynamicContext( e0 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e0 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasGRAD( dynamicContext( e0 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e0 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e0 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e0 ) ), false );
    double expected_v = t.evaluateWithSympy( "x" );
    t.check( e0, expected_v );
    
    
    auto e1 = expr( "u*x:x:u", "u", idv( u ) );
    BOOST_CHECK_EQUAL( hasDynamicContext( e1 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasGRAD( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e1 ) ), true );
    BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e1 ) ), false );
    expected_v = t.evaluateWithSympy( "x*x" );
    t.check( e1, expected_v );
    
    auto e2 = expr( "u*x:x:u", "u", gradv( u )( 0, 0 ) );
    BOOST_CHECK_EQUAL( hasDynamicContext( e2 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e2 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasGRAD( dynamicContext( e2 ) ), true );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e2 ) ), true );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e2 ) ), true );
    BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e2 ) ), false );
    expected_v = t.evaluateWithSympy( "x" );
    t.check( e2, expected_v );
    
    auto e3 = expr( "2*x*u*v:x:u:v" );
    auto e3v = expr( e3, symbolExpr( "u", dxv( u ) ), symbolExpr( "v", Nx() ) );
    BOOST_CHECK_EQUAL( hasDynamicContext( e3v ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e3v ) ), false );
    BOOST_CHECK_EQUAL( vm::hasGRAD( dynamicContext( e3v ) ), true );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e3v ) ), true );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e3v ) ), true );
    BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e3v ) ), true );
    expected_v = t.evaluateOnBoundary( 2*Px()*dxv(u)*Nx() );
    t.checkOnBoundary( e3v, expected_v );

    auto e40 = expr("2*x*nx+3*y*ny+4*z*nz:x:y:z:nx:ny:nz");
    BOOST_CHECK_EQUAL( hasDynamicContext( e40 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e40 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasGRAD( dynamicContext( e40 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e40 ) ), true );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e40 ) ), true );
    BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e40 ) ), true );
    expected_v = t.evaluateWithSympy( "9" );
    t.checkOnBoundary( e40, expected_v );
    
    u.on(_range=elements(mesh), _expr=expr( "x+y+z:x:y:z" ) );
    auto e4 = expr( "u*nx+u*ny+u*nz:u:nx:ny:nz", "u", idv( u ) );
    BOOST_CHECK_EQUAL( hasDynamicContext( e4 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e4 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasGRAD( dynamicContext( e4 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e4 ) ), true );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e4 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e4 ) ), true );
    expected_v = t.evaluateWithSympy( "3" );
    t.checkOnBoundary( e4, expected_v );
}
#endif
BOOST_AUTO_TEST_SUITE_END()

