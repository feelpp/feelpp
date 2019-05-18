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
//! @date 05 Mar 2019
//! @copyright 2019 Feel++ Consortium
//!

#define BOOST_TEST_MODULE test_expr
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
    AboutData about( "test_expr",
                     "test_expr",
                     "0.1",
                     "test expr",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2019 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( context_suite )

template<int Dim = 3>
class TestIntegrate
{
public:
    using mesh_type=Mesh<Simplex<Dim,1>>;
    TestIntegrate( double meshSize = 0.1 )
        {
            if constexpr( Dim == 3 )
            {
                GeoTool::Node x1( 0, 0, 0);
                GeoTool::Node x2( 1, 0, 0 );
                GeoTool::Node x3( 0, 1, 0 );
                GeoTool::Node x4( 0, 0, 1 );
                GeoTool::Tetrahedron C( meshSize, "OMEGA", x1, x2, x3, x4 );
                mesh = C.createMesh( _mesh = new mesh_type,
                                     _name = "omega_" + mesh_type::shape_type::name(),
                                     _partitions = Environment::worldComm().localSize() );
            }
            else if constexpr ( Dim ==1 )
                mesh = unitSegment();
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
    std::shared_ptr<mesh_type> mesh;
};

    
BOOST_AUTO_TEST_CASE( t0 )
{
    using namespace Feel::vf;
    TestIntegrate<3> t;
    
    auto e = Px();
    BOOST_CHECK_EQUAL( hasDynamicContext( e ), false );
    BOOST_CHECK_EQUAL( hasValue( staticContext( e ), vm::POINT ), true );

    auto e1 = expr( "1." );
    BOOST_CHECK_EQUAL( hasDynamicContext( e1 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e1 ) ), false );
    double v1 = t.evaluateInDomain( cst(1.) );
    t.check( e1, v1 );
    
    for ( std::string s : {"x", "y", "z", "x+y", "x+y+z"} )
    {
        auto e3 = expr( s + ":x:y:z" );
        BOOST_CHECK_EQUAL( hasDynamicContext( e3 ), true );
        BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e3 ) ), false );
        BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e3 ) ), true );
        double expected_v = t.evaluateWithSympy( s );
        t.check( e3, expected_v );
    }
}

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

    //TestIntegrate<3> t;
    for ( auto s : {"nx:nx", "ny:ny", "nz:nz", "nx+ny:nx:ny", "nx+ny+nz:nx:ny:nz"} )
    {
        auto e3 = expr( s );
        BOOST_CHECK_EQUAL( hasDynamicContext( e3 ), true );
        BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e3 ) ), false );
        BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e3 ) ), true );
        BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e3 ) ), true );
        //expected_v = t.evaluateOnBoundary( Nx() );
        //t.checkOnBoundary( e3, expected_v, 1e-2 );
    }
}

BOOST_AUTO_TEST_CASE( t3 )
{
    using namespace Feel::vf;
    TestIntegrate<1> t;
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
}
BOOST_AUTO_TEST_SUITE_END()
