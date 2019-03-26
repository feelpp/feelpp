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

/** use Feel namespace */
using namespace Feel;

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

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() )
BOOST_AUTO_TEST_SUITE( context_suite )

BOOST_AUTO_TEST_CASE( t0 )
{
    using namespace Feel::vf;

    auto e = Px();
    BOOST_CHECK_EQUAL( hasDynamicContext( e ), false );
    BOOST_CHECK_EQUAL( hasValue( staticContext( e ), vm::POINT ), true );

    auto e1 = expr( "1." );
    BOOST_CHECK_EQUAL( hasDynamicContext( e1 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e1 ) ), false );

    for ( auto s : {"x:x", "y:y", "z:z", "x+y:x:y", "x+y+z:x:y:z"} )
    {
        auto e3 = expr( s );
        BOOST_CHECK_EQUAL( hasDynamicContext( e3 ), true );
        BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e3 ) ), false );
        BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e3 ) ), true );
    }
}

BOOST_AUTO_TEST_CASE( t1 )
{
    using namespace Feel::vf;

    auto e00 = cst( 1.0 ) * expr( "1." ) - expr( "2." );
    BOOST_CHECK_EQUAL( hasDynamicContext( e00 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e00 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e00 ) ), false );
    BOOST_CHECK_EQUAL( hasValue( staticContext( e00 ), vm::POINT ), false );

    auto e01 = Px() + expr( "1." );
    BOOST_CHECK_EQUAL( hasDynamicContext( e01 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e01 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e01 ) ), true );
    BOOST_CHECK_EQUAL( hasValue( staticContext( e01 ), vm::POINT ), true );

    auto e1 = sin( cst( 1. ) + expr( "z:x:y:z" ) );
    BOOST_CHECK_EQUAL( hasDynamicContext( e1 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e1 ) ), true );
    BOOST_CHECK_EQUAL( hasValue( staticContext( e1 ), vm::POINT ), false );
}

BOOST_AUTO_TEST_CASE( t2 )
{
    using namespace Feel::vf;
    for ( auto s : {"nx:nx", "ny:ny", "nz:nz", "nx+ny:nx:ny", "nx+ny+nz:nx:ny:nz"} )
    {
        auto e3 = expr( s );
        BOOST_CHECK_EQUAL( hasDynamicContext( e3 ), true );
        BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e3 ) ), false );
        BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e3 ) ), true );
        BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e3 ) ), true );
    }
}

BOOST_AUTO_TEST_CASE( t3 )
{
    using namespace Feel::vf;
    auto mesh = unitSegment();
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();

    auto e0 = expr( "u:u", "u", idv( u ) );
    BOOST_CHECK_EQUAL( hasDynamicContext( e0 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e0 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasGRAD( dynamicContext( e0 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e0 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e0 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e0 ) ), false );
    
    auto e1 = expr( "u*x:x:u", "u", idv( u ) );
    BOOST_CHECK_EQUAL( hasDynamicContext( e1 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasGRAD( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e1 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e1 ) ), true );
    BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e1 ) ), false );
    
    auto e2 = expr( "u*x:x:u", "u", gradv( u )( 0, 0 ) );
    BOOST_CHECK_EQUAL( hasDynamicContext( e2 ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e2 ) ), false );
    BOOST_CHECK_EQUAL( vm::hasGRAD( dynamicContext( e2 ) ), true );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e2 ) ), true );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e2 ) ), true );
    BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e2 ) ), false );
    
    auto e3 = expr( "2*x*u*v:x:u:v" );
    auto e3v = expr( e3, symbolExpr( "u", dxv( u ) ), symbolExpr( "v", N() ) );
    BOOST_CHECK_EQUAL( hasDynamicContext( e3v ), true );
    BOOST_CHECK_EQUAL( vm::hasJACOBIAN( dynamicContext( e3v ) ), false );
    BOOST_CHECK_EQUAL( vm::hasGRAD( dynamicContext( e3v ) ), true );
    BOOST_CHECK_EQUAL( vm::hasKB( dynamicContext( e3v ) ), true );
    BOOST_CHECK_EQUAL( vm::hasPOINT( dynamicContext( e3v ) ), true );
    BOOST_CHECK_EQUAL( vm::hasNORMAL( dynamicContext( e3v ) ), true );
}
BOOST_AUTO_TEST_SUITE_END()
