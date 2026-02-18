/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-01-05

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
 * @file test_pch_dynamic.cpp
 * @brief Test dynamic order Pch function spaces
 * @author Christophe Prud'homme
 * @date 2026-01-05
 */
#define BOOST_TEST_MODULE test_pch_dynamic
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelpoly/dubiner.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <sstream>
#include <array>
#include <limits>
#include <cmath>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
using oneelt_mesh_type = Mesh<Simplex<2,1>>;
using oneelt_mesh_ptrtype = std::shared_ptr<oneelt_mesh_type>;

inline oneelt_mesh_ptrtype
makeSingleSimplexElementMesh( std::string const& name, Gmsh::DomainType domainType )
{
    auto desc = std::make_shared<GmshSimplexDomain>( 2, 1, domainType );
    desc->setPrefix( name );
    desc->usePhysicalNames( true );

    if ( domainType == Gmsh::GMSH_REFERENCE_DOMAIN )
    {
        desc->setCharacteristicLength( 2.0 );
    }
    else
    {
        // Keep one affine element, but not the reference simplex.
        desc->setX( { 0.2, 1.4 } );
        desc->setY( { -0.1, 1.1 } );
        desc->setCharacteristicLength( 1.2 );
    }

    auto mesh = createGMSHMesh( _mesh = new oneelt_mesh_type,
                                _desc = desc,
                                _h = desc->h() );
    return mesh;
}

inline std::array<double,4>
singleElementBounds2D( oneelt_mesh_ptrtype const& mesh )
{
    auto const& elt = mesh->beginElement()->second;
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();

    for ( uint16_type i = 0; i < oneelt_mesh_type::element_type::numVertices; ++i )
    {
        auto const& node = elt.point( i ).node();
        xmin = std::min( xmin, node[0] );
        xmax = std::max( xmax, node[0] );
        ymin = std::min( ymin, node[1] );
        ymax = std::max( ymax, node[1] );
    }
    return { xmin, xmax, ymin, ymax };
}

template<int Order, typename ExprT>
void
checkExactManufacturedPolynomial( oneelt_mesh_ptrtype const& mesh,
                                  ExprT const& expr,
                                  std::string const& context,
                                  double tol )
{
    auto VhStatic = Pch<Order>( mesh );
    auto VhDynamic = Pch<Dynamic>( mesh, RuntimeOrder( Order ) );

    BOOST_REQUIRE_EQUAL( VhDynamic->runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( VhStatic->nDof(), VhDynamic->nDof() );

    auto uStatic = project( _space = VhStatic, _expr = expr );
    auto uDynamic = project( _space = VhDynamic, _expr = expr );

    const double errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - expr );
    const double errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - expr );
    const double errStaticVsDynamic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ) );

    BOOST_TEST_MESSAGE( context << " P" << Order
                        << ": nDof(static)=" << VhStatic->nDof()
                        << ", nDof(dynamic)=" << VhDynamic->nDof()
                        << ", errStatic=" << errStatic
                        << ", errDynamic=" << errDynamic
                        << ", errStaticVsDynamic=" << errStaticVsDynamic );

    BOOST_CHECK_SMALL( errStatic, tol );
    BOOST_CHECK_SMALL( errDynamic, tol );
    BOOST_CHECK_SMALL( errStaticVsDynamic, tol );
}

template<int Order>
auto makeManufacturedVectorExpression()
{
    if constexpr ( Order == 1 )
        return vec( cst( 1.0 ) + 2.0 * Px() - Py(),
                    cst( -0.5 ) + Px() + 3.0 * Py() );
    else if constexpr ( Order == 2 )
        return vec( cst( 1.0 ) + Px() + Py() + Px() * Px(),
                    cst( -0.5 ) + 2.0 * Px() * Py() + Py() * Py() );
    else
        return vec( cst( 1.0 ) + Px() * Px() * Px() + Px() * Py(),
                    cst( -0.5 ) + Py() * Py() * Py() + Px() * Py() * Py() );
}

template<int Order>
auto makeManufacturedTensor2Expression()
{
    if constexpr ( Order == 1 )
        return mat<2,2>( cst( 1.0 ) + Px(),
                         Py(),
                         cst( -1.0 ) + 2.0 * Px(),
                         cst( 2.0 ) - Py() );
    else if constexpr ( Order == 2 )
        return mat<2,2>( cst( 1.0 ) + Px() * Px() + Py(),
                         Px() * Py(),
                         Py() * Py() + Px(),
                         cst( 2.0 ) - Px() * Py() );
    else
        return mat<2,2>( cst( 1.0 ) + Px() * Px() * Px(),
                         Px() * Px() * Py(),
                         Px() * Py() * Py(),
                         cst( 2.0 ) + Py() * Py() * Py() );
}

template<int Order>
void checkPdhvStaticDynamic( oneelt_mesh_ptrtype const& mesh, double tol )
{
    auto VhStatic = Pdhv<Order>( mesh );
    auto VhDynamic = Pdhv<Dynamic>( mesh, RuntimeOrder( Order ) );
    auto expr = makeManufacturedVectorExpression<Order>();

    BOOST_REQUIRE_EQUAL( VhDynamic->runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( VhStatic->nDof(), VhDynamic->nDof() );

    auto uStatic = VhStatic->element();
    auto uDynamic = VhDynamic->element();
    uStatic.on( _range = elements( mesh ), _expr = expr );
    uDynamic.on( _range = elements( mesh ), _expr = expr );

    const double errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - expr, _quad = _Q<12>() );
    const double errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - expr, _quad = _Q<12>() );
    const double errStaticVsDynamic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<12>() );

    BOOST_TEST_MESSAGE( "Pdhv static vs dynamic P" << Order
                        << ": nDof(static)=" << VhStatic->nDof()
                        << ", nDof(dynamic)=" << VhDynamic->nDof()
                        << ", errStatic=" << errStatic
                        << ", errDynamic=" << errDynamic
                        << ", errStaticVsDynamic=" << errStaticVsDynamic );

    BOOST_CHECK_SMALL( errStatic, tol );
    BOOST_CHECK_SMALL( errDynamic, tol );
    BOOST_CHECK_SMALL( errStaticVsDynamic, tol );
}

template<int Order>
void checkPdhmStaticDynamic( oneelt_mesh_ptrtype const& mesh, double tol )
{
    auto MhStatic = Pdhm<Order>( mesh );
    auto MhDynamic = Pdhm<Dynamic>( mesh, RuntimeOrder( Order ) );
    auto expr = makeManufacturedTensor2Expression<Order>();

    BOOST_REQUIRE_EQUAL( MhDynamic->runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( MhStatic->nDof(), MhDynamic->nDof() );

    auto mStatic = MhStatic->element();
    auto mDynamic = MhDynamic->element();
    mStatic.on( _range = elements( mesh ), _expr = expr );
    mDynamic.on( _range = elements( mesh ), _expr = expr );

    const double errStatic = normL2( _range = elements( mesh ), _expr = idv( mStatic ) - expr, _quad = _Q<12>() );
    const double errDynamic = normL2( _range = elements( mesh ), _expr = idv( mDynamic ) - expr, _quad = _Q<12>() );
    const double errStaticVsDynamic = normL2( _range = elements( mesh ), _expr = idv( mStatic ) - idv( mDynamic ), _quad = _Q<12>() );

    BOOST_TEST_MESSAGE( "Pdhm static vs dynamic P" << Order
                        << ": nDof(static)=" << MhStatic->nDof()
                        << ", nDof(dynamic)=" << MhDynamic->nDof()
                        << ", errStatic=" << errStatic
                        << ", errDynamic=" << errDynamic
                        << ", errStaticVsDynamic=" << errStaticVsDynamic );

    BOOST_CHECK_SMALL( errStatic, tol );
    BOOST_CHECK_SMALL( errDynamic, tol );
    BOOST_CHECK_SMALL( errStaticVsDynamic, tol );
}

template<int Order>
void checkPchmStaticDynamic( oneelt_mesh_ptrtype const& mesh, double tol )
{
    auto MhStatic = Pchm<Order>( mesh );
    auto MhDynamic = Pchm<Dynamic>( mesh, RuntimeOrder( Order ) );
    auto expr = makeManufacturedTensor2Expression<Order>();

    BOOST_REQUIRE_EQUAL( MhDynamic->runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( MhStatic->nDof(), MhDynamic->nDof() );

    auto mStatic = MhStatic->element();
    auto mDynamic = MhDynamic->element();
    mStatic.on( _range = elements( mesh ), _expr = expr );
    mDynamic.on( _range = elements( mesh ), _expr = expr );

    const double errStatic = normL2( _range = elements( mesh ), _expr = idv( mStatic ) - expr, _quad = _Q<12>() );
    const double errDynamic = normL2( _range = elements( mesh ), _expr = idv( mDynamic ) - expr, _quad = _Q<12>() );
    const double errStaticVsDynamic = normL2( _range = elements( mesh ), _expr = idv( mStatic ) - idv( mDynamic ), _quad = _Q<12>() );

    BOOST_TEST_MESSAGE( "Pchm static vs dynamic P" << Order
                        << ": nDof(static)=" << MhStatic->nDof()
                        << ", nDof(dynamic)=" << MhDynamic->nDof()
                        << ", errStatic=" << errStatic
                        << ", errDynamic=" << errDynamic
                        << ", errStaticVsDynamic=" << errStaticVsDynamic );

    BOOST_CHECK_SMALL( errStatic, tol );
    BOOST_CHECK_SMALL( errDynamic, tol );
    BOOST_CHECK_SMALL( errStaticVsDynamic, tol );
}

template<int Order>
std::string
makeManufacturedPolynomialString()
{
    static_assert( Order >= 1, "Order must be >= 1" );
    std::ostringstream os;
    // Build a total-degree Order polynomial with mixed x/y monomials.
    // Example (Order=4): 1*x^4 + 2*x^3*y + ... + 5*y^4
    for ( int i = 0; i <= Order; ++i )
    {
        const int px = Order - i;
        const int py = i;
        if ( i > 0 )
            os << " + ";
        os << ( i + 1 ) << "*";

        bool hasVar = false;
        if ( px > 0 )
        {
            os << "x";
            if ( px > 1 )
                os << "^" << px;
            hasVar = true;
        }
        if ( py > 0 )
        {
            if ( hasVar )
                os << "*";
            os << "y";
            if ( py > 1 )
                os << "^" << py;
            hasVar = true;
        }
        if ( !hasVar )
            os << "1";
    }
    os << ":x:y";
    return os.str();
}

template<typename SpacePtrType, typename GExprT, typename LapExprT>
auto solveDirichletLaplacian( SpacePtrType const& Vh, GExprT const& gexpr, LapExprT const& lapexpr )
{
    auto u = Vh->element();
    auto v = Vh->element();
    auto b = backend( _name = "", _kind = soption( _name = "backend" ), _rebuild = true );

    auto l = form1( _test = Vh, _backend = b );
    l = integrate( _range = elements( Vh->mesh() ),
                   _expr = -lapexpr * id( v ),
                   _quad = _Q<8>() );

    auto a = form2( _trial = Vh, _test = Vh, _backend = b );
    a = integrate( _range = elements( Vh->mesh() ),
                   _expr = gradt( u ) * trans( grad( v ) ),
                   _quad = _Q<8>() );
    a += on( _range = boundaryfaces( Vh->mesh() ), _rhs = l, _element = u, _expr = gexpr, _type = "elimination" );
    a.solveb( _rhs = l, _solution = u, _backend = b );

    return u;
}

template<typename SpacePtrType, typename GExprT>
auto solveMassProjection( SpacePtrType const& Vh, GExprT const& gexpr )
{
    auto u = Vh->element();
    auto v = Vh->element();
    auto b = backend( _name = "", _kind = soption( _name = "backend" ), _rebuild = true );

    auto l = form1( _test = Vh, _backend = b );
    l = integrate( _range = elements( Vh->mesh() ),
                   _expr = gexpr * id( v ),
                   _quad = _Q<8>() );

    auto a = form2( _trial = Vh, _test = Vh, _backend = b );
    a = integrate( _range = elements( Vh->mesh() ),
                   _expr = idt( u ) * id( v ),
                   _quad = _Q<8>() );
    a.solveb( _rhs = l, _solution = u, _backend = b );

    return u;
}

template<int Order>
void
checkSimpleLapEndToEndStaticDynamic( oneelt_mesh_ptrtype const& mesh, double tol )
{
    auto VhStatic = Pch<Order>( mesh );
    auto VhDynamic = Pch<Dynamic>( mesh, RuntimeOrder( Order ) );

    BOOST_REQUIRE_EQUAL( VhDynamic->runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( VhStatic->nDof(), VhDynamic->nDof() );

    if constexpr ( Order == 1 )
    {
        auto gexpr = cst( 1.0 ) + 2.0 * Px() - 3.0 * Py();
        auto ggrad = trans( vec( cst( 2.0 ), cst( -3.0 ) ) );
        auto lapexpr = cst( 0.0 );
        auto uStatic = solveDirichletLaplacian( VhStatic, gexpr, lapexpr );
        auto uDynamic = solveDirichletLaplacian( VhDynamic, gexpr, lapexpr );

        const double errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - gexpr, _quad = _Q<8>() );
        const double errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - gexpr, _quad = _Q<8>() );
        const double errStaticVsDynamic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<8>() );
        const double errGradStatic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - ggrad, _quad = _Q<8>() );
        const double errGradDynamic = normL2( _range = elements( mesh ), _expr = gradv( uDynamic ) - ggrad, _quad = _Q<8>() );
        const double errGradStaticVsDynamic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - gradv( uDynamic ), _quad = _Q<8>() );

        BOOST_TEST_MESSAGE( "end-to-end P1 (laplacian): errStatic=" << errStatic
                            << ", errDynamic=" << errDynamic
                            << ", errStaticVsDynamic=" << errStaticVsDynamic
                            << ", gradErrStatic=" << errGradStatic
                            << ", gradErrDynamic=" << errGradDynamic
                            << ", gradErrStaticVsDynamic=" << errGradStaticVsDynamic );

        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errStaticVsDynamic, tol );
        BOOST_CHECK_SMALL( errGradStatic, tol );
        BOOST_CHECK_SMALL( errGradDynamic, tol );
        BOOST_CHECK_SMALL( errGradStaticVsDynamic, tol );
    }
    else if constexpr ( Order == 2 )
    {
        auto gexpr = cst( 1.0 ) + Px() + 2.0 * Py() + Px() * Px() + Px() * Py() + 2.0 * Py() * Py();
        auto ggrad = trans( vec( cst( 1.0 ) + 2.0 * Px() + Py(),
                                 cst( 2.0 ) + Px() + 4.0 * Py() ) );
        auto lapexpr = cst( 6.0 );
        auto uStatic = solveDirichletLaplacian( VhStatic, gexpr, lapexpr );
        auto uDynamic = solveDirichletLaplacian( VhDynamic, gexpr, lapexpr );

        const double errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - gexpr, _quad = _Q<8>() );
        const double errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - gexpr, _quad = _Q<8>() );
        const double errStaticVsDynamic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<8>() );
        const double errGradStatic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - ggrad, _quad = _Q<8>() );
        const double errGradDynamic = normL2( _range = elements( mesh ), _expr = gradv( uDynamic ) - ggrad, _quad = _Q<8>() );
        const double errGradStaticVsDynamic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - gradv( uDynamic ), _quad = _Q<8>() );

        BOOST_TEST_MESSAGE( "end-to-end P2 (laplacian): errStatic=" << errStatic
                            << ", errDynamic=" << errDynamic
                            << ", errStaticVsDynamic=" << errStaticVsDynamic
                            << ", gradErrStatic=" << errGradStatic
                            << ", gradErrDynamic=" << errGradDynamic
                            << ", gradErrStaticVsDynamic=" << errGradStaticVsDynamic );

        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errStaticVsDynamic, tol );
        BOOST_CHECK_SMALL( errGradStatic, tol );
        BOOST_CHECK_SMALL( errGradDynamic, tol );
        BOOST_CHECK_SMALL( errGradStaticVsDynamic, tol );
    }
    else if constexpr ( Order == 3 )
    {
        auto gexpr = cst( 1.0 ) + Px() + Py() + Px() * Px() + Px() * Py() + Py() * Py()
                     + Px() * Px() * Px() + Px() * Px() * Py() + Px() * Py() * Py() + Py() * Py() * Py();
        auto ggrad = trans( vec( cst( 1.0 ) + 2.0 * Px() + Py() + 3.0 * Px() * Px() + 2.0 * Px() * Py() + Py() * Py(),
                                 cst( 1.0 ) + Px() + 2.0 * Py() + Px() * Px() + 2.0 * Px() * Py() + 3.0 * Py() * Py() ) );
        auto lapexpr = cst( 4.0 ) + 8.0 * Px() + 8.0 * Py();
        auto uStatic = solveDirichletLaplacian( VhStatic, gexpr, lapexpr );
        auto uDynamic = solveDirichletLaplacian( VhDynamic, gexpr, lapexpr );

        const double errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - gexpr, _quad = _Q<8>() );
        const double errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - gexpr, _quad = _Q<8>() );
        const double errStaticVsDynamic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<8>() );
        const double errGradStatic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - ggrad, _quad = _Q<8>() );
        const double errGradDynamic = normL2( _range = elements( mesh ), _expr = gradv( uDynamic ) - ggrad, _quad = _Q<8>() );
        const double errGradStaticVsDynamic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - gradv( uDynamic ), _quad = _Q<8>() );

        BOOST_TEST_MESSAGE( "end-to-end P3 (laplacian): errStatic=" << errStatic
                            << ", errDynamic=" << errDynamic
                            << ", errStaticVsDynamic=" << errStaticVsDynamic
                            << ", gradErrStatic=" << errGradStatic
                            << ", gradErrDynamic=" << errGradDynamic
                            << ", gradErrStaticVsDynamic=" << errGradStaticVsDynamic );

        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errStaticVsDynamic, tol );
        BOOST_CHECK_SMALL( errGradStatic, tol );
        BOOST_CHECK_SMALL( errGradDynamic, tol );
        BOOST_CHECK_SMALL( errGradStaticVsDynamic, tol );
    }
}

template<int Order>
void
checkAnalyticalGradientIntegrals( oneelt_mesh_ptrtype const& mesh, double tol )
{
    auto VhStatic = Pch<Order>( mesh );
    auto VhDynamic = Pch<Dynamic>( mesh, RuntimeOrder( Order ) );

    BOOST_REQUIRE_EQUAL( VhDynamic->runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( VhStatic->nDof(), VhDynamic->nDof() );

    auto runCase = [&]( auto const& polyExpr,
                        auto const& gradExpr,
                        double expectedGradX,
                        double expectedGradY,
                        std::string const& label )
    {
        auto uStatic = project( _space = VhStatic, _expr = polyExpr );
        auto uDynamic = project( _space = VhDynamic, _expr = polyExpr );

        const double intGradXStatic = integrate( _range = elements( mesh ), _expr = gradv( uStatic )( 0, 0 ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double intGradYStatic = integrate( _range = elements( mesh ), _expr = gradv( uStatic )( 0, 1 ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double intGradXDynamic = integrate( _range = elements( mesh ), _expr = gradv( uDynamic )( 0, 0 ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double intGradYDynamic = integrate( _range = elements( mesh ), _expr = gradv( uDynamic )( 0, 1 ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - polyExpr, _quad = _Q<10>() );
        const double errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - polyExpr, _quad = _Q<10>() );
        const double errStaticVsDynamic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<10>() );
        const double errGradStatic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - gradExpr, _quad = _Q<10>() );
        const double errGradDynamic = normL2( _range = elements( mesh ), _expr = gradv( uDynamic ) - gradExpr, _quad = _Q<10>() );
        const double errGradStaticVsDynamic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - gradv( uDynamic ), _quad = _Q<10>() );

        BOOST_TEST_MESSAGE( "analytical grad integral " << label
                            << " P" << Order
                            << ": static=(" << intGradXStatic << "," << intGradYStatic << ")"
                            << ", dynamic=(" << intGradXDynamic << "," << intGradYDynamic << ")"
                            << ", expected=(" << expectedGradX << "," << expectedGradY << ")"
                            << ", errStatic=" << errStatic
                            << ", errDynamic=" << errDynamic
                            << ", errStaticVsDynamic=" << errStaticVsDynamic
                            << ", errGradStatic=" << errGradStatic
                            << ", errGradDynamic=" << errGradDynamic
                            << ", errGradStaticVsDynamic=" << errGradStaticVsDynamic );

        BOOST_CHECK_SMALL( std::abs( intGradXStatic - expectedGradX ), tol );
        BOOST_CHECK_SMALL( std::abs( intGradYStatic - expectedGradY ), tol );
        BOOST_CHECK_SMALL( std::abs( intGradXDynamic - expectedGradX ), tol );
        BOOST_CHECK_SMALL( std::abs( intGradYDynamic - expectedGradY ), tol );
        BOOST_CHECK_SMALL( std::abs( intGradXStatic - intGradXDynamic ), tol );
        BOOST_CHECK_SMALL( std::abs( intGradYStatic - intGradYDynamic ), tol );
        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errStaticVsDynamic, tol );
        BOOST_CHECK_SMALL( errGradStatic, tol );
        BOOST_CHECK_SMALL( errGradDynamic, tol );
        BOOST_CHECK_SMALL( errGradStaticVsDynamic, tol );
    };

    if constexpr ( Order == 1 )
    {
        // On reference simplex with vertices (-1,-1), (1,-1), (-1,1), area = 2:
        // grad(x+y) = (1,1), hence integral = (2,2).
        runCase( Px() + Py(),
                 trans( vec( cst( 1.0 ), cst( 1.0 ) ) ),
                 2.0,
                 2.0,
                 "u=x+y" );
    }
    else if constexpr ( Order == 2 )
    {
        // grad(x^2+y^2) = (2x,2y), and over the same simplex:
        // integral 2x = integral 2y = -4/3.
        runCase( Px() * Px() + Py() * Py(),
                 trans( vec( 2.0 * Px(), 2.0 * Py() ) ),
                 -4.0 / 3.0,
                 -4.0 / 3.0,
                 "u=x^2+y^2" );
    }
    else if constexpr ( Order == 3 )
    {
        // grad(x^3+y^3) = (3x^2,3y^2), and over the same simplex:
        // integral 3x^2 = integral 3y^2 = 2.
        runCase( Px() * Px() * Px() + Py() * Py() * Py(),
                 trans( vec( 3.0 * Px() * Px(), 3.0 * Py() * Py() ) ),
                 2.0,
                 2.0,
                 "u=x^3+y^3" );
    }
}

template<int Order>
void
checkAnalyticalSemiH1( oneelt_mesh_ptrtype const& mesh, double tol )
{
    auto VhStatic = Pch<Order>( mesh );
    auto VhDynamic = Pch<Dynamic>( mesh, RuntimeOrder( Order ) );

    BOOST_REQUIRE_EQUAL( VhDynamic->runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( VhStatic->nDof(), VhDynamic->nDof() );

    auto runCase = [&]( auto const& polyExpr,
                        auto const& gradExpr,
                        double expectedSemiH1Sq,
                        std::string const& label )
    {
        auto uStatic = project( _space = VhStatic, _expr = polyExpr );
        auto uDynamic = project( _space = VhDynamic, _expr = polyExpr );

        const double semiH1SqStatic = integrate( _range = elements( mesh ), _expr = inner( gradv( uStatic ) ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double semiH1SqDynamic = integrate( _range = elements( mesh ), _expr = inner( gradv( uDynamic ) ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - polyExpr, _quad = _Q<10>() );
        const double errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - polyExpr, _quad = _Q<10>() );
        const double errStaticVsDynamic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<10>() );
        const double errGradStatic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - gradExpr, _quad = _Q<10>() );
        const double errGradDynamic = normL2( _range = elements( mesh ), _expr = gradv( uDynamic ) - gradExpr, _quad = _Q<10>() );
        const double errGradStaticVsDynamic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - gradv( uDynamic ), _quad = _Q<10>() );

        const double semiH1Static = normSemiH1( _range = elements( mesh ), _grad_expr = gradv( uStatic ), _quad = _Q<10>() );
        const double semiH1Dynamic = normSemiH1( _range = elements( mesh ), _grad_expr = gradv( uDynamic ), _quad = _Q<10>() );
        const double expectedSemiH1 = std::sqrt( expectedSemiH1Sq );

        BOOST_TEST_MESSAGE( "analytical semi-H1 " << label
                            << " P" << Order
                            << ": static=" << semiH1Static
                            << ", dynamic=" << semiH1Dynamic
                            << ", expected=" << expectedSemiH1
                            << ", errStatic=" << errStatic
                            << ", errDynamic=" << errDynamic
                            << ", errStaticVsDynamic=" << errStaticVsDynamic
                            << ", errGradStatic=" << errGradStatic
                            << ", errGradDynamic=" << errGradDynamic
                            << ", errGradStaticVsDynamic=" << errGradStaticVsDynamic );

        BOOST_CHECK_SMALL( std::abs( semiH1SqStatic - expectedSemiH1Sq ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1SqDynamic - expectedSemiH1Sq ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1Static - expectedSemiH1 ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1Dynamic - expectedSemiH1 ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1Static - semiH1Dynamic ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1Static - std::sqrt( semiH1SqStatic ) ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1Dynamic - std::sqrt( semiH1SqDynamic ) ), tol );
        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errStaticVsDynamic, tol );
        BOOST_CHECK_SMALL( errGradStatic, tol );
        BOOST_CHECK_SMALL( errGradDynamic, tol );
        BOOST_CHECK_SMALL( errGradStaticVsDynamic, tol );
    };

    if constexpr ( Order == 1 )
    {
        // |grad(x+y)|^2 = 2, integral over area 2 is 4.
        runCase( Px() + Py(),
                 trans( vec( cst( 1.0 ), cst( 1.0 ) ) ),
                 4.0,
                 "u=x+y" );
    }
    else if constexpr ( Order == 2 )
    {
        // |grad(x^2+y^2)|^2 = 4(x^2+y^2), integral is 16/3.
        runCase( Px() * Px() + Py() * Py(),
                 trans( vec( 2.0 * Px(), 2.0 * Py() ) ),
                 16.0 / 3.0,
                 "u=x^2+y^2" );
    }
    else if constexpr ( Order == 3 )
    {
        // |grad(x^3+y^3)|^2 = 9(x^4+y^4), integral is 36/5.
        runCase( Px() * Px() * Px() + Py() * Py() * Py(),
                 trans( vec( 3.0 * Px() * Px(), 3.0 * Py() * Py() ) ),
                 36.0 / 5.0,
                 "u=x^3+y^3" );
    }
}

template<int Order, typename MeshPtrType>
void
checkAnalyticalGradientIntegralsUnitSquare( MeshPtrType const& mesh, double tol )
{
    auto VhStatic = Pch<Order>( mesh );
    auto VhDynamic = Pch<Dynamic>( mesh, RuntimeOrder( Order ) );

    BOOST_REQUIRE_EQUAL( VhDynamic->runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( VhStatic->nDof(), VhDynamic->nDof() );

    auto runCase = [&]( auto const& polyExpr,
                        auto const& gradExpr,
                        double expectedGradX,
                        double expectedGradY,
                        std::string const& label )
    {
        auto uStatic = project( _space = VhStatic, _expr = polyExpr );
        auto uDynamic = project( _space = VhDynamic, _expr = polyExpr );

        const double intGradXStatic = integrate( _range = elements( mesh ), _expr = gradv( uStatic )( 0, 0 ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double intGradYStatic = integrate( _range = elements( mesh ), _expr = gradv( uStatic )( 0, 1 ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double intGradXDynamic = integrate( _range = elements( mesh ), _expr = gradv( uDynamic )( 0, 0 ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double intGradYDynamic = integrate( _range = elements( mesh ), _expr = gradv( uDynamic )( 0, 1 ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - polyExpr, _quad = _Q<10>() );
        const double errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - polyExpr, _quad = _Q<10>() );
        const double errStaticVsDynamic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<10>() );
        const double errGradStatic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - gradExpr, _quad = _Q<10>() );
        const double errGradDynamic = normL2( _range = elements( mesh ), _expr = gradv( uDynamic ) - gradExpr, _quad = _Q<10>() );
        const double errGradStaticVsDynamic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - gradv( uDynamic ), _quad = _Q<10>() );

        BOOST_TEST_MESSAGE( "analytical grad integral (unitSquare multi-elt) " << label
                            << " P" << Order
                            << ": static=(" << intGradXStatic << "," << intGradYStatic << ")"
                            << ", dynamic=(" << intGradXDynamic << "," << intGradYDynamic << ")"
                            << ", expected=(" << expectedGradX << "," << expectedGradY << ")"
                            << ", errStatic=" << errStatic
                            << ", errDynamic=" << errDynamic
                            << ", errStaticVsDynamic=" << errStaticVsDynamic
                            << ", errGradStatic=" << errGradStatic
                            << ", errGradDynamic=" << errGradDynamic
                            << ", errGradStaticVsDynamic=" << errGradStaticVsDynamic );

        BOOST_CHECK_SMALL( std::abs( intGradXStatic - expectedGradX ), tol );
        BOOST_CHECK_SMALL( std::abs( intGradYStatic - expectedGradY ), tol );
        BOOST_CHECK_SMALL( std::abs( intGradXDynamic - expectedGradX ), tol );
        BOOST_CHECK_SMALL( std::abs( intGradYDynamic - expectedGradY ), tol );
        BOOST_CHECK_SMALL( std::abs( intGradXStatic - intGradXDynamic ), tol );
        BOOST_CHECK_SMALL( std::abs( intGradYStatic - intGradYDynamic ), tol );
        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errStaticVsDynamic, tol );
        BOOST_CHECK_SMALL( errGradStatic, tol );
        BOOST_CHECK_SMALL( errGradDynamic, tol );
        BOOST_CHECK_SMALL( errGradStaticVsDynamic, tol );
    };

    if constexpr ( Order == 1 )
    {
        // On unit square [0,1]^2: grad(x+y)=(1,1), so integrals are (1,1).
        runCase( Px() + Py(),
                 trans( vec( cst( 1.0 ), cst( 1.0 ) ) ),
                 1.0,
                 1.0,
                 "u=x+y" );
    }
    else if constexpr ( Order == 2 )
    {
        // grad(x^2+y^2)=(2x,2y), so integrals are (1,1).
        runCase( Px() * Px() + Py() * Py(),
                 trans( vec( 2.0 * Px(), 2.0 * Py() ) ),
                 1.0,
                 1.0,
                 "u=x^2+y^2" );
    }
    else if constexpr ( Order == 3 )
    {
        // grad(x^3+y^3)=(3x^2,3y^2), so integrals are (1,1).
        runCase( Px() * Px() * Px() + Py() * Py() * Py(),
                 trans( vec( 3.0 * Px() * Px(), 3.0 * Py() * Py() ) ),
                 1.0,
                 1.0,
                 "u=x^3+y^3" );
    }
}

template<int Order, typename MeshPtrType>
void
checkAnalyticalSemiH1UnitSquare( MeshPtrType const& mesh, double tol )
{
    auto VhStatic = Pch<Order>( mesh );
    auto VhDynamic = Pch<Dynamic>( mesh, RuntimeOrder( Order ) );

    BOOST_REQUIRE_EQUAL( VhDynamic->runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( VhStatic->nDof(), VhDynamic->nDof() );

    auto runCase = [&]( auto const& polyExpr,
                        auto const& gradExpr,
                        double expectedSemiH1Sq,
                        std::string const& label )
    {
        auto uStatic = project( _space = VhStatic, _expr = polyExpr );
        auto uDynamic = project( _space = VhDynamic, _expr = polyExpr );

        const double semiH1SqStatic = integrate( _range = elements( mesh ), _expr = inner( gradv( uStatic ) ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double semiH1SqDynamic = integrate( _range = elements( mesh ), _expr = inner( gradv( uDynamic ) ), _quad = _Q<10>() ).evaluate()( 0, 0 );
        const double errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - polyExpr, _quad = _Q<10>() );
        const double errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - polyExpr, _quad = _Q<10>() );
        const double errStaticVsDynamic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<10>() );
        const double errGradStatic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - gradExpr, _quad = _Q<10>() );
        const double errGradDynamic = normL2( _range = elements( mesh ), _expr = gradv( uDynamic ) - gradExpr, _quad = _Q<10>() );
        const double errGradStaticVsDynamic = normL2( _range = elements( mesh ), _expr = gradv( uStatic ) - gradv( uDynamic ), _quad = _Q<10>() );

        const double semiH1Static = normSemiH1( _range = elements( mesh ), _grad_expr = gradv( uStatic ), _quad = _Q<10>() );
        const double semiH1Dynamic = normSemiH1( _range = elements( mesh ), _grad_expr = gradv( uDynamic ), _quad = _Q<10>() );
        const double expectedSemiH1 = std::sqrt( expectedSemiH1Sq );

        BOOST_TEST_MESSAGE( "analytical semi-H1 (unitSquare multi-elt) " << label
                            << " P" << Order
                            << ": static=" << semiH1Static
                            << ", dynamic=" << semiH1Dynamic
                            << ", expected=" << expectedSemiH1
                            << ", errStatic=" << errStatic
                            << ", errDynamic=" << errDynamic
                            << ", errStaticVsDynamic=" << errStaticVsDynamic
                            << ", errGradStatic=" << errGradStatic
                            << ", errGradDynamic=" << errGradDynamic
                            << ", errGradStaticVsDynamic=" << errGradStaticVsDynamic );

        BOOST_CHECK_SMALL( std::abs( semiH1SqStatic - expectedSemiH1Sq ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1SqDynamic - expectedSemiH1Sq ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1Static - expectedSemiH1 ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1Dynamic - expectedSemiH1 ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1Static - semiH1Dynamic ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1Static - std::sqrt( semiH1SqStatic ) ), tol );
        BOOST_CHECK_SMALL( std::abs( semiH1Dynamic - std::sqrt( semiH1SqDynamic ) ), tol );
        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errStaticVsDynamic, tol );
        BOOST_CHECK_SMALL( errGradStatic, tol );
        BOOST_CHECK_SMALL( errGradDynamic, tol );
        BOOST_CHECK_SMALL( errGradStaticVsDynamic, tol );
    };

    if constexpr ( Order == 1 )
    {
        // |grad(x+y)|^2 = 2 on unit square, integral is 2.
        runCase( Px() + Py(),
                 trans( vec( cst( 1.0 ), cst( 1.0 ) ) ),
                 2.0,
                 "u=x+y" );
    }
    else if constexpr ( Order == 2 )
    {
        // |grad(x^2+y^2)|^2 = 4(x^2+y^2), integral is 8/3.
        runCase( Px() * Px() + Py() * Py(),
                 trans( vec( 2.0 * Px(), 2.0 * Py() ) ),
                 8.0 / 3.0,
                 "u=x^2+y^2" );
    }
    else if constexpr ( Order == 3 )
    {
        // |grad(x^3+y^3)|^2 = 9(x^4+y^4), integral is 18/5.
        runCase( Px() * Px() * Px() + Py() * Py() * Py(),
                 trans( vec( 3.0 * Px() * Px(), 3.0 * Py() * Py() ) ),
                 18.0 / 5.0,
                 "u=x^3+y^3" );
    }
}
} // namespace

BOOST_AUTO_TEST_SUITE( test_pch_dynamic )

BOOST_AUTO_TEST_CASE( test_pch_dynamic_end_to_end_simple_lap )
{
    BOOST_TEST_MESSAGE( "End-to-end solve (simple_lap style): P1/P2/P3 Laplacian, static vs dynamic" );

    auto mesh = makeSingleSimplexElementMesh( "pch-dynamic-simple-lap", Gmsh::GMSH_REFERENCE_DOMAIN );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_EQUAL( mesh->numElements(), 1 );

    constexpr double tol = 1e-10;
    checkSimpleLapEndToEndStaticDynamic<1>( mesh, tol );
    checkSimpleLapEndToEndStaticDynamic<2>( mesh, tol );
    checkSimpleLapEndToEndStaticDynamic<3>( mesh, tol );
}

BOOST_AUTO_TEST_CASE( test_pch_dynamic_end_to_end_simple_lap_unitsquare_l2 )
{
    BOOST_TEST_MESSAGE( "End-to-end solve on unitSquare (simple_lap style): P1/P2/P3 Laplacian, static vs dynamic (L2 consistency)" );

    auto mesh = unitSquare();
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE( mesh->numElements() > 0 );

    constexpr double tol = 1e-10;
    checkSimpleLapEndToEndStaticDynamic<1>( mesh, tol );
    checkSimpleLapEndToEndStaticDynamic<2>( mesh, tol );
    checkSimpleLapEndToEndStaticDynamic<3>( mesh, tol );
}

BOOST_AUTO_TEST_CASE( test_pch_dynamic_analytical_grad_integrals_reference_simplex )
{
    BOOST_TEST_MESSAGE( "Analytical gradient integral checks on one reference simplex: static vs dynamic P1/P2/P3" );

    auto mesh = makeSingleSimplexElementMesh( "pch-dynamic-analytical-grad", Gmsh::GMSH_REFERENCE_DOMAIN );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_EQUAL( mesh->numElements(), 1 );

    auto const [xmin, xmax, ymin, ymax] = singleElementBounds2D( mesh );
    BOOST_CHECK_SMALL( std::abs( xmin + 1.0 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( xmax - 1.0 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( ymin + 1.0 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( ymax - 1.0 ), 1e-12 );

    constexpr double tol = 1e-10;
    checkAnalyticalGradientIntegrals<1>( mesh, tol );
    checkAnalyticalGradientIntegrals<2>( mesh, tol );
    checkAnalyticalGradientIntegrals<3>( mesh, tol );
}

BOOST_AUTO_TEST_CASE( test_pch_dynamic_analytical_semih1_reference_simplex )
{
    BOOST_TEST_MESSAGE( "Analytical H1 seminorm checks on one reference simplex: static vs dynamic P1/P2/P3" );

    auto mesh = makeSingleSimplexElementMesh( "pch-dynamic-analytical-semih1", Gmsh::GMSH_REFERENCE_DOMAIN );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_EQUAL( mesh->numElements(), 1 );

    auto const [xmin, xmax, ymin, ymax] = singleElementBounds2D( mesh );
    BOOST_CHECK_SMALL( std::abs( xmin + 1.0 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( xmax - 1.0 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( ymin + 1.0 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( ymax - 1.0 ), 1e-12 );

    constexpr double tol = 1e-10;
    checkAnalyticalSemiH1<1>( mesh, tol );
    checkAnalyticalSemiH1<2>( mesh, tol );
    checkAnalyticalSemiH1<3>( mesh, tol );
}

BOOST_AUTO_TEST_CASE( test_pch_dynamic_analytical_grad_integrals_unitsquare_multielements )
{
    BOOST_TEST_MESSAGE( "Analytical gradient integral checks on multi-element unitSquare: static vs dynamic P1/P2/P3" );

    auto mesh = unitSquare( 0.2 );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE( mesh->numElements() > 1 );

    constexpr double tol = 1e-10;
    checkAnalyticalGradientIntegralsUnitSquare<1>( mesh, tol );
    checkAnalyticalGradientIntegralsUnitSquare<2>( mesh, tol );
    checkAnalyticalGradientIntegralsUnitSquare<3>( mesh, tol );
}

BOOST_AUTO_TEST_CASE( test_pch_dynamic_analytical_semih1_unitsquare_multielements )
{
    BOOST_TEST_MESSAGE( "Analytical H1 seminorm checks on multi-element unitSquare: static vs dynamic P1/P2/P3" );

    auto mesh = unitSquare( 0.2 );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE( mesh->numElements() > 1 );

    constexpr double tol = 1e-10;
    checkAnalyticalSemiH1UnitSquare<1>( mesh, tol );
    checkAnalyticalSemiH1UnitSquare<2>( mesh, tol );
    checkAnalyticalSemiH1UnitSquare<3>( mesh, tol );
}

BOOST_AUTO_TEST_CASE( test_pch_dynamic_one_element_reference_identity_geomap )
{
    BOOST_TEST_MESSAGE( "Single-element reference simplex: exact manufactured polynomial checks for P1/P2/P3" );

    auto mesh = makeSingleSimplexElementMesh( "pch-dynamic-oneelt-ref", Gmsh::GMSH_REFERENCE_DOMAIN );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_EQUAL( mesh->numElements(), 1 );

    auto const [xmin, xmax, ymin, ymax] = singleElementBounds2D( mesh );
    BOOST_CHECK_SMALL( std::abs( xmin + 1.0 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( xmax - 1.0 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( ymin + 1.0 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( ymax - 1.0 ), 1e-12 );

    constexpr double tol = 1e-10;

    checkExactManufacturedPolynomial<1>(
        mesh,
        cst( 1.0 ) + 2.0 * Px() - 3.0 * Py(),
        "reference-identity",
        tol );
    checkExactManufacturedPolynomial<2>(
        mesh,
        cst( 1.0 ) + Px() + 2.0 * Py()
        + Px() * Px() + Px() * Py() + 2.0 * Py() * Py(),
        "reference-identity",
        tol );
    checkExactManufacturedPolynomial<3>(
        mesh,
        cst( 1.0 ) + Px() + Py()
        + Px() * Px() + Px() * Py() + Py() * Py()
        + Px() * Px() * Px()
        + Px() * Px() * Py()
        + Px() * Py() * Py()
        + Py() * Py() * Py(),
        "reference-identity",
        tol );
}

BOOST_AUTO_TEST_CASE( test_pch_dynamic_one_element_non_identity_geomap )
{
    BOOST_TEST_MESSAGE( "Single-element affine non-reference simplex: exact manufactured polynomial checks for P1/P2/P3" );

    auto mesh = makeSingleSimplexElementMesh( "pch-dynamic-oneelt-real", Gmsh::GMSH_REAL_DOMAIN );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_EQUAL( mesh->numElements(), 1 );

    auto const [xmin, xmax, ymin, ymax] = singleElementBounds2D( mesh );
    BOOST_CHECK_SMALL( std::abs( xmin - 0.2 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( xmax - 1.4 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( ymin + 0.1 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( ymax - 1.1 ), 1e-12 );

    constexpr double tol = 1e-10;

    checkExactManufacturedPolynomial<1>(
        mesh,
        cst( 1.0 ) + 2.0 * Px() - 3.0 * Py(),
        "real-affine-nonidentity",
        tol );
    checkExactManufacturedPolynomial<2>(
        mesh,
        cst( 1.0 ) + Px() + 2.0 * Py()
        + Px() * Px() + Px() * Py() + 2.0 * Py() * Py(),
        "real-affine-nonidentity",
        tol );
    checkExactManufacturedPolynomial<3>(
        mesh,
        cst( 1.0 ) + Px() + Py()
        + Px() * Px() + Px() * Py() + Py() * Py()
        + Px() * Px() * Px()
        + Px() * Px() * Py()
        + Px() * Py() * Py()
        + Py() * Py() * Py(),
        "real-affine-nonidentity",
        tol );
}

BOOST_AUTO_TEST_CASE( test_pch_dynamic_runtime_orders_up_to_10_reference_simplex )
{
    BOOST_TEST_MESSAGE( "Single-element reference simplex: static vs dynamic consistency for RuntimeOrder 1..10" );

    auto mesh = makeSingleSimplexElementMesh( "pch-dynamic-runtime-1-to-10", Gmsh::GMSH_REFERENCE_DOMAIN );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_EQUAL( mesh->numElements(), 1 );

    constexpr double tol = 1e-10;

    auto exprP1 = expr( makeManufacturedPolynomialString<1>() );
    auto exprP2 = expr( makeManufacturedPolynomialString<2>() );
    auto exprP3 = expr( makeManufacturedPolynomialString<3>() );
    auto exprP4 = expr( makeManufacturedPolynomialString<4>() );
    auto exprP5 = expr( makeManufacturedPolynomialString<5>() );
    auto exprP6 = expr( makeManufacturedPolynomialString<6>() );
    auto exprP7 = expr( makeManufacturedPolynomialString<7>() );
    auto exprP8 = expr( makeManufacturedPolynomialString<8>() );
    auto exprP9 = expr( makeManufacturedPolynomialString<9>() );
    auto exprP10 = expr( makeManufacturedPolynomialString<10>() );

    checkExactManufacturedPolynomial<1>( mesh, exprP1, "runtime-orders-1-to-10", tol );
    checkExactManufacturedPolynomial<2>( mesh, exprP2, "runtime-orders-1-to-10", tol );
    checkExactManufacturedPolynomial<3>( mesh, exprP3, "runtime-orders-1-to-10", tol );
    checkExactManufacturedPolynomial<4>( mesh, exprP4, "runtime-orders-1-to-10", tol );
    checkExactManufacturedPolynomial<5>( mesh, exprP5, "runtime-orders-1-to-10", tol );
    checkExactManufacturedPolynomial<6>( mesh, exprP6, "runtime-orders-1-to-10", tol );
    checkExactManufacturedPolynomial<7>( mesh, exprP7, "runtime-orders-1-to-10", tol );
    checkExactManufacturedPolynomial<8>( mesh, exprP8, "runtime-orders-1-to-10", tol );
    checkExactManufacturedPolynomial<9>( mesh, exprP9, "runtime-orders-1-to-10", tol );
    checkExactManufacturedPolynomial<10>( mesh, exprP10, "runtime-orders-1-to-10", tol );
}

BOOST_AUTO_TEST_CASE( test_pdhv_dynamic_one_element_reference_simplex )
{
    BOOST_TEST_MESSAGE( "Single-element reference simplex: static vs dynamic checks for Pdhv at P1/P2/P3" );

    auto mesh = makeSingleSimplexElementMesh( "pch-dynamic-pdhv", Gmsh::GMSH_REFERENCE_DOMAIN );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_EQUAL( mesh->numElements(), 1 );

    constexpr double tol = 1e-10;
    checkPdhvStaticDynamic<1>( mesh, tol );
    checkPdhvStaticDynamic<2>( mesh, tol );
    checkPdhvStaticDynamic<3>( mesh, tol );
}

BOOST_AUTO_TEST_CASE( test_pdhm_dynamic_one_element_reference_simplex )
{
    BOOST_TEST_MESSAGE( "Single-element reference simplex: static vs dynamic checks for Pdhm at P1/P2/P3" );

    auto mesh = makeSingleSimplexElementMesh( "pch-dynamic-pdhm", Gmsh::GMSH_REFERENCE_DOMAIN );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_EQUAL( mesh->numElements(), 1 );

    constexpr double tol = 1e-10;
    checkPdhmStaticDynamic<1>( mesh, tol );
    checkPdhmStaticDynamic<2>( mesh, tol );
    checkPdhmStaticDynamic<3>( mesh, tol );
}

BOOST_AUTO_TEST_CASE( test_pchm_dynamic_one_element_reference_simplex )
{
    BOOST_TEST_MESSAGE( "Single-element reference simplex: static vs dynamic checks for Pchm at P1/P2/P3" );

    auto mesh = makeSingleSimplexElementMesh( "pch-dynamic-pchm", Gmsh::GMSH_REFERENCE_DOMAIN );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_EQUAL( mesh->numElements(), 1 );

    constexpr double tol = 1e-10;
    checkPchmStaticDynamic<1>( mesh, tol );
    checkPchmStaticDynamic<2>( mesh, tol );
    checkPchmStaticDynamic<3>( mesh, tol );
}

/**
 * @brief Test that Lagrange basis has proper is_order_static/dynamic flags
 */
BOOST_AUTO_TEST_CASE( test_lagrange_order_flags )
{
    BOOST_TEST_MESSAGE( "Test Lagrange order type flags" );

    // Static order P2
    using lagrange_p2_t = Lagrange<2, Scalar>;
    static_assert( !lagrange_p2_t::is_order_dynamic, "P2 should not be dynamic" );
    static_assert( lagrange_p2_t::is_order_static, "P2 should be static" );

    // Dynamic order
    using lagrange_dyn_t = Lagrange<Dynamic, Scalar>;
    static_assert( lagrange_dyn_t::is_order_dynamic, "Dynamic should be dynamic" );
    static_assert( !lagrange_dyn_t::is_order_static, "Dynamic should not be static" );

    BOOST_CHECK( true );
}

/**
 * @brief Test static order Pch spaces (baseline)
 */
BOOST_AUTO_TEST_CASE( test_pch_static_orders )
{
    BOOST_TEST_MESSAGE( "Test static order Pch spaces" );

    auto mesh = unitSquare();

    // Test P1
    {
        auto Vh = Pch<1>( mesh );
        BOOST_TEST_MESSAGE( "  P1 nDof = " << Vh->nDof() );
        BOOST_CHECK( Vh->nDof() > 0 );

        auto u = Vh->element();
        u = project( _space = Vh, _expr = Px() + Py() );
        double l2norm = normL2( _range = elements( mesh ), _expr = idv( u ) );
        BOOST_CHECK( l2norm > 0 );
    }

    // Test P2
    {
        auto Vh = Pch<2>( mesh );
        BOOST_TEST_MESSAGE( "  P2 nDof = " << Vh->nDof() );
        BOOST_CHECK( Vh->nDof() > 0 );

        auto u = Vh->element();
        u = project( _space = Vh, _expr = Px() * Py() );
        double l2norm = normL2( _range = elements( mesh ), _expr = idv( u ) );
        BOOST_CHECK( l2norm > 0 );
    }

    // Test P3
    {
        auto Vh = Pch<3>( mesh );
        BOOST_TEST_MESSAGE( "  P3 nDof = " << Vh->nDof() );
        BOOST_CHECK( Vh->nDof() > 0 );
    }
}

/**
 * @brief Test meta::Pch is_dynamic flag
 */
BOOST_AUTO_TEST_CASE( test_meta_pch_dynamic_flag )
{
    BOOST_TEST_MESSAGE( "Test meta::Pch is_dynamic flag" );

    using mesh_type = Mesh<Simplex<2, 1>>;

    // Static order should not be dynamic
    static_assert( !meta::Pch<mesh_type, 2>::is_dynamic,
                   "Static order Pch should not be dynamic" );

    // Dynamic order should be dynamic
    static_assert( meta::Pch<mesh_type, Dynamic>::is_dynamic,
                   "Dynamic order Pch should be dynamic" );

    BOOST_CHECK( true );
}

/**
 * @brief Test runtime DOF accessor methods on fem::Lagrange
 */
BOOST_AUTO_TEST_CASE( test_lagrange_runtime_dof_accessors )
{
    BOOST_TEST_MESSAGE( "Test Lagrange runtime DOF accessor methods" );

    // Create a P2 basis and verify runtime accessors match static values
    using basis_type = fem::Lagrange<2, 2, 2, Scalar, Continuous, double, Simplex, PointSetFekete, 0>;

    basis_type basis;

    BOOST_CHECK_EQUAL( basis.runtimeOrder(), 2 );
    BOOST_CHECK_EQUAL( basis.runtimeLocalDof(), basis_type::nLocalDof );
    BOOST_CHECK_EQUAL( basis.runtimeDofPerVertex(), basis_type::nDofPerVertex );
    BOOST_CHECK_EQUAL( basis.runtimeDofPerEdge(), basis_type::nDofPerEdge );
    BOOST_CHECK_EQUAL( basis.runtimeDofPerFace(), basis_type::nDofPerFace );

    BOOST_TEST_MESSAGE( "  P2 runtimeLocalDof = " << basis.runtimeLocalDof() );
    BOOST_TEST_MESSAGE( "  P2 nLocalDof = " << basis_type::nLocalDof );
}

/**
 * @brief Test Lagrange constructor with RuntimeOrder
 */
BOOST_AUTO_TEST_CASE( test_lagrange_runtime_order_constructor )
{
    BOOST_TEST_MESSAGE( "Test Lagrange RuntimeOrder constructor" );

    // Static order - RuntimeOrder is ignored
    {
        using basis_type = fem::Lagrange<2, 2, 2, Scalar, Continuous, double, Simplex, PointSetFekete, 0>;
        basis_type basis( RuntimeOrder( 5 ) );  // Should be ignored
        BOOST_CHECK_EQUAL( basis.runtimeOrder(), 2 );  // Still P2
    }

    BOOST_CHECK( true );
}

/**
 * @brief Test Pch with unified RuntimeOrder interface
 *
 * Both static and dynamic Pch spaces can be created with RuntimeOrder.
 * For static spaces, the RuntimeOrder is stored but the compile-time
 * order is used for DOF calculations.
 */
BOOST_AUTO_TEST_CASE( test_pch_unified_interface )
{
    BOOST_TEST_MESSAGE( "Test Pch unified RuntimeOrder interface" );

    auto mesh = unitSquare();

    // Static order with explicit RuntimeOrder (should be same as Pch<2>(mesh))
    {
        auto Vh_static = Pch<2>( mesh );
        auto Vh_with_ro = Pch<2>( mesh, RuntimeOrder( 2 ) );

        BOOST_CHECK_EQUAL( Vh_static->nDof(), Vh_with_ro->nDof() );
        BOOST_CHECK_EQUAL( Vh_static->runtimeOrder(), Vh_with_ro->runtimeOrder() );
        BOOST_CHECK_EQUAL( Vh_static->runtimeOrder(), 2 );

        BOOST_TEST_MESSAGE( "  Static P2: nDof = " << Vh_static->nDof()
                            << ", runtimeOrder = " << Vh_static->runtimeOrder() );
    }

    BOOST_CHECK( true );
}

/**
 * @brief Test Pch<Dynamic> creates function space with correct order
 */
BOOST_AUTO_TEST_CASE( test_pch_dynamic_creation )
{
    BOOST_TEST_MESSAGE( "Test Pch<Dynamic> function space creation" );

    auto mesh = unitSquare();

    // Create dynamic order P1
    {
        auto Vh = Pch<Dynamic>( mesh, RuntimeOrder( 1 ) );
        BOOST_CHECK( Vh->nDof() > 0 );
        BOOST_CHECK_EQUAL( Vh->runtimeOrder(), 1 );
        BOOST_TEST_MESSAGE( "  Dynamic P1: nDof = " << Vh->nDof()
                            << ", runtimeOrder = " << Vh->runtimeOrder() );
    }

    // Create dynamic order P2
    {
        auto Vh = Pch<Dynamic>( mesh, RuntimeOrder( 2 ) );
        BOOST_CHECK( Vh->nDof() > 0 );
        BOOST_CHECK_EQUAL( Vh->runtimeOrder(), 2 );
        BOOST_TEST_MESSAGE( "  Dynamic P2: nDof = " << Vh->nDof()
                            << ", runtimeOrder = " << Vh->runtimeOrder() );
    }

    // Create dynamic order P3
    {
        auto Vh = Pch<Dynamic>( mesh, RuntimeOrder( 3 ) );
        BOOST_CHECK( Vh->nDof() > 0 );
        BOOST_CHECK_EQUAL( Vh->runtimeOrder(), 3 );
        BOOST_TEST_MESSAGE( "  Dynamic P3: nDof = " << Vh->nDof()
                            << ", runtimeOrder = " << Vh->runtimeOrder() );
    }
}

/**
 * @brief Test that runtimeOrder() works correctly for static and dynamic Pch
 *
 * Verifies that runtimeOrder() returns the correct polynomial order
 * for both static and dynamic function spaces.
 *
 * @note **Phase 3 Pending**: Full nDof() equality between static and dynamic
 *       spaces requires updating DofTable to use runtime DOF counts from the
 *       basis class. Currently, dynamic spaces use a placeholder DOF structure.
 */
BOOST_AUTO_TEST_CASE( test_pch_static_vs_dynamic_runtime_order )
{
    BOOST_TEST_MESSAGE( "Test static vs dynamic Pch runtimeOrder equality" );

    auto mesh = unitSquare();

    // Compare P1
    {
        auto Vh_static = Pch<1>( mesh );
        auto Vh_dynamic = Pch<Dynamic>( mesh, RuntimeOrder( 1 ) );

        // runtimeOrder should match
        BOOST_CHECK_EQUAL( Vh_static->runtimeOrder(), Vh_dynamic->runtimeOrder() );
        BOOST_CHECK_EQUAL( Vh_static->runtimeOrder(), 1 );

        BOOST_TEST_MESSAGE( "  P1 static runtimeOrder = " << Vh_static->runtimeOrder()
                            << ", dynamic runtimeOrder = " << Vh_dynamic->runtimeOrder() );
    }

    // Compare P2
    {
        auto Vh_static = Pch<2>( mesh );
        auto Vh_dynamic = Pch<Dynamic>( mesh, RuntimeOrder( 2 ) );

        // runtimeOrder should match
        BOOST_CHECK_EQUAL( Vh_static->runtimeOrder(), Vh_dynamic->runtimeOrder() );
        BOOST_CHECK_EQUAL( Vh_static->runtimeOrder(), 2 );

        BOOST_TEST_MESSAGE( "  P2 static runtimeOrder = " << Vh_static->runtimeOrder()
                            << ", dynamic runtimeOrder = " << Vh_dynamic->runtimeOrder() );
    }

    // Compare P3
    {
        auto Vh_static = Pch<3>( mesh );
        auto Vh_dynamic = Pch<Dynamic>( mesh, RuntimeOrder( 3 ) );

        // runtimeOrder should match
        BOOST_CHECK_EQUAL( Vh_static->runtimeOrder(), Vh_dynamic->runtimeOrder() );
        BOOST_CHECK_EQUAL( Vh_static->runtimeOrder(), 3 );

        BOOST_TEST_MESSAGE( "  P3 static runtimeOrder = " << Vh_static->runtimeOrder()
                            << ", dynamic runtimeOrder = " << Vh_dynamic->runtimeOrder() );
    }
}

/**
 * @brief Test projection with dynamic order Pch
 *
 * @note **Phase 3 Pending**: Currently dynamic order spaces use P1 DOF structure
 *       as a placeholder. Full dynamic order support requires DofTable updates.
 */
BOOST_AUTO_TEST_CASE( test_pch_dynamic_projection )
{
    BOOST_TEST_MESSAGE( "Test projection with dynamic order Pch" );

    auto mesh = unitSquare();

    // Create dynamic P1 space (P1 works correctly since placeholder DOF is P1)
    auto Vh = Pch<Dynamic>( mesh, RuntimeOrder( 1 ) );
    BOOST_CHECK( Vh->nDof() > 0 );
    BOOST_CHECK_EQUAL( Vh->runtimeOrder(), 1 );

    // Project a function
    auto u = Vh->element();
    u = project( _space = Vh, _expr = Px() * Py() );

    // Compute L2 norm
    double l2norm = normL2( _range = elements( mesh ), _expr = idv( u ) );
    BOOST_CHECK( l2norm > 0 );

    BOOST_TEST_MESSAGE( "  Dynamic P1 projection: L2 norm = " << l2norm );
}

/**
 * @brief Test that static and dynamic P1 projection give same result
 *
 * P1 is the placeholder order for dynamic spaces, so this test verifies
 * the dynamic infrastructure works correctly for P1.
 */
BOOST_AUTO_TEST_CASE( test_pch_static_vs_dynamic_p1_projection )
{
    BOOST_TEST_MESSAGE( "Test static vs dynamic P1 Pch projection gives same result" );

    auto mesh = unitSquare();

    // Create static and dynamic P1 spaces
    auto Vh_static = Pch<1>( mesh );
    auto Vh_dynamic = Pch<Dynamic>( mesh, RuntimeOrder( 1 ) );

    // Both should have same nDof for P1
    BOOST_CHECK_EQUAL( Vh_static->nDof(), Vh_dynamic->nDof() );

    // Project same expression
    auto u_static = Vh_static->element();
    auto u_dynamic = Vh_dynamic->element();

    u_static = project( _space = Vh_static, _expr = Px() * Px() + Py() * Py() );
    u_dynamic = project( _space = Vh_dynamic, _expr = Px() * Px() + Py() * Py() );

    // Compute L2 norms - should be identical for P1
    double l2_static = normL2( _range = elements( mesh ), _expr = idv( u_static ) );
    double l2_dynamic = normL2( _range = elements( mesh ), _expr = idv( u_dynamic ) );

    BOOST_CHECK_CLOSE( l2_static, l2_dynamic, 1e-10 );

    BOOST_TEST_MESSAGE( "  Static P1 L2 norm = " << l2_static );
    BOOST_TEST_MESSAGE( "  Dynamic P1 L2 norm = " << l2_dynamic );
}

/**
 * @brief Test that static and dynamic P2 projection give same result
 *
 * This verifies the dynamic infrastructure works correctly for P2,
 * which requires runtime polynomial evaluation (order > compile-time order 1).
 */
BOOST_AUTO_TEST_CASE( test_pch_static_vs_dynamic_p2_projection )
{
    BOOST_TEST_MESSAGE( "Test static vs dynamic P2 Pch projection gives same result" );

    auto mesh = unitSquare();

    // Create static and dynamic P2 spaces
    auto Vh_static = Pch<2>( mesh );
    auto Vh_dynamic = Pch<Dynamic>( mesh, RuntimeOrder( 2 ) );

    // Both should have same nDof for P2
    BOOST_CHECK_EQUAL( Vh_static->nDof(), Vh_dynamic->nDof() );

    // Project same expression (quadratic - should be exactly representable in P2)
    auto u_static = Vh_static->element();
    auto u_dynamic = Vh_dynamic->element();

    u_static = project( _space = Vh_static, _expr = Px() * Px() + Py() * Py() );
    u_dynamic = project( _space = Vh_dynamic, _expr = Px() * Px() + Py() * Py() );

    // Compute L2 norms - should be identical for P2
    double l2_static = normL2( _range = elements( mesh ), _expr = idv( u_static ) );
    double l2_dynamic = normL2( _range = elements( mesh ), _expr = idv( u_dynamic ) );

    BOOST_CHECK_CLOSE( l2_static, l2_dynamic, 1e-10 );

    BOOST_TEST_MESSAGE( "  Static P2 L2 norm = " << l2_static );
    BOOST_TEST_MESSAGE( "  Dynamic P2 L2 norm = " << l2_dynamic );
}

/**
 * @brief Test that static and dynamic P3 projection give same result
 *
 * This verifies the dynamic infrastructure works correctly for P3,
 * which requires runtime polynomial evaluation (order > compile-time order 1).
 */
BOOST_AUTO_TEST_CASE( test_pch_static_vs_dynamic_p3_projection )
{
    BOOST_TEST_MESSAGE( "Test static vs dynamic P3 Pch projection gives same result" );

    auto mesh = unitSquare();

    // Create static and dynamic P3 spaces
    auto Vh_static = Pch<3>( mesh );
    auto Vh_dynamic = Pch<Dynamic>( mesh, RuntimeOrder( 3 ) );

    // Both should have same nDof for P3
    BOOST_CHECK_EQUAL( Vh_static->nDof(), Vh_dynamic->nDof() );

    // Project same expression (cubic - should be exactly representable in P3)
    auto u_static = Vh_static->element();
    auto u_dynamic = Vh_dynamic->element();

    u_static = project( _space = Vh_static, _expr = Px() * Px() * Px() + Py() * Py() * Py() );
    u_dynamic = project( _space = Vh_dynamic, _expr = Px() * Px() * Px() + Py() * Py() * Py() );

    // Compute L2 norms - should be identical for P3
    double l2_static = normL2( _range = elements( mesh ), _expr = idv( u_static ) );
    double l2_dynamic = normL2( _range = elements( mesh ), _expr = idv( u_dynamic ) );

    BOOST_CHECK_CLOSE( l2_static, l2_dynamic, 1e-10 );

    BOOST_TEST_MESSAGE( "  Static P3 L2 norm = " << l2_static );
    BOOST_TEST_MESSAGE( "  Dynamic P3 L2 norm = " << l2_dynamic );
}

/**
 * @brief Debug test: Compare Dubiner static vs runtime evaluation directly
 *
 * This test isolates whether the bug is in Dubiner's evaluateRuntime.
 */
BOOST_AUTO_TEST_CASE( test_dubiner_static_vs_runtime_evaluation )
{
    BOOST_TEST_MESSAGE( "Debug: Compare Dubiner static vs runtime evaluation" );

    using namespace Feel;

    // Test points - use equispaced points for P2 on reference triangle
    using points_type = ublas::matrix<double>;
    points_type pts( 2, 6 );  // 2D, 6 points for P2

    // P2 equispaced nodes on reference simplex [-1,1]^2 with y <= -x
    // Vertices
    pts( 0, 0 ) = -1.0; pts( 1, 0 ) = -1.0;  // (−1,−1)
    pts( 0, 1 ) =  1.0; pts( 1, 1 ) = -1.0;  // ( 1,−1)
    pts( 0, 2 ) = -1.0; pts( 1, 2 ) =  1.0;  // (−1, 1)
    // Edge midpoints
    pts( 0, 3 ) =  0.0; pts( 1, 3 ) = -1.0;  // ( 0,−1)
    pts( 0, 4 ) =  0.0; pts( 1, 4 ) =  0.0;  // ( 0, 0)
    pts( 0, 5 ) = -1.0; pts( 1, 5 ) =  0.0;  // (−1, 0)

    // Static P2 Dubiner evaluation
    using dubiner_p2_t = Dubiner<2, 2, 2, Normalized<true>, double>;
    auto eval_static = dubiner_p2_t::evaluate( pts );

    BOOST_TEST_MESSAGE( "  Static P2 Dubiner evaluation shape: "
                        << eval_static.size1() << " x " << eval_static.size2() );

    // Dynamic Dubiner: use CompileTimeOrder=1 with runtime order=2
    // This is how OrthonormalPolynomialSet<...,Dynamic,...> uses Dubiner
    // When order(2) > nOrder(1), evaluate() calls evaluateRuntime()
    using dubiner_dyn_t = Dubiner<2, 2, 1, Normalized<true>, double>;
    auto eval_runtime = dubiner_dyn_t::evaluate( pts, uint16_type( 2 ) );

    BOOST_TEST_MESSAGE( "  Runtime P2 Dubiner evaluation shape: "
                        << eval_runtime.size1() << " x " << eval_runtime.size2() );

    // Check dimensions match
    BOOST_CHECK_EQUAL( eval_static.size1(), eval_runtime.size1() );
    BOOST_CHECK_EQUAL( eval_static.size2(), eval_runtime.size2() );

    // Compare values
    bool all_match = true;
    double max_diff = 0.0;
    for ( size_t i = 0; i < eval_static.size1(); ++i )
    {
        for ( size_t j = 0; j < eval_static.size2(); ++j )
        {
            double diff = std::abs( eval_static( i, j ) - eval_runtime( i, j ) );
            max_diff = std::max( max_diff, diff );
            if ( diff > 1e-10 )
            {
                all_match = false;
                BOOST_TEST_MESSAGE( "  Mismatch at (" << i << "," << j << "): "
                                    << "static=" << eval_static( i, j )
                                    << " runtime=" << eval_runtime( i, j )
                                    << " diff=" << diff );
            }
        }
    }

    BOOST_TEST_MESSAGE( "  Max difference: " << max_diff );
    BOOST_CHECK_MESSAGE( all_match, "Dubiner static vs runtime evaluation should match" );

    // Print full matrices for debugging
    BOOST_TEST_MESSAGE( "  Static evaluation matrix:" );
    for ( size_t i = 0; i < eval_static.size1(); ++i )
    {
        std::ostringstream oss;
        oss << "    [" << i << "]: ";
        for ( size_t j = 0; j < eval_static.size2(); ++j )
            oss << eval_static( i, j ) << " ";
        BOOST_TEST_MESSAGE( oss.str() );
    }

    BOOST_TEST_MESSAGE( "  Runtime evaluation matrix:" );
    for ( size_t i = 0; i < eval_runtime.size1(); ++i )
    {
        std::ostringstream oss;
        oss << "    [" << i << "]: ";
        for ( size_t j = 0; j < eval_runtime.size2(); ++j )
            oss << eval_runtime( i, j ) << " ";
        BOOST_TEST_MESSAGE( oss.str() );
    }
}

/**
 * @brief Debug test: Test Pch function space construction step by step
 *
 * This test traces through the Pch construction to identify memory corruption.
 */
BOOST_AUTO_TEST_CASE( test_pch_construction_debug )
{
    BOOST_TEST_MESSAGE( "Debug: Test Pch function space construction with element and projection" );

    using namespace Feel;

    // Create mesh
    BOOST_TEST_MESSAGE( "  Step 1: Creating mesh..." );
    auto mesh = unitSquare();
    BOOST_TEST_MESSAGE( "  Mesh created" );

    // Create static P2 space
    BOOST_TEST_MESSAGE( "  Step 2: Creating static P2 space..." );
    auto Vh_static = Pch<2>( mesh );
    BOOST_TEST_MESSAGE( "  Static P2 space created, nDof = " << Vh_static->nDof() );

    // Create dynamic P2 space
    BOOST_TEST_MESSAGE( "  Step 3: Creating dynamic P2 space..." );
    auto Vh_dynamic = Pch<Dynamic>( mesh, RuntimeOrder( 2 ) );
    BOOST_TEST_MESSAGE( "  Dynamic P2 space created, nDof = " << Vh_dynamic->nDof() );

    // Compare nDof
    BOOST_CHECK_EQUAL( Vh_static->nDof(), Vh_dynamic->nDof() );

    // Create elements
    BOOST_TEST_MESSAGE( "  Step 4: Creating static element..." );
    auto u_static = Vh_static->element();
    BOOST_TEST_MESSAGE( "  Static element created" );

    BOOST_TEST_MESSAGE( "  Step 5: Creating dynamic element..." );
    auto u_dynamic = Vh_dynamic->element();
    BOOST_TEST_MESSAGE( "  Dynamic element created" );

    // Project - static
    BOOST_TEST_MESSAGE( "  Step 6: Projecting to static element..." );
    u_static = project( _space = Vh_static, _expr = Px() * Px() + Py() * Py() );
    BOOST_TEST_MESSAGE( "  Static projection done" );

    // Project - dynamic
    BOOST_TEST_MESSAGE( "  Step 7: Projecting to dynamic element..." );
    u_dynamic = project( _space = Vh_dynamic, _expr = Px() * Px() + Py() * Py() );
    BOOST_TEST_MESSAGE( "  Dynamic projection done" );

    // Compute L2 norms - separate steps to find where crash occurs
    BOOST_TEST_MESSAGE( "  Step 8a: Computing static L2 norm..." );
    double l2_static = normL2( _range = elements( mesh ), _expr = idv( u_static ) );
    BOOST_TEST_MESSAGE( "  Static L2 norm = " << l2_static );

    BOOST_TEST_MESSAGE( "  Step 8b: Computing dynamic L2 norm..." );
    double l2_dynamic = normL2( _range = elements( mesh ), _expr = idv( u_dynamic ) );
    BOOST_TEST_MESSAGE( "  Dynamic L2 norm = " << l2_dynamic );

    BOOST_CHECK_CLOSE( l2_static, l2_dynamic, 1e-10 );

    BOOST_TEST_MESSAGE( "  Test completed successfully!" );
}

/**
 * @brief Debug test: Test full Lagrange finite element construction
 *
 * This test traces through the Lagrange construction step by step to identify
 * where memory corruption occurs.
 */
BOOST_AUTO_TEST_CASE( test_lagrange_full_construction_debug )
{
    BOOST_TEST_MESSAGE( "Debug: Test full Lagrange finite element construction" );

    using namespace Feel;

    // Static P2 Lagrange - should work
    BOOST_TEST_MESSAGE( "  Step 1: Creating static P2 Lagrange..." );
    using lagrange_p2_static_t = fem::Lagrange<2, 2, 2, Scalar, Continuous, double, Simplex, PointSetEquiSpaced, 0>;
    lagrange_p2_static_t lagrange_static;
    BOOST_TEST_MESSAGE( "  Static P2 Lagrange created successfully" );
    BOOST_TEST_MESSAGE( "  Static P2 nDof = " << lagrange_static.nDof );
    BOOST_TEST_MESSAGE( "  Static P2 runtimeLocalDof() = " << lagrange_static.runtimeLocalDof() );
    BOOST_TEST_MESSAGE( "  Static P2 coeff() shape: " << lagrange_static.coeff().size1() << " x " << lagrange_static.coeff().size2() );

    // Dynamic P2 Lagrange - this is where we expect the problem
    BOOST_TEST_MESSAGE( "  Step 2: Creating dynamic P2 Lagrange (RuntimeOrder(2))..." );
    using lagrange_dyn_t = fem::Lagrange<2, 2, Dynamic, Scalar, Continuous, double, Simplex, PointSetEquiSpaced, 0>;

    BOOST_TEST_MESSAGE( "  About to construct..." );
    lagrange_dyn_t lagrange_dynamic( RuntimeOrder( 2 ) );
    BOOST_TEST_MESSAGE( "  Dynamic P2 Lagrange created successfully!" );
    BOOST_TEST_MESSAGE( "  Dynamic P2 runtimeLocalDof() = " << lagrange_dynamic.runtimeLocalDof() );
    BOOST_TEST_MESSAGE( "  Dynamic P2 runtimeOrder() = " << lagrange_dynamic.runtimeOrder() );
    BOOST_TEST_MESSAGE( "  Dynamic P2 coeff() shape: " << lagrange_dynamic.coeff().size1() << " x " << lagrange_dynamic.coeff().size2() );

    // Check that both have the same number of DOFs
    BOOST_CHECK_EQUAL( lagrange_static.runtimeLocalDof(), lagrange_dynamic.runtimeLocalDof() );
    BOOST_CHECK_EQUAL( lagrange_static.coeff().size1(), lagrange_dynamic.coeff().size1() );
    BOOST_CHECK_EQUAL( lagrange_static.coeff().size2(), lagrange_dynamic.coeff().size2() );

    // Try evaluating at some points
    using points_type = ublas::matrix<double>;
    points_type pts( 2, 6 );
    pts( 0, 0 ) = -1.0; pts( 1, 0 ) = -1.0;
    pts( 0, 1 ) =  1.0; pts( 1, 1 ) = -1.0;
    pts( 0, 2 ) = -1.0; pts( 1, 2 ) =  1.0;
    pts( 0, 3 ) =  0.0; pts( 1, 3 ) = -1.0;
    pts( 0, 4 ) =  0.0; pts( 1, 4 ) =  0.0;
    pts( 0, 5 ) = -1.0; pts( 1, 5 ) =  0.0;

    // First, let's look at the raw basis evaluation
    BOOST_TEST_MESSAGE( "  Step 3a: Evaluating static basis at points..." );
    auto basis_eval_static = lagrange_static.basis().evaluate( pts );
    BOOST_TEST_MESSAGE( "  Static basis eval shape: " << basis_eval_static.size1() << " x " << basis_eval_static.size2() );
    BOOST_TEST_MESSAGE( "  Static basis row 0: " );
    {
        std::ostringstream oss;
        for ( size_t j = 0; j < basis_eval_static.size2(); ++j )
            oss << basis_eval_static( 0, j ) << " ";
        BOOST_TEST_MESSAGE( "    " << oss.str() );
    }

    BOOST_TEST_MESSAGE( "  Step 3b: Evaluating dynamic basis at points with runtimeOrder..." );
    auto basis_eval_dynamic = lagrange_dynamic.basis().evaluate( pts, lagrange_dynamic.runtimeOrder() );
    BOOST_TEST_MESSAGE( "  Dynamic basis eval shape: " << basis_eval_dynamic.size1() << " x " << basis_eval_dynamic.size2() );
    BOOST_TEST_MESSAGE( "  Dynamic basis row 0: " );
    {
        std::ostringstream oss;
        for ( size_t j = 0; j < basis_eval_dynamic.size2(); ++j )
            oss << basis_eval_dynamic( 0, j ) << " ";
        BOOST_TEST_MESSAGE( "    " << oss.str() );
    }

    BOOST_TEST_MESSAGE( "  Step 3: Evaluating static Lagrange at points..." );
    auto eval_static = lagrange_static.evaluate( pts );
    BOOST_TEST_MESSAGE( "  Static evaluation shape: " << eval_static.size1() << " x " << eval_static.size2() );

    BOOST_TEST_MESSAGE( "  Step 4: Evaluating dynamic Lagrange at points..." );
    auto eval_dynamic = lagrange_dynamic.evaluate( pts );
    BOOST_TEST_MESSAGE( "  Dynamic evaluation shape: " << eval_dynamic.size1() << " x " << eval_dynamic.size2() );

    // Check dimensions match
    BOOST_CHECK_EQUAL( eval_static.size1(), eval_dynamic.size1() );
    BOOST_CHECK_EQUAL( eval_static.size2(), eval_dynamic.size2() );

    // Compare actual values
    double max_diff = 0.0;
    for ( size_t i = 0; i < eval_static.size1(); ++i )
    {
        for ( size_t j = 0; j < eval_static.size2(); ++j )
        {
            double diff = std::abs( eval_static( i, j ) - eval_dynamic( i, j ) );
            if ( diff > 1e-10 )
            {
                BOOST_TEST_MESSAGE( "  Mismatch at (" << i << "," << j << "): "
                                    << "static=" << eval_static( i, j )
                                    << " dynamic=" << eval_dynamic( i, j ) );
            }
            max_diff = std::max( max_diff, diff );
        }
    }
    BOOST_TEST_MESSAGE( "  Max evaluation difference: " << max_diff );

    // Also compare coefficient matrices
    BOOST_TEST_MESSAGE( "  Static coeff matrix:" );
    for ( size_t i = 0; i < std::min<size_t>( 3, lagrange_static.coeff().size1() ); ++i )
    {
        std::ostringstream oss;
        oss << "    row " << i << ": ";
        for ( size_t j = 0; j < lagrange_static.coeff().size2(); ++j )
            oss << lagrange_static.coeff()( i, j ) << " ";
        BOOST_TEST_MESSAGE( oss.str() );
    }

    BOOST_TEST_MESSAGE( "  Dynamic coeff matrix:" );
    for ( size_t i = 0; i < std::min<size_t>( 3, lagrange_dynamic.coeff().size1() ); ++i )
    {
        std::ostringstream oss;
        oss << "    row " << i << ": ";
        for ( size_t j = 0; j < lagrange_dynamic.coeff().size2(); ++j )
            oss << lagrange_dynamic.coeff()( i, j ) << " ";
        BOOST_TEST_MESSAGE( oss.str() );
    }

    BOOST_CHECK_SMALL( max_diff, 1e-10 );

    BOOST_TEST_MESSAGE( "  Test completed successfully!" );
}

/**
 * @brief Debug test: Compare Lagrange primal space static vs dynamic basisEvaluate
 *
 * This test checks if Lagrange's primal_space_type basisEvaluate gives same results
 * for static P2 vs dynamic P2.
 */
BOOST_AUTO_TEST_CASE( test_lagrange_primal_static_vs_dynamic_basisEvaluate )
{
    BOOST_TEST_MESSAGE( "Debug: Compare Lagrange primal basisEvaluate" );

    using namespace Feel;

    // Test points - use equispaced points for P2 on reference triangle
    using points_type = ublas::matrix<double>;
    points_type pts( 2, 6 );  // 2D, 6 points for P2

    // P2 equispaced nodes on reference simplex [-1,1]^2 with y <= -x
    pts( 0, 0 ) = -1.0; pts( 1, 0 ) = -1.0;
    pts( 0, 1 ) =  1.0; pts( 1, 1 ) = -1.0;
    pts( 0, 2 ) = -1.0; pts( 1, 2 ) =  1.0;
    pts( 0, 3 ) =  0.0; pts( 1, 3 ) = -1.0;
    pts( 0, 4 ) =  0.0; pts( 1, 4 ) =  0.0;
    pts( 0, 5 ) = -1.0; pts( 1, 5 ) =  0.0;

    // Static P2 Lagrange
    using lagrange_p2_static_t = fem::Lagrange<2, 2, 2, Scalar, Continuous, double, Simplex, PointSetEquiSpaced, 0>;
    using primal_static_t = typename lagrange_p2_static_t::primal_space_type;
    primal_static_t primal_static;
    auto eval_static = primal_static.basisEvaluate( pts );

    BOOST_TEST_MESSAGE( "  Static P2 primal basisEvaluate shape: "
                        << eval_static.size1() << " x " << eval_static.size2() );

    // Dynamic P2 Lagrange: <Dim, RealDim, Order, PolySetType, ...>
    using lagrange_dyn_t = fem::Lagrange<2, 2, Dynamic, Scalar, Continuous, double, Simplex, PointSetEquiSpaced, 0>;
    using primal_dyn_t = typename lagrange_dyn_t::primal_space_type;
    primal_dyn_t primal_dynamic( RuntimeOrder( 2 ) );

    BOOST_TEST_MESSAGE( "  Dynamic primal runtimeOrder() = " << primal_dynamic.runtimeOrder() );
    BOOST_TEST_MESSAGE( "  Dynamic primal runtimeLocalDof() = " << primal_dynamic.runtimeLocalDof() );
    BOOST_TEST_MESSAGE( "  Dynamic primal is_order_dynamic = " << primal_dyn_t::is_order_dynamic );

    auto eval_dynamic = primal_dynamic.basisEvaluate( pts );

    BOOST_TEST_MESSAGE( "  Dynamic P2 primal basisEvaluate shape: "
                        << eval_dynamic.size1() << " x " << eval_dynamic.size2() );

    // Check dimensions match
    BOOST_CHECK_EQUAL( eval_static.size1(), eval_dynamic.size1() );
    BOOST_CHECK_EQUAL( eval_static.size2(), eval_dynamic.size2() );

    // Compare values
    bool all_match = true;
    double max_diff = 0.0;
    for ( size_t i = 0; i < eval_static.size1(); ++i )
    {
        for ( size_t j = 0; j < eval_static.size2(); ++j )
        {
            double diff = std::abs( eval_static( i, j ) - eval_dynamic( i, j ) );
            max_diff = std::max( max_diff, diff );
            if ( diff > 1e-10 )
            {
                all_match = false;
                BOOST_TEST_MESSAGE( "  Mismatch at (" << i << "," << j << "): "
                                    << "static=" << eval_static( i, j )
                                    << " dynamic=" << eval_dynamic( i, j )
                                    << " diff=" << diff );
            }
        }
    }

    BOOST_TEST_MESSAGE( "  Max difference: " << max_diff );
    BOOST_CHECK_MESSAGE( all_match, "Lagrange primal basisEvaluate static vs dynamic should match" );
}

BOOST_AUTO_TEST_SUITE_END()
