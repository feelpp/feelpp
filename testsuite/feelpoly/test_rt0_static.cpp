/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-02-17

  Copyright (C) 2026 Feel++ Consortium

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
/**
 * @file test_rt0_static.cpp
 * @brief Compile-time/runtime sanity checks for RT simplex elements
 */
#define BOOST_TEST_MODULE test_rt0_static
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/raviartthomas.hpp>

#include <array>
#include <cmath>
#include <stdexcept>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
struct FakeRtGeom2D
{
    uint16_type faceId() const
    {
        return invalid_uint16_type_value;
    }

    void faceNormal( int face, ublas::vector<double>& normal, bool ) const
    {
        normal.resize( 2, false );
        switch ( face )
        {
        case 0:
            normal( 0 ) = -1.0;
            normal( 1 ) = -1.0;
            break;
        case 1:
            normal( 0 ) = 1.0;
            normal( 1 ) = 0.0;
            break;
        default:
            normal( 0 ) = 0.0;
            normal( 1 ) = 1.0;
            break;
        }
    }
};

struct ConstantVectorExpr2D
{
    struct shape
    {
        static constexpr int M = 2;
    };

    FakeRtGeom2D* geom()
    {
        return &M_geom;
    }

    FakeRtGeom2D const* geom() const
    {
        return &M_geom;
    }

    double evalq( int c1, int, int ) const
    {
        return ( c1 == 0 ) ? 1.25 : -0.75;
    }

private:
    FakeRtGeom2D M_geom;
};

template<typename RTType>
struct ReferenceRtGeom
{
    using value_type = typename RTType::value_type;
    using reference_convex_type = typename RTType::reference_convex_type;

    uint16_type faceId() const
    {
        return invalid_uint16_type_value;
    }

    void faceNormal( int face, ublas::vector<value_type>& normal, bool ) const
    {
        normal.resize( RTType::nDim, false );
        em_node_type<value_type> mappedNormal( normal.data().begin(), normal.size() );
        if constexpr ( RTType::nDim == 2 )
        {
            static constexpr std::array<value_type, 3> scaling = { value_type( 2.8284271247461903 ), value_type( 2.0 ), value_type( 2.0 ) };
            mappedNormal = M_ref.normal( face ) * scaling[face];
        }
        else
        {
            static constexpr std::array<value_type, 4> scaling = { value_type( 3.464101615137754 ), value_type( 2.0 ), value_type( 2.0 ), value_type( 2.0 ) };
            mappedNormal = M_ref.normal( face ) * scaling[face];
        }
    }

private:
    reference_convex_type M_ref;
};

template<typename RTType, int Order>
struct ReferenceRtExpr
{
    using value_type = typename RTType::value_type;
    using points_type = typename RTType::points_type;

    struct shape
    {
        static constexpr int M = RTType::nDim;
    };

    explicit ReferenceRtExpr( points_type const& pts )
        :
        M_pts( pts )
    {}

    ReferenceRtGeom<RTType>* geom()
    {
        return &M_geom;
    }

    ReferenceRtGeom<RTType> const* geom() const
    {
        return &M_geom;
    }

    int nPoints() const
    {
        return M_pts.size2();
    }

    value_type evalq( int c1, int, int q ) const
    {
        return value( c1, q );
    }

    value_type value( int c1, int q ) const
    {
        const value_type x = M_pts( 0, q );
        const value_type y = ( RTType::nDim > 1 ) ? M_pts( 1, q ) : value_type( 0 );
        const value_type z = ( RTType::nDim > 2 ) ? M_pts( 2, q ) : value_type( 0 );
        const value_type coord = ( c1 == 0 ) ? x : ( c1 == 1 ? y : z );

        if constexpr ( Order == 0 )
            return ( c1 == 0 ? value_type( 1.25 ) : ( c1 == 1 ? value_type( -0.75 ) : value_type( 0.5 ) ) ) + value_type( 0.125 ) * coord;
        else if constexpr ( Order == 1 )
        {
            const value_type h = x;
            if ( c1 == 0 )
                return value_type( 1.0 ) + value_type( 2.0 )*x - value_type( 0.25 )*y + coord*h;
            if ( c1 == 1 )
                return value_type( -0.5 ) + x + value_type( 3.0 )*y + coord*h;
            return value_type( 0.25 ) - y + value_type( 0.5 )*z + coord*h;
        }
        else if constexpr ( Order == 2 )
        {
            const value_type h = x*y;
            if ( c1 == 0 )
                return value_type( 1.0 ) + x + y + x*y + coord*h;
            if ( c1 == 1 )
                return value_type( -0.5 ) + x*x + y*y + coord*h;
            return value_type( 0.25 ) + z + x*z + coord*h;
        }
        else
        {
            const value_type h = x*x*x;
            if ( c1 == 0 )
                return value_type( 1.0 ) + x + y + x*y + x*x*y + coord*h;
            if ( c1 == 1 )
                return value_type( -0.5 ) + x*x + y*y + x*y*y + coord*h;
            return value_type( 0.25 ) + z + x*z + x*y*z + coord*h;
        }
    }

private:
    points_type M_pts;
    ReferenceRtGeom<RTType> M_geom;
};

template<typename RTType>
void
checkFiniteEvaluationAtNodes( RTType const& rt )
{
    auto const pts = rt.points();
    auto const values = rt.evaluate( pts );
    for ( std::size_t i = 0; i < values.size1(); ++i )
        for ( std::size_t j = 0; j < values.size2(); ++j )
            BOOST_CHECK( std::isfinite( values( i, j ) ) );
}

template<typename RTType, int Order>
void
checkReferenceRtInterpolationExactness( RTType const& rt )
{
    ReferenceRtExpr<RTType, Order> dofExpr( rt.points() );
    auto Ihloc = rt.localInterpolant();
    rt.interpolate( dofExpr, Ihloc );

    IMGeneral<RTType::nDim, typename RTType::value_type, Simplex> im( 12 );
    ReferenceRtExpr<RTType, Order> exactExpr( im.points() );
    auto const basisAtQuad = rt.evaluate( im.points() );

    BOOST_TEST_CONTEXT( "RT reference interpolation dim=" << RTType::nDim << " order=" << Order )
    {
        for ( int q = 0; q < im.nPoints(); ++q )
        {
            for ( int c = 0; c < RTType::nDim; ++c )
            {
                typename RTType::value_type reconstructed = 0;
                for ( uint16_type l = 0; l < rt.localDofPerComponent(); ++l )
                    reconstructed += Ihloc( l ) * basisAtQuad( RTType::nDim*l + c, q );
                BOOST_CHECK_SMALL( reconstructed - exactExpr.value( c, q ), typename RTType::value_type( 1e-8 ) );
            }
        }
    }
}

template<typename RTType>
uint16_type
runtimeLocalDof( RTType const& rt )
{
    if constexpr ( requires { rt.localDofCount(); } )
        return rt.localDofCount();
    else
        return RTType::nLocalDof;
}

template<typename RTType>
uint16_type
runtimeDofPerVertex( RTType const& rt )
{
    if constexpr ( requires { rt.dofPerVertex(); } )
        return rt.dofPerVertex();
    else
        return RTType::nDofPerVertex;
}

template<typename RTType>
uint16_type
runtimeDofPerEdge( RTType const& rt )
{
    if constexpr ( requires { rt.dofPerEdge(); } )
        return rt.dofPerEdge();
    else
        return RTType::nDofPerEdge;
}

template<typename RTType>
uint16_type
runtimeDofPerFace( RTType const& rt )
{
    if constexpr ( requires { rt.dofPerFace(); } )
        return rt.dofPerFace();
    else
        return RTType::nDofPerFace;
}

template<typename RTType>
uint16_type
runtimeDofPerVolume( RTType const& rt )
{
    if constexpr ( requires { rt.dofPerVolume(); } )
        return rt.dofPerVolume();
    else
        return RTType::nDofPerVolume;
}

template<typename RTType>
void
checkDofAttachmentContract( RTType const& rt )
{
    static_assert( !RTType::is_product );

    std::array<uint16_type, 4> countByDim = { 0, 0, 0, 0 };

    for ( uint16_type ldof = 0; ldof < runtimeLocalDof( rt ); ++ldof )
    {
        BOOST_CHECK_EQUAL( rt.component( ldof ), 0 );
        BOOST_CHECK_EQUAL( rt.dofParent( ldof ), ldof );
        BOOST_CHECK_EQUAL( rt.localDofId( rt.dofParent( ldof ), 0 ), ldof );

        auto const attachment = rt.dofAttachment( ldof );
        BOOST_REQUIRE( attachment.isValid() );
        BOOST_REQUIRE_GE( attachment.entityDim, 0 );
        BOOST_REQUIRE_LE( attachment.entityDim, 3 );
        BOOST_CHECK_EQUAL( attachment.kind, rt.dofType( ldof ) );

        countByDim[attachment.entityDim] += 1;

        switch ( attachment.entityDim )
        {
        case 0:
            BOOST_CHECK_LT( attachment.entityId, RTType::reference_convex_type::numVertices );
            BOOST_CHECK_LT( attachment.ordinal, runtimeDofPerVertex( rt ) );
            break;
        case 1:
            BOOST_CHECK_LT( attachment.entityId, RTType::reference_convex_type::numEdges );
            BOOST_CHECK_LT( attachment.ordinal, runtimeDofPerEdge( rt ) );
            break;
        case 2:
            BOOST_CHECK_LT( attachment.entityId, RTType::reference_convex_type::numFaces );
            BOOST_CHECK_LT( attachment.ordinal, runtimeDofPerFace( rt ) );
            break;
        case 3:
            BOOST_CHECK_EQUAL( attachment.entityId, 0 );
            BOOST_CHECK_LT( attachment.ordinal, runtimeDofPerVolume( rt ) );
            break;
        default:
            BOOST_FAIL( "invalid dof attachment entity dimension" );
        }
    }

    BOOST_CHECK_EQUAL( countByDim[0],
                       static_cast<uint16_type>( RTType::reference_convex_type::numVertices * runtimeDofPerVertex( rt ) ) );
    BOOST_CHECK_EQUAL( countByDim[1],
                       static_cast<uint16_type>( RTType::reference_convex_type::numEdges * runtimeDofPerEdge( rt ) ) );
    BOOST_CHECK_EQUAL( countByDim[2],
                       static_cast<uint16_type>( RTType::reference_convex_type::numFaces * runtimeDofPerFace( rt ) ) );
    BOOST_CHECK_EQUAL( countByDim[3], runtimeDofPerVolume( rt ) );
}

template<int Order>
constexpr uint16_type rtSimplex2dDof()
{
    return static_cast<uint16_type>( ( Order + 1 ) * ( Order + 3 ) );
}

template<int Order>
constexpr uint16_type rtSimplex3dDof()
{
    return static_cast<uint16_type>( ( Order + 1 ) * ( Order + 2 ) * ( Order + 4 ) / 2 );
}

template<int Order>
void
checkRtSimplex2d()
{
    using rt_t = typename RaviartThomas<Order>::template apply<2>::type;

    static_assert( rt_t::nDim == 2 );
    static_assert( rt_t::nOrder == Order + 1 );
    static_assert( rt_t::nLocalDof == rtSimplex2dDof<Order>() );
    static_assert( rt_t::nDofPerVertex == 0 );
    static_assert( rt_t::nDofPerEdge == Order + 1 );
    static_assert( rt_t::nDofPerFace == Order * ( Order + 1 ) );
    static_assert( rt_t::nDofPerVolume == 0 );

    rt_t rt;
    BOOST_CHECK_EQUAL( rt.familyName(), "raviartthomas" );
    BOOST_CHECK_EQUAL( rt.runtimeOrder(), Order + 1 );
    BOOST_CHECK_EQUAL( rt.localDofPerComponent(), rtSimplex2dDof<Order>() );
    if constexpr ( Order == 0 )
        BOOST_CHECK_EQUAL( rt.nbPoints(), rtSimplex2dDof<Order>() );
    else
    {
        IM<2, 2*( rt_t::nOrder + 1 ), double, Simplex> im;
        constexpr uint16_type firstInternalDof = rt_t::reference_convex_type::numEdges * rt_t::nDofPerEdge;
        BOOST_CHECK_GE( rt.nbPoints(), static_cast<uint16_type>( firstInternalDof + im.nPoints() ) );
    }

    checkFiniteEvaluationAtNodes( rt );
    checkDofAttachmentContract( rt );
    checkReferenceRtInterpolationExactness<rt_t, Order>( rt );
}

template<int Order>
void
checkRtSimplex3d()
{
    using rt_t = typename RaviartThomas<Order>::template apply<3>::type;

    static_assert( rt_t::nDim == 3 );
    static_assert( rt_t::nOrder == Order + 1 );
    static_assert( rt_t::nLocalDof == rtSimplex3dDof<Order>() );
    static_assert( rt_t::nDofPerVertex == 0 );
    static_assert( rt_t::nDofPerEdge == 0 );
    static_assert( rt_t::nDofPerFace == ( Order + 1 ) * ( Order + 2 ) / 2 );
    static_assert( rt_t::nDofPerVolume ==
                   rtSimplex3dDof<Order>() - 4 * ( ( Order + 1 ) * ( Order + 2 ) / 2 ) );

    rt_t rt;
    BOOST_CHECK_EQUAL( rt.familyName(), "raviartthomas" );
    BOOST_CHECK_EQUAL( rt.runtimeOrder(), Order + 1 );
    BOOST_CHECK_EQUAL( rt.localDofPerComponent(), rtSimplex3dDof<Order>() );
    if constexpr ( Order == 0 )
        BOOST_CHECK_EQUAL( rt.nbPoints(), rtSimplex3dDof<Order>() );
    else
    {
        IM<3, 2*( rt_t::nOrder + 1 ), double, Simplex> im;
        constexpr uint16_type firstInternalDof = rt_t::reference_convex_type::numTopologicalFaces * rt_t::nDofPerFace;
        BOOST_CHECK_GE( rt.nbPoints(), static_cast<uint16_type>( firstInternalDof + im.nPoints() ) );
    }

    checkFiniteEvaluationAtNodes( rt );
    checkDofAttachmentContract( rt );
    checkReferenceRtInterpolationExactness<rt_t, Order>( rt );
}

template<int Order>
void
checkRtSimplex2dInteriorInterpolation()
{
    using rt_t = typename RaviartThomas<Order>::template apply<2>::type;
    static_assert( Order > 0 );
    static_assert( rt_t::nDofPerFace > 0 );

    rt_t rt;
    ConstantVectorExpr2D expr;
    auto Ihloc = rt.localInterpolant();
    rt.interpolate( expr, Ihloc );

    constexpr uint16_type firstInternalDof = rt_t::reference_convex_type::numEdges * rt_t::nDofPerEdge;
    double interiorNorm = 0.0;
    for ( uint16_type l = 0; l < rt_t::nDofPerFace; ++l )
        interiorNorm += std::abs( Ihloc( firstInternalDof + l ) );

    BOOST_CHECK_GT( interiorNorm, 1e-12 );
}

template<int Order>
void
checkRtDynamicSimplex2d()
{
    using rt_t = typename RaviartThomas<Dynamic>::template apply<2>::type;

    rt_t rt{ RuntimeOrder{ Order } };
    BOOST_CHECK_EQUAL( rt.familyName(), "raviartthomas" );
    BOOST_CHECK_EQUAL( rt.order(), Order );
    BOOST_CHECK_EQUAL( rt.runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( rt.internalOrder(), Order + 1 );
    BOOST_CHECK_EQUAL( rt.localDofPerComponent(), rtSimplex2dDof<Order>() );
    BOOST_CHECK_EQUAL( rt.localDofCount(), rtSimplex2dDof<Order>() );
    BOOST_CHECK_EQUAL( rt.dofPerVertex(), 0 );
    BOOST_CHECK_EQUAL( rt.dofPerEdge(), Order + 1 );
    BOOST_CHECK_EQUAL( rt.dofPerFace(), Order * ( Order + 1 ) );
    BOOST_CHECK_EQUAL( rt.dofPerVolume(), 0 );

    if constexpr ( Order == 0 )
        BOOST_CHECK_EQUAL( rt.nbPoints(), rtSimplex2dDof<Order>() );
    else
    {
        IMGeneral<2, double, Simplex> im( static_cast<uint16_type>( 2*( rt.internalOrder() + 1 ) ) );
        const uint16_type firstInternalDof = static_cast<uint16_type>( rt_t::reference_convex_type::numEdges * rt.dofPerEdge() );
        BOOST_CHECK_GE( rt.nbPoints(), static_cast<uint16_type>( firstInternalDof + im.nPoints() ) );
    }

    checkFiniteEvaluationAtNodes( rt );
    checkDofAttachmentContract( rt );
    checkReferenceRtInterpolationExactness<rt_t, Order>( rt );
}

template<int Order>
void
checkRtDynamicSimplex3d()
{
    using rt_t = typename RaviartThomas<Dynamic>::template apply<3>::type;

    rt_t rt{ RuntimeOrder{ Order } };
    BOOST_CHECK_EQUAL( rt.familyName(), "raviartthomas" );
    BOOST_CHECK_EQUAL( rt.order(), Order );
    BOOST_CHECK_EQUAL( rt.runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( rt.internalOrder(), Order + 1 );
    BOOST_CHECK_EQUAL( rt.localDofPerComponent(), rtSimplex3dDof<Order>() );
    BOOST_CHECK_EQUAL( rt.localDofCount(), rtSimplex3dDof<Order>() );
    BOOST_CHECK_EQUAL( rt.dofPerVertex(), 0 );
    BOOST_CHECK_EQUAL( rt.dofPerEdge(), 0 );
    BOOST_CHECK_EQUAL( rt.dofPerFace(), ( Order + 1 ) * ( Order + 2 ) / 2 );
    BOOST_CHECK_EQUAL( rt.dofPerVolume(),
                       rtSimplex3dDof<Order>() - 4 * ( ( Order + 1 ) * ( Order + 2 ) / 2 ) );

    if constexpr ( Order == 0 )
        BOOST_CHECK_EQUAL( rt.nbPoints(), rtSimplex3dDof<Order>() );
    else
    {
        IMGeneral<3, double, Simplex> im( static_cast<uint16_type>( 2*( rt.internalOrder() + 1 ) ) );
        const uint16_type firstInternalDof = static_cast<uint16_type>( rt_t::reference_convex_type::numTopologicalFaces * rt.dofPerFace() );
        BOOST_CHECK_GE( rt.nbPoints(), static_cast<uint16_type>( firstInternalDof + im.nPoints() ) );
    }

    checkFiniteEvaluationAtNodes( rt );
    checkDofAttachmentContract( rt );
    checkReferenceRtInterpolationExactness<rt_t, Order>( rt );
}
} // namespace

BOOST_AUTO_TEST_SUITE( rt0_static_suite )

BOOST_AUTO_TEST_CASE( rt0_2d_compile_time_sanity )
{
    using rt0_2d_t = RaviartThomas<0>::apply<2>::type;

    static_assert( rt0_2d_t::nDim == 2 );
    static_assert( rt0_2d_t::nOrder == 1 );
    static_assert( rt0_2d_t::nLocalDof == 3 );

    rt0_2d_t rt0;
    BOOST_CHECK_EQUAL( rt0.familyName(), "raviartthomas" );
    BOOST_CHECK_EQUAL( rt0.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( rt0.localDofPerComponent(), 3 );

    checkFiniteEvaluationAtNodes( rt0 );
    checkDofAttachmentContract( rt0 );
}

BOOST_AUTO_TEST_CASE( rt0_3d_compile_time_sanity )
{
    using rt0_3d_t = RaviartThomas<0>::apply<3>::type;

    static_assert( rt0_3d_t::nDim == 3 );
    static_assert( rt0_3d_t::nOrder == 1 );
    static_assert( rt0_3d_t::nLocalDof == 4 );

    rt0_3d_t rt0;
    BOOST_CHECK_EQUAL( rt0.familyName(), "raviartthomas" );
    BOOST_CHECK_EQUAL( rt0.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( rt0.localDofPerComponent(), 4 );

    checkFiniteEvaluationAtNodes( rt0 );
    checkDofAttachmentContract( rt0 );
}

BOOST_AUTO_TEST_CASE( rt1_2d_type_contract )
{
    using rt1_2d_t = RaviartThomas<1>::apply<2>::type;
    static_assert( rt1_2d_t::nDim == 2 );
    static_assert( !rt1_2d_t::is_product );
}

BOOST_AUTO_TEST_CASE( rt1_2d_compile_time_sanity )
{
    checkRtSimplex2d<1>();
}

BOOST_AUTO_TEST_CASE( rt2_2d_compile_time_sanity )
{
    checkRtSimplex2d<2>();
}

BOOST_AUTO_TEST_CASE( rt3_2d_compile_time_sanity )
{
    checkRtSimplex2d<3>();
}

BOOST_AUTO_TEST_CASE( rt1_2d_interior_moment_interpolation )
{
    checkRtSimplex2dInteriorInterpolation<1>();
}

BOOST_AUTO_TEST_CASE( rt2_2d_interior_moment_interpolation )
{
    checkRtSimplex2dInteriorInterpolation<2>();
}

BOOST_AUTO_TEST_CASE( rt3_2d_interior_moment_interpolation )
{
    checkRtSimplex2dInteriorInterpolation<3>();
}

BOOST_AUTO_TEST_CASE( rt1_3d_compile_time_sanity )
{
    checkRtSimplex3d<1>();
}

BOOST_AUTO_TEST_CASE( rt2_3d_compile_time_sanity )
{
    checkRtSimplex3d<2>();
}

BOOST_AUTO_TEST_CASE( rt3_3d_compile_time_sanity )
{
    checkRtSimplex3d<3>();
}

BOOST_AUTO_TEST_CASE( rt_dynamic_simplex_orders_sanity )
{
    using rt0_static_2d_t = RaviartThomas<0>::apply<2>::type;
    using rt0_dynamic_2d_t = RaviartThomas<Dynamic>::apply<2>::type;
    using rt0_static_3d_t = RaviartThomas<0>::apply<3>::type;
    using rt0_dynamic_3d_t = RaviartThomas<Dynamic>::apply<3>::type;

    static_assert( rt0_dynamic_2d_t::is_order_dynamic );
    static_assert( rt0_dynamic_2d_t::nOrder_v == Dynamic );
    static_assert( rt0_dynamic_2d_t::nLocalDof == rt0_static_2d_t::nLocalDof );
    static_assert( rt0_dynamic_3d_t::nLocalDof == rt0_static_3d_t::nLocalDof );

    rt0_dynamic_2d_t rt0_2d{ RuntimeOrder{ 0 } };
    rt0_dynamic_3d_t rt0_3d{ RuntimeOrder{ 0 } };

    BOOST_CHECK_EQUAL( rt0_2d.familyName(), "raviartthomas" );
    BOOST_CHECK_EQUAL( rt0_2d.runtimeOrder(), 0 );
    BOOST_CHECK_EQUAL( rt0_2d.localDofPerComponent(), rt0_static_2d_t::nLocalDof );
    BOOST_CHECK_EQUAL( rt0_3d.runtimeOrder(), 0 );
    BOOST_CHECK_EQUAL( rt0_3d.localDofPerComponent(), rt0_static_3d_t::nLocalDof );

    checkFiniteEvaluationAtNodes( rt0_2d );
    checkFiniteEvaluationAtNodes( rt0_3d );
    checkDofAttachmentContract( rt0_2d );
    checkDofAttachmentContract( rt0_3d );

    checkRtDynamicSimplex2d<1>();
    checkRtDynamicSimplex2d<2>();
    checkRtDynamicSimplex2d<3>();
    checkRtDynamicSimplex3d<1>();
    checkRtDynamicSimplex3d<2>();
    checkRtDynamicSimplex3d<3>();
}

BOOST_AUTO_TEST_SUITE_END()
