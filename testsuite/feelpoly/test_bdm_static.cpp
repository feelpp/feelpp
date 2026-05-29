/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-02-18

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
 * @file test_bdm_static.cpp
 * @brief Compile-time/runtime sanity checks for BDM metadata and topology.
 */
#define BOOST_TEST_MODULE test_bdm_static
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/brezzidouglasmarini.hpp>

#include <array>
#include <cmath>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
struct FakeBdmGeom2D
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

struct IndexedVectorExpr2D
{
    struct shape
    {
        static constexpr int M = 2;
    };

    FakeBdmGeom2D* geom()
    {
        return &M_geom;
    }

    FakeBdmGeom2D const* geom() const
    {
        return &M_geom;
    }

    int nPoints() const
    {
        return 512;
    }

    double evalq( int c1, int, int q ) const
    {
        return ( c1 == 0 ) ? 1.0 + 0.125*q : -0.5 + 0.25*q;
    }

private:
    FakeBdmGeom2D M_geom;
};

template<typename BDMType>
struct ReferenceBdmGeom
{
    using value_type = typename BDMType::value_type;
    using reference_convex_type = typename BDMType::reference_convex_type;

    uint16_type faceId() const
    {
        return invalid_uint16_type_value;
    }

    void faceNormal( int face, ublas::vector<value_type>& normal, bool ) const
    {
        normal.resize( BDMType::nDim, false );
        em_node_type<value_type> mappedNormal( normal.data().begin(), normal.size() );
        if constexpr ( BDMType::nDim == 2 )
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

template<typename BDMType, int Order>
struct ReferenceBdmExpr
{
    using value_type = typename BDMType::value_type;
    using points_type = typename BDMType::points_type;

    struct shape
    {
        static constexpr int M = BDMType::nDim;
    };

    explicit ReferenceBdmExpr( points_type const& pts )
        :
        M_pts( pts )
    {}

    ReferenceBdmGeom<BDMType>* geom()
    {
        return &M_geom;
    }

    ReferenceBdmGeom<BDMType> const* geom() const
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
        const value_type y = ( BDMType::nDim > 1 ) ? M_pts( 1, q ) : value_type( 0 );
        const value_type z = ( BDMType::nDim > 2 ) ? M_pts( 2, q ) : value_type( 0 );

        if constexpr ( Order == 0 )
        {
            if ( c1 == 0 )
                return value_type( 1.0 ) + value_type( 0.5 )*x - y;
            if ( c1 == 1 )
                return value_type( -0.25 ) + value_type( 2.0 )*y + value_type( 0.125 )*x;
            return value_type( 0.5 ) - x + value_type( 1.5 )*z;
        }
        else if constexpr ( Order == 1 )
        {
            if ( c1 == 0 )
                return value_type( 0.75 ) + x + x*y;
            if ( c1 == 1 )
                return value_type( -1.25 ) + y + x*x;
            return value_type( 0.5 ) + z + x*z;
        }
        else if constexpr ( Order == 2 )
        {
            if ( c1 == 0 )
                return value_type( 0.75 ) + x + y + x*y + x*x*y;
            if ( c1 == 1 )
                return value_type( -1.25 ) + y + x*x + y*y + x*y*y;
            return value_type( 0.5 ) + z + x*z + y*z + x*y*z;
        }
        else
        {
            if ( c1 == 0 )
                return value_type( 0.75 ) + x + y + x*y + x*x*y + x*x*y*y;
            if ( c1 == 1 )
                return value_type( -1.25 ) + y + x*x + y*y + x*y*y + x*x*x*y;
            return value_type( 0.5 ) + z + x*z + y*z + x*y*z + x*x*y*z;
        }
    }

private:
    points_type M_pts;
    ReferenceBdmGeom<BDMType> M_geom;
};

template<typename BDMType>
uint16_type
runtimeLocalDof( BDMType const& bdm )
{
    if constexpr ( requires { bdm.localDofCount(); } )
        return bdm.localDofCount();
    else
        return BDMType::nLocalDof;
}

template<typename BDMType>
uint16_type
runtimeDofPerVertex( BDMType const& bdm )
{
    if constexpr ( requires { bdm.dofPerVertex(); } )
        return bdm.dofPerVertex();
    else
        return BDMType::nDofPerVertex;
}

template<typename BDMType>
uint16_type
runtimeDofPerEdge( BDMType const& bdm )
{
    if constexpr ( requires { bdm.dofPerEdge(); } )
        return bdm.dofPerEdge();
    else
        return BDMType::nDofPerEdge;
}

template<typename BDMType>
uint16_type
runtimeDofPerFace( BDMType const& bdm )
{
    if constexpr ( requires { bdm.dofPerFace(); } )
        return bdm.dofPerFace();
    else
        return BDMType::nDofPerFace;
}

template<typename BDMType>
uint16_type
runtimeDofPerVolume( BDMType const& bdm )
{
    if constexpr ( requires { bdm.dofPerVolume(); } )
        return bdm.dofPerVolume();
    else
        return BDMType::nDofPerVolume;
}

template<typename BDMType>
void
checkFiniteEvaluationAtNodes( BDMType const& bdm )
{
    auto const pts = bdm.points();
    auto const values = bdm.evaluate( pts );
    for ( std::size_t i = 0; i < values.size1(); ++i )
        for ( std::size_t j = 0; j < values.size2(); ++j )
            BOOST_CHECK( std::isfinite( values( i, j ) ) );
}

template<typename BDMType, int Order>
void
checkReferenceBdmInterpolationExactness( BDMType const& bdm )
{
    ReferenceBdmExpr<BDMType, Order> dofExpr( bdm.points() );
    auto Ihloc = bdm.localInterpolant();
    bdm.interpolate( dofExpr, Ihloc );

    IMGeneral<BDMType::nDim, typename BDMType::value_type, Simplex> im( 12 );
    ReferenceBdmExpr<BDMType, Order> exactExpr( im.points() );
    auto const basisAtQuad = bdm.evaluate( im.points() );

    BOOST_TEST_CONTEXT( "BDM reference interpolation dim=" << BDMType::nDim << " order=" << Order )
    {
        for ( int q = 0; q < im.nPoints(); ++q )
        {
            for ( int c = 0; c < BDMType::nDim; ++c )
            {
                typename BDMType::value_type reconstructed = 0;
                for ( uint16_type l = 0; l < bdm.localDofPerComponent(); ++l )
                    reconstructed += Ihloc( l ) * basisAtQuad( BDMType::nDim*l + c, q );
                BOOST_CHECK_SMALL( reconstructed - exactExpr.value( c, q ), typename BDMType::value_type( 1e-8 ) );
            }
        }
    }
}

template<typename BDMType>
void
checkDofAttachmentContract( BDMType const& bdm )
{
    static_assert( !BDMType::is_product );

    std::array<uint16_type, 4> countByDim = { 0, 0, 0, 0 };

    for ( uint16_type ldof = 0; ldof < runtimeLocalDof( bdm ); ++ldof )
    {
        BOOST_CHECK_EQUAL( bdm.component( ldof ), 0 );
        BOOST_CHECK_EQUAL( bdm.dofParent( ldof ), ldof );
        BOOST_CHECK_EQUAL( bdm.localDofId( bdm.dofParent( ldof ), 0 ), ldof );

        auto const attachment = bdm.dofAttachment( ldof );
        BOOST_REQUIRE( attachment.isValid() );
        BOOST_REQUIRE_GE( attachment.entityDim, 0 );
        BOOST_REQUIRE_LE( attachment.entityDim, 3 );
        BOOST_CHECK_EQUAL( attachment.kind, bdm.dofType( ldof ) );

        countByDim[attachment.entityDim] += 1;

        switch ( attachment.entityDim )
        {
        case 0:
            BOOST_CHECK_LT( attachment.entityId, BDMType::reference_convex_type::numVertices );
            BOOST_CHECK_LT( attachment.ordinal, runtimeDofPerVertex( bdm ) );
            break;
        case 1:
            BOOST_CHECK_LT( attachment.entityId, BDMType::reference_convex_type::numEdges );
            BOOST_CHECK_LT( attachment.ordinal, runtimeDofPerEdge( bdm ) );
            break;
        case 2:
            BOOST_CHECK_LT( attachment.entityId, BDMType::reference_convex_type::numFaces );
            BOOST_CHECK_LT( attachment.ordinal, runtimeDofPerFace( bdm ) );
            break;
        case 3:
            BOOST_CHECK_EQUAL( attachment.entityId, 0 );
            BOOST_CHECK_LT( attachment.ordinal, runtimeDofPerVolume( bdm ) );
            break;
        default:
            BOOST_FAIL( "invalid dof attachment entity dimension" );
        }
    }

    BOOST_CHECK_EQUAL( countByDim[0],
                       static_cast<uint16_type>( BDMType::reference_convex_type::numVertices * runtimeDofPerVertex( bdm ) ) );
    BOOST_CHECK_EQUAL( countByDim[1],
                       static_cast<uint16_type>( BDMType::reference_convex_type::numEdges * runtimeDofPerEdge( bdm ) ) );
    BOOST_CHECK_EQUAL( countByDim[2],
                       static_cast<uint16_type>( BDMType::reference_convex_type::numFaces * runtimeDofPerFace( bdm ) ) );
    BOOST_CHECK_EQUAL( countByDim[3], runtimeDofPerVolume( bdm ) );
}

template<int Order>
constexpr uint16_type bdmSimplex2dDof()
{
    return static_cast<uint16_type>( ( Order + 2 ) * ( Order + 3 ) );
}

template<int Order>
constexpr uint16_type bdmSimplex3dDof()
{
    return static_cast<uint16_type>( ( Order + 2 ) * ( Order + 3 ) * ( Order + 4 ) / 2 );
}

template<int Order>
void
checkBdmDynamicSimplex2d()
{
    using bdm_t = typename BrezziDouglasMarini<Dynamic>::template apply<2>::type;

    bdm_t bdm{ RuntimeOrder{ Order } };
    BOOST_CHECK_EQUAL( bdm.familyName(), "brezzidouglasmarini" );
    BOOST_CHECK_EQUAL( bdm.order(), Order );
    BOOST_CHECK_EQUAL( bdm.runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( bdm.internalOrder(), Order + 1 );
    BOOST_CHECK_EQUAL( bdm.localDofPerComponent(), bdmSimplex2dDof<Order>() );
    BOOST_CHECK_EQUAL( bdm.dofPerEdge(), Order + 2 );
    BOOST_CHECK_EQUAL( bdm.dofPerFace(), ( Order == 0 ) ? 0 : ( Order + 1 ) * ( Order + 1 ) - 1 );
    BOOST_CHECK_EQUAL( bdm.dofPerVolume(), 0 );

    checkFiniteEvaluationAtNodes( bdm );
    checkDofAttachmentContract( bdm );
    checkReferenceBdmInterpolationExactness<bdm_t, Order>( bdm );
}

template<int Order>
void
checkBdmDynamicSimplex3d()
{
    using bdm_t = typename BrezziDouglasMarini<Dynamic>::template apply<3>::type;

    bdm_t bdm{ RuntimeOrder{ Order } };
    BOOST_CHECK_EQUAL( bdm.familyName(), "brezzidouglasmarini" );
    BOOST_CHECK_EQUAL( bdm.order(), Order );
    BOOST_CHECK_EQUAL( bdm.runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( bdm.internalOrder(), Order + 1 );
    BOOST_CHECK_EQUAL( bdm.localDofPerComponent(), bdmSimplex3dDof<Order>() );
    BOOST_CHECK_EQUAL( bdm.dofPerEdge(), 0 );
    BOOST_CHECK_EQUAL( bdm.dofPerFace(), ( Order + 2 ) * ( Order + 3 ) / 2 );
    BOOST_CHECK_EQUAL( bdm.dofPerVolume(),
                       ( Order == 0 ) ? 0 : ( ( Order + 2 ) * ( Order + 3 ) * Order / 2 ) );

    checkFiniteEvaluationAtNodes( bdm );
    checkDofAttachmentContract( bdm );
    checkReferenceBdmInterpolationExactness<bdm_t, Order>( bdm );
}

template<typename BDMType>
void
checkBdm2dInteriorInterpolationOffset( BDMType const& bdm )
{
    IndexedVectorExpr2D expr;
    auto Ihloc = bdm.localInterpolant();
    bdm.interpolate( expr, Ihloc );

    const uint16_type firstInternalDof = static_cast<uint16_type>( BDMType::reference_convex_type::numEdges * runtimeDofPerEdge( bdm ) );
    const uint16_type nInternalDof = runtimeDofPerFace( bdm );
    BOOST_REQUIRE_GT( nInternalDof, 0 );

    double interiorNorm = 0.0;
    for ( uint16_type l = 0; l < nInternalDof; ++l )
        interiorNorm += std::abs( Ihloc( firstInternalDof + l ) );
    BOOST_CHECK_GT( interiorNorm, 1e-12 );
}
} // namespace

BOOST_AUTO_TEST_SUITE( bdm_static_suite )

BOOST_AUTO_TEST_CASE( bdm1_2d_compile_time_sanity )
{
    // BDM<Order=0> corresponds to BDM_k with k=1.
    using bdm1_2d_t = BrezziDouglasMarini<0>::apply<2>::type;

    static_assert( bdm1_2d_t::nDim == 2 );
    static_assert( bdm1_2d_t::nOrder == 1 );
    static_assert( bdm1_2d_t::nLocalDof == 6 );
    static_assert( bdm1_2d_t::nDofPerEdge == 2 );
    static_assert( bdm1_2d_t::nDofPerFace == 0 );
    static_assert( bdm1_2d_t::nDofPerVolume == 0 );

    bdm1_2d_t bdm1;
    BOOST_CHECK_EQUAL( bdm1.familyName(), "brezzidouglasmarini" );
    BOOST_CHECK_EQUAL( bdm1.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( bdm1.localDofPerComponent(), bdm1_2d_t::nLocalDof );

    checkFiniteEvaluationAtNodes( bdm1 );
    checkDofAttachmentContract( bdm1 );
    checkReferenceBdmInterpolationExactness<bdm1_2d_t, 0>( bdm1 );
}

BOOST_AUTO_TEST_CASE( bdm1_3d_compile_time_sanity )
{
    // BDM<Order=0> corresponds to BDM_k with k=1.
    using bdm1_3d_t = BrezziDouglasMarini<0>::apply<3>::type;

    static_assert( bdm1_3d_t::nDim == 3 );
    static_assert( bdm1_3d_t::nOrder == 1 );
    static_assert( bdm1_3d_t::nLocalDof == 12 );
    static_assert( bdm1_3d_t::nDofPerEdge == 0 );
    static_assert( bdm1_3d_t::nDofPerFace == 3 );
    static_assert( bdm1_3d_t::nDofPerVolume == 0 );

    bdm1_3d_t bdm1;
    BOOST_CHECK_EQUAL( bdm1.familyName(), "brezzidouglasmarini" );
    BOOST_CHECK_EQUAL( bdm1.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( bdm1.localDofPerComponent(), bdm1_3d_t::nLocalDof );

    checkFiniteEvaluationAtNodes( bdm1 );
    checkDofAttachmentContract( bdm1 );
    checkReferenceBdmInterpolationExactness<bdm1_3d_t, 0>( bdm1 );
}

BOOST_AUTO_TEST_CASE( bdm2_high_order_interior_dofs )
{
    // BDM<Order=1> corresponds to BDM_k with k=2.
    using bdm2_2d_t = BrezziDouglasMarini<1>::apply<2>::type;
    using bdm2_3d_t = BrezziDouglasMarini<1>::apply<3>::type;

    static_assert( bdm2_2d_t::nOrder == 2 );
    static_assert( bdm2_2d_t::nDofPerFace == 3 );
    static_assert( bdm2_2d_t::nDofPerVolume == 0 );

    static_assert( bdm2_3d_t::nOrder == 2 );
    static_assert( bdm2_3d_t::nDofPerFace == 6 );
    static_assert( bdm2_3d_t::nDofPerVolume == 6 );
    static_assert( bdm2_3d_t::nLocalDof == 30 );

    bdm2_2d_t bdm2_2d;
    bdm2_3d_t bdm2_3d;
    checkFiniteEvaluationAtNodes( bdm2_2d );
    checkFiniteEvaluationAtNodes( bdm2_3d );
    checkDofAttachmentContract( bdm2_2d );
    checkDofAttachmentContract( bdm2_3d );
    checkBdm2dInteriorInterpolationOffset( bdm2_2d );
    checkReferenceBdmInterpolationExactness<bdm2_2d_t, 1>( bdm2_2d );
    checkReferenceBdmInterpolationExactness<bdm2_3d_t, 1>( bdm2_3d );
}

BOOST_AUTO_TEST_CASE( bdm3_high_order_interior_dofs )
{
    // BDM<Order=2> corresponds to BDM_k with k=3.
    using bdm3_2d_t = BrezziDouglasMarini<2>::apply<2>::type;
    using bdm3_3d_t = BrezziDouglasMarini<2>::apply<3>::type;

    static_assert( bdm3_2d_t::nOrder == 3 );
    static_assert( bdm3_2d_t::nDofPerFace == 8 );
    static_assert( bdm3_2d_t::nDofPerVolume == 0 );

    static_assert( bdm3_3d_t::nOrder == 3 );
    static_assert( bdm3_3d_t::nDofPerFace == 10 );
    static_assert( bdm3_3d_t::nDofPerVolume == 20 );
    static_assert( bdm3_3d_t::nLocalDof == 60 );

    bdm3_2d_t bdm3_2d;
    bdm3_3d_t bdm3_3d;
    checkFiniteEvaluationAtNodes( bdm3_2d );
    checkFiniteEvaluationAtNodes( bdm3_3d );
    checkDofAttachmentContract( bdm3_2d );
    checkDofAttachmentContract( bdm3_3d );
    checkBdm2dInteriorInterpolationOffset( bdm3_2d );
    checkReferenceBdmInterpolationExactness<bdm3_2d_t, 2>( bdm3_2d );
    checkReferenceBdmInterpolationExactness<bdm3_3d_t, 2>( bdm3_3d );
}

BOOST_AUTO_TEST_CASE( bdm_dynamic_simplex_orders_sanity )
{
    checkBdmDynamicSimplex2d<0>();
    checkBdmDynamicSimplex2d<1>();
    checkBdmDynamicSimplex2d<2>();
    checkBdmDynamicSimplex2d<3>();

    checkBdmDynamicSimplex3d<0>();
    checkBdmDynamicSimplex3d<1>();
    checkBdmDynamicSimplex3d<2>();
    checkBdmDynamicSimplex3d<3>();

    using bdm_dynamic_2d_t = BrezziDouglasMarini<Dynamic>::apply<2>::type;
    bdm_dynamic_2d_t bdm2d{ RuntimeOrder{ 1 } };
    checkBdm2dInteriorInterpolationOffset( bdm2d );

    bdm_dynamic_2d_t bdm3dof{ RuntimeOrder{ 2 } };
    bdm_dynamic_2d_t bdm4dof{ RuntimeOrder{ 3 } };
    checkBdm2dInteriorInterpolationOffset( bdm3dof );
    checkBdm2dInteriorInterpolationOffset( bdm4dof );
    checkReferenceBdmInterpolationExactness<bdm_dynamic_2d_t, 1>( bdm2d );
    checkReferenceBdmInterpolationExactness<bdm_dynamic_2d_t, 2>( bdm3dof );
    checkReferenceBdmInterpolationExactness<bdm_dynamic_2d_t, 3>( bdm4dof );

    using bdm_dynamic_3d_t = BrezziDouglasMarini<Dynamic>::apply<3>::type;
    bdm_dynamic_3d_t bdm3d_order1{ RuntimeOrder{ 1 } };
    bdm_dynamic_3d_t bdm3d_order2{ RuntimeOrder{ 2 } };
    bdm_dynamic_3d_t bdm3d_order3{ RuntimeOrder{ 3 } };
    checkReferenceBdmInterpolationExactness<bdm_dynamic_3d_t, 1>( bdm3d_order1 );
    checkReferenceBdmInterpolationExactness<bdm_dynamic_3d_t, 2>( bdm3d_order2 );
    checkReferenceBdmInterpolationExactness<bdm_dynamic_3d_t, 3>( bdm3d_order3 );
}

BOOST_AUTO_TEST_SUITE_END()
