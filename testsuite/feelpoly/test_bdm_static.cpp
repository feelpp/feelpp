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

template<typename BDMType>
void
checkDofAttachmentContract( BDMType const& bdm )
{
    static_assert( !BDMType::is_product );

    std::array<uint16_type, 4> countByDim = { 0, 0, 0, 0 };

    for ( uint16_type ldof = 0; ldof < BDMType::nLocalDof; ++ldof )
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
            BOOST_CHECK_LT( attachment.ordinal, BDMType::nDofPerVertex );
            break;
        case 1:
            BOOST_CHECK_LT( attachment.entityId, BDMType::reference_convex_type::numEdges );
            BOOST_CHECK_LT( attachment.ordinal, BDMType::nDofPerEdge );
            break;
        case 2:
            BOOST_CHECK_LT( attachment.entityId, BDMType::reference_convex_type::numFaces );
            BOOST_CHECK_LT( attachment.ordinal, BDMType::nDofPerFace );
            break;
        case 3:
            BOOST_CHECK_EQUAL( attachment.entityId, 0 );
            BOOST_CHECK_LT( attachment.ordinal, BDMType::nDofPerVolume );
            break;
        default:
            BOOST_FAIL( "invalid dof attachment entity dimension" );
        }
    }

    BOOST_CHECK_EQUAL( countByDim[0],
                       static_cast<uint16_type>( BDMType::reference_convex_type::numVertices * BDMType::nDofPerVertex ) );
    BOOST_CHECK_EQUAL( countByDim[1],
                       static_cast<uint16_type>( BDMType::reference_convex_type::numEdges * BDMType::nDofPerEdge ) );
    BOOST_CHECK_EQUAL( countByDim[2],
                       static_cast<uint16_type>( BDMType::reference_convex_type::numFaces * BDMType::nDofPerFace ) );
    BOOST_CHECK_EQUAL( countByDim[3], BDMType::nDofPerVolume );
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
}

BOOST_AUTO_TEST_SUITE_END()
