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
 * @brief Compile-time/runtime sanity checks for RT0 only
 */
#define BOOST_TEST_MODULE test_rt0_static
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/raviartthomas.hpp>

#include <array>
#include <cmath>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
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

template<typename RTType>
void
checkDofAttachmentContract( RTType const& rt )
{
    static_assert( !RTType::is_product );

    std::array<uint16_type, 4> countByDim = { 0, 0, 0, 0 };

    for ( uint16_type ldof = 0; ldof < RTType::nLocalDof; ++ldof )
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
            BOOST_CHECK_LT( attachment.ordinal, RTType::nDofPerVertex );
            break;
        case 1:
            BOOST_CHECK_LT( attachment.entityId, RTType::reference_convex_type::numEdges );
            BOOST_CHECK_LT( attachment.ordinal, RTType::nDofPerEdge );
            break;
        case 2:
            BOOST_CHECK_LT( attachment.entityId, RTType::reference_convex_type::numFaces );
            BOOST_CHECK_LT( attachment.ordinal, RTType::nDofPerFace );
            break;
        case 3:
            BOOST_CHECK_EQUAL( attachment.entityId, 0 );
            BOOST_CHECK_LT( attachment.ordinal, RTType::nDofPerVolume );
            break;
        default:
            BOOST_FAIL( "invalid dof attachment entity dimension" );
        }
    }

    BOOST_CHECK_EQUAL( countByDim[0],
                       static_cast<uint16_type>( RTType::reference_convex_type::numVertices * RTType::nDofPerVertex ) );
    BOOST_CHECK_EQUAL( countByDim[1],
                       static_cast<uint16_type>( RTType::reference_convex_type::numEdges * RTType::nDofPerEdge ) );
    BOOST_CHECK_EQUAL( countByDim[2],
                       static_cast<uint16_type>( RTType::reference_convex_type::numFaces * RTType::nDofPerFace ) );
    BOOST_CHECK_EQUAL( countByDim[3], RTType::nDofPerVolume );
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

BOOST_AUTO_TEST_SUITE_END()
