/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-01-04

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
   \file test_traits.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-04
   \brief Tests for feelpoly traits and concepts
 */

#define BOOST_TEST_MODULE test_feelpoly_traits
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/policy.hpp>
#include <feel/feelpoly/traits.hpp>

#include <type_traits>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
struct OrderPoly
{
    static constexpr int order() { return 2; }
    static constexpr bool is_linear = false;
};

struct LinearPoly
{
    static constexpr int order() { return 1; }
    static constexpr bool is_linear = true;
};

struct NoOrder
{
    static constexpr bool is_linear = false;
};
} // namespace

BOOST_AUTO_TEST_SUITE( traits_suite )

BOOST_AUTO_TEST_CASE( test_polynomial_type_traits )
{
    static_assert( is_scalar_polynomial<Scalar<2>>::value, "Scalar should be scalar polynomial" );
    static_assert( is_vector_polynomial<Vectorial<3>>::value, "Vectorial should be vector polynomial" );
    static_assert( is_tensor2_polynomial<Tensor2<2>>::value, "Tensor2 should be tensor2 polynomial" );
    static_assert( !is_scalar_polynomial<Vectorial<2>>::value, "Vectorial should not be scalar polynomial" );

    static_assert( ScalarPolynomialType<Scalar<2>>, "ScalarPolynomialType should accept Scalar" );
    static_assert( VectorPolynomialType<Vectorial<2>>, "VectorPolynomialType should accept Vectorial" );
    static_assert( Tensor2PolynomialType<Tensor2<3>>, "Tensor2PolynomialType should accept Tensor2" );
    static_assert( !ScalarPolynomialType<Vectorial<2>>, "ScalarPolynomialType should reject Vectorial" );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_polynomial_order_traits )
{
    static_assert( polynomial_order<OrderPoly>::value == 2, "polynomial_order should use order()" );
    static_assert( polynomial_order_v<OrderPoly>, "polynomial_order_v should be true for non-zero order" );
    static_assert( is_linear_polynomial_v<LinearPoly>, "is_linear_polynomial should use is_linear" );
    static_assert( !is_linear_polynomial_v<OrderPoly>, "non-linear polynomial should be false" );

    static_assert( LinearPolynomialType<LinearPoly>, "LinearPolynomialType should accept linear type" );
    static_assert( !LinearPolynomialType<OrderPoly>, "LinearPolynomialType should reject non-linear type" );
    static_assert( PolynomialOrderable<OrderPoly>, "PolynomialOrderable should accept order() types" );
    static_assert( !PolynomialOrderable<NoOrder>, "PolynomialOrderable should reject missing order()" );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
