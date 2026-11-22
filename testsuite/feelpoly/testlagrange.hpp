/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2006 EPFL
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#ifndef FEELPP_TEST_LAGRANGE_HPP
#define FEELPP_TEST_LAGRANGE_HPP

#include <feel/feelpoly/lagrange.hpp>

/** @brief Shared nodal-identity and differentiation assertions for Lagrange bases.
 * @tparam FE Concrete static finite element selected by a Boost.Test dataset.
 */
template<typename FE>
class TestLagrange
{
public:
    /** Finite-element scalar type used for residuals and tolerances. */
    using value_type = typename FE::value_type;
    /** Owning matrix type for sampled and coefficient-based derivatives. */
    using matrix_type = typename FE::matrix_type;

    /** @brief Set the absolute Frobenius-norm tolerance for this configuration.
     * @param tolerance Maximum accepted nodal and derivative residual.
     */
    explicit TestLagrange( value_type tolerance ) : M_tolerance( tolerance ) {}

    /** @brief Check nodal identity and sampled/coefficient derivative agreement. */
    void operator()() const
    {
        FE fe;
        matrix_type const values( FE::polyset_type::toMatrix( fe.evaluate( fe.points() ) ) );
        BOOST_REQUIRE_EQUAL( values.size1(), values.size2() );
        BOOST_REQUIRE_GT( values.size1(), 0 );
        value_type const identityError = Feel::ublas::norm_frobenius(
            values - Feel::ublas::identity_matrix<value_type>( values.size1() ) );
        BOOST_CHECK_SMALL( identityError, M_tolerance );

        auto const derivatives = fe.derivate( fe.points() );
        BOOST_REQUIRE_EQUAL( derivatives.size(), FE::nDim );
        for ( unsigned direction = 0; direction < FE::nDim; ++direction )
        {
            matrix_type const fromCoefficients = fe.derivate( direction ).evaluate( fe.points() );
            BOOST_REQUIRE_EQUAL( derivatives[direction].size1(), fromCoefficients.size1() );
            BOOST_REQUIRE_EQUAL( derivatives[direction].size2(), fromCoefficients.size2() );
            value_type const derivativeError = Feel::ublas::norm_frobenius(
                derivatives[direction] - fromCoefficients );
            BOOST_TEST_CONTEXT( "derivative direction=" << direction )
                BOOST_CHECK_SMALL( derivativeError, M_tolerance );
        }
    }

private:
    /** Absolute tolerance selected according to the dataset's polynomial order. */
    value_type M_tolerance;
};
#endif
