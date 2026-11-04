// SPDX-FileCopyrightText: 2005-2006 EPFL
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
     \file test_poly.cpp
     \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
     \date 2005-08-17
 */

#define BOOST_TEST_MODULE test_poly
#include <boost/test/data/test_case.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpoly/polynomial.hpp>
#include <feel/feelpoly/operations.hpp>
#include <cmath>
#include <array>

namespace bdata = boost::unit_test::data;
namespace ublas = boost::numeric::ublas;
using namespace Feel;
/** @brief Weighted quadratic with a mixed term in two dimensions. */
struct QuadraticPolynomial
{
    /** Scalar value type used by polynomial projection. */
    typedef double value_type;
    /** @brief Evaluate the quadratic at each point column. */
    ublas::vector<double>
    operator()( ublas::matrix<double,ublas::column_major> const& m ) const
    {
        ublas::vector<double> v( m.size2() );
        v.clear();

        for ( int i = 0; i < m.size2(); ++i )
        {
            for ( int j = 0; j < m.size1(); ++j )
                v( i ) += ( j+2 )*m( j,i )*m( j,i );

            if ( m.size1() == 2 )
                v( i ) += 100*m( 0,i )*m( 1,i )+200*m( 1,i )*m( 0,i );

            //v( i ) = m(0,i)*m(0,i);
        }

        return v;
    }
};

/** @brief Asymmetric probes on the reference cell and its boundary. */
template<int Dim>
ublas::matrix<double,ublas::column_major> makePoints();

template<>
ublas::matrix<double,ublas::column_major> makePoints<1>()
{
    ublas::matrix<double,ublas::column_major> pts( 1, 3 );
    pts( 0,0 ) = 0;
    pts( 0,1 ) = -1;
    pts( 0,2 ) = 1;
    return pts;
}

template<>
ublas::matrix<double,ublas::column_major> makePoints<2>()
{
    ublas::matrix<double,ublas::column_major> pts( 2, 4 );
    pts( 0,0 ) = -1; pts( 1,0 ) = -1;
    pts( 0,1 ) =  1; pts( 1,1 ) = -1;
    pts( 0,2 ) = -1; pts( 1,2 ) =  1;
    pts( 0,3 ) =  0; pts( 1,3 ) =  0;
    return pts;
}

template<>
ublas::matrix<double,ublas::column_major> makePoints<3>()
{
    ublas::matrix<double,ublas::column_major> pts( 3, 4 );
    pts( 0,0 ) = -1; pts( 1,0 ) = -1; pts( 2,0 ) = -1;
    pts( 0,1 ) =  1; pts( 1,1 ) = -1; pts( 2,1 ) = -1;
    pts( 0,2 ) = -1; pts( 1,2 ) =  1; pts( 2,2 ) = -1;
    pts( 0,3 ) = -0.5; pts( 1,3 ) = -0.5; pts( 2,3 ) = -0.5;
    return pts;
}

/** @brief Check projection, gradients and Hessians against analytic derivatives. */
template<int Dim, int QuadOrder, int PolyOrder>
void checkPolynomialProjection()
{
    using polyset_type = Feel::detail::OrthonormalPolynomialSet<Dim,PolyOrder,Dim,Scalar>;

    auto pts = makePoints<Dim>();
    polyset_type ps;
    Polynomial<polyset_type, Scalar> p1 = project( ps, QuadraticPolynomial(), IM<Dim,QuadOrder>() );
    PolynomialSet<polyset_type, Scalar> p( ps );
    p.insert( p1.toSet( true ) );
    p.insert( p1.toSet( true ) );

    const auto values = p1.evaluate(pts);
    const auto expected = QuadraticPolynomial()(pts);
    const auto grad_vals = p.gradient().evaluate( pts );
    const auto hess_vals = p.gradient().gradient().evaluate( pts );

    BOOST_TEST_CONTEXT( "dim=" << Dim << " quad=" << QuadOrder )
    {
        BOOST_CHECK_EQUAL( grad_vals.size1(), p.polynomialDimension()*Dim );
        BOOST_CHECK_EQUAL( hess_vals.size1(), p.polynomialDimension()*Dim*Dim );
        BOOST_CHECK_EQUAL( grad_vals.size2(), pts.size2() );
        BOOST_CHECK_EQUAL( hess_vals.size2(), pts.size2() );
        auto const last = p.polynomialDimension()-1;
        for ( std::size_t q=0; q<pts.size2(); ++q )
        {
            BOOST_CHECK_SMALL(values(0,q)-expected(q),2e-9);
            for ( unsigned i=0; i<Dim; ++i )
            {
                auto expectedGradient=2*(i+2)*pts(i,q);
                if constexpr (Dim==2) expectedGradient+=300*pts(1-i,q);
                BOOST_CHECK_SMALL(grad_vals(last*Dim+i,q)-expectedGradient,2e-9);
                for(unsigned j=0;j<Dim;++j)
                {
                    auto expectedHessian=i==j ? 2.0*(i+2) : (Dim==2 ? 300.0 : 0.0);
                    BOOST_CHECK_SMALL(hess_vals(last*Dim*Dim+i*Dim+j,q)-expectedHessian,2e-9);
                }
            }
        }
    }
}

/** Shared signature of the dimension/order configurations. */
using Runner = void(*)();
/** Projection configurations registered through Boost.Test data cases. */
static const std::array<Runner,3> polynomialCases = {
    &checkPolynomialProjection<1,4,2>,
    &checkPolynomialProjection<2,4,2>,
    &checkPolynomialProjection<3,6,3>
};

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( polynomial_suite )

BOOST_DATA_TEST_CASE( polynomialVariants, bdata::make( polynomialCases ), runCase )
{
    runCase();
}

BOOST_AUTO_TEST_SUITE_END()
