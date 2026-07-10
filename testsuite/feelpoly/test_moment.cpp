/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2011 Université Joseph Fourier (Grenoble I)
    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
   \file test_moment.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-12-11
 */
#define BOOST_TEST_MODULE test_moment
#include <feel/feelcore/testsuite.hpp>

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <feel/feelpoly/moment.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( moment )

typedef boost::mpl::list<boost::mpl::int_<2> > test_types;

template<typename FE>
void
checkCrouzeixRaviartSimplexElement( typename FE::value_type tol )
{
    using namespace Feel;

    FE fe;
    typename FE::reference_convex_type ref;

    static_assert( FE::nDim == 2 || FE::nDim == 3 );
    BOOST_CHECK_EQUAL( FE::nOrder, 1 );
    BOOST_CHECK_EQUAL( FE::nLocalDof, FE::nDim + 1 );
    BOOST_CHECK_EQUAL( fe.nbPoints(), FE::nDim + 1 );
    BOOST_CHECK_EQUAL( FE::nDofPerVertex, 0 );
    BOOST_CHECK_EQUAL( FE::nDofPerVolume, 0 );

    if constexpr ( FE::nDim == 2 )
    {
        BOOST_CHECK_EQUAL( FE::nDofPerEdge, 1 );
        BOOST_CHECK_EQUAL( FE::nDofPerFace, 0 );
    }
    else
    {
        BOOST_CHECK_EQUAL( FE::nDofPerEdge, 0 );
        BOOST_CHECK_EQUAL( FE::nDofPerFace, 1 );
    }

    for ( uint16_type f = 0; f < FE::nDim + 1; ++f )
    {
        auto const& ptsOnFace = fe.points( f );
        BOOST_REQUIRE_EQUAL( ptsOnFace.size2(), 1 );
        for ( uint16_type c = 0; c < FE::nDim; ++c )
        {
            BOOST_CHECK_SMALL( std::abs( ptsOnFace( c, 0 ) - ref.faceBarycenter( f )( c ) ), tol );
            BOOST_CHECK_SMALL( std::abs( fe.points()( c, f ) - ref.faceBarycenter( f )( c ) ), tol );
        }
    }

    auto evalAtDofs = FE::polyset_type::toMatrix( fe.evaluate( fe.points() ) );
    auto dofIdentity = ublas::identity_matrix<typename FE::value_type>( evalAtDofs.size1() );
    BOOST_CHECK_EQUAL( evalAtDofs.size1(), FE::nDim + 1 );
    BOOST_CHECK_EQUAL( evalAtDofs.size2(), FE::nDim + 1 );
    typename FE::value_type dofIdentityError = ublas::norm_frobenius( evalAtDofs - dofIdentity );
    BOOST_CHECK_SMALL( dofIdentityError, tol );

    typename FE::points_type samplePts( FE::nDim, FE::nDim + 2 );
    ublas::subrange( samplePts, 0, FE::nDim, 0, FE::nDim + 1 ) = ref.vertices();
    ublas::column( samplePts, FE::nDim + 1 ) = ref.barycenter();

    auto const linearFunction = []( auto const& pts, uint16_type q )
        {
            typename FE::value_type value = 1.25;
            for ( uint16_type c = 0; c < FE::nDim; ++c )
                value += ( ( c % 2 ) ? -1.0 : 1.0 ) * ( c + 2.0 ) * pts( c, q );
            return value;
        };

    std::vector<typename FE::value_type> dofValues( FE::nLocalDof );
    for ( uint16_type i = 0; i < FE::nLocalDof; ++i )
        dofValues[i] = linearFunction( fe.points(), i );

    auto basisAtSamples = FE::polyset_type::toMatrix( fe.evaluate( samplePts ) );
    BOOST_CHECK_EQUAL( basisAtSamples.size1(), FE::nLocalDof );
    BOOST_CHECK_EQUAL( basisAtSamples.size2(), FE::nDim + 2 );

    for ( uint16_type q = 0; q < samplePts.size2(); ++q )
    {
        typename FE::value_type interpolated = 0;
        for ( uint16_type i = 0; i < FE::nLocalDof; ++i )
            interpolated += dofValues[i] * basisAtSamples( i, q );
        BOOST_CHECK_SMALL( std::abs( interpolated - linearFunction( samplePts, q ) ), tol );
    }

    ublas::vector<typename FE::matrix_type> der = fe.derivate( fe.points() );
    for ( uint16_type c = 0; c < FE::nDim; ++c )
    {
        typename FE::value_type derivativeError = ublas::norm_frobenius( der[c] - fe.derivate( c ).evaluate( fe.points() ) );
        BOOST_CHECK_SMALL( derivativeError, tol );
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( test_QK, T, test_types )
{
    using namespace Feel;

    Moment<2,T::value,Hypercube<2> > m;
#if 0
    std::cout << "1  :" << m.template pick<Scalar>( 0 ).evaluate( m.points() ) << "\n";
    std::cout << "x  :" << m.template pick<Scalar>( 1 ).evaluate( m.points() ) << "\n";
    std::cout << "x^2:" << m.template pick<Scalar>( 2 ).evaluate( m.points() ) << "\n";
    std::cout << "y  :" << m.template pick<Scalar>( 3 ).evaluate( m.points() ) << "\n";
    std::cout << "y^2:" << m.template pick<Scalar>( 6 ).evaluate( m.points() ) << "\n";

    auto p  = m.template pick<Scalar>( 2 ) - m.template pick<Scalar>( 6 );
    std::cout << "x^2-y^2:" << p.evaluate( m.points() ) << "\n";

    std::cout << "d 1/dx  :" << m.template pick<Scalar>( 0 ).derivate( 0,m.points() ) << "\n";
    std::cout << "d x/dx  :" << m.template pick<Scalar>( 1 ).derivate( 0,m.points() ) << "\n";
    std::cout << "d x^2/dx:" << m.template pick<Scalar>( 2 ).derivate( 0,m.points() ) << "\n";
    std::cout << "d y/dx   :" << m.template pick<Scalar>( 3 ).derivate( 0,m.points() ) << "\n";
    std::cout << "d y^2 /dx:" << m.template pick<Scalar>( 6 ).derivate( 0,m.points() ) << "\n";
    std::cout << "dx^2-y^2/x :" << p.derivate( 0, m.points() ) << "\n";
    std::cout << "\n";
    std::cout << "d 1/dy  :" << m.template pick<Scalar>( 0 ).derivate( 1,m.points() ) << "\n";
    std::cout << "d x/dy  :" << m.template pick<Scalar>( 1 ).derivate( 1,m.points() ) << "\n";
    std::cout << "d x^2/dy:" << m.template pick<Scalar>( 2 ).derivate( 1,m.points() ) << "\n";
    std::cout << "d y/dy   :" << m.template pick<Scalar>( 3 ).derivate( 1,m.points() ) << "\n";
    std::cout << "d y^2 /dy:" << m.template pick<Scalar>( 6 ).derivate( 1,m.points() ) << "\n";
    std::cout << "d x^2-y^2 /dy:" << p.derivate( 1,m.points() ) << "\n";

    PolynomialSet<Moment<2,T::value,Hypercube<2> >, Scalar> pset( m );
    pset.insert( m.template pick<Scalar>( 0 ).toSet(), true );
    pset.insert( m.template pick<Scalar>( 1 ).toSet() );
    pset.insert( m.template pick<Scalar>( 3 ).toSet() );
    pset.insert( p.toSet() );

    std::cout << "pset :" << pset.evaluate( m.points() ) << "\n";
    std::cout << "d pset/dx :" << pset.derivate( 0, m.points() ) << "\n";
    std::cout << "d pset/dy :" << pset.derivate( 1, m.points() ) << "\n";

    fem::detail::RannacherTurekPolynomialSet<2,Scalar> RQ;
    std::cout << "rq :" << RQ.evaluate( m.points() ) << "\n";
    std::cout << "d rq/dx :" << RQ.derivate( 0, m.points() ) << "\n";
    std::cout << "d rq/dy :" << RQ.derivate( 1, m.points() ) << "\n";

    fem::CrouzeixRaviart<2,2,Scalar,double,Hypercube> cr;
    std::cout << "cr :" << cr.evaluate( m.points() ) << "\n";

    fem::detail::RannacherTurekPolynomialSet<2,Vectorial> RQv;
    std::cout << "rqv :" << RQv.evaluate( m.points() ) << "\n";
    std::cout << "d rqv/dx :" << RQv.derivate( 0, m.points() ) << "\n";
    std::cout << "d rqv/dy :" << RQv.derivate( 1, m.points() ) << "\n";
#endif
    fem::CrouzeixRaviart<2,2,Scalar,double,Hypercube> cr;
    auto crEval = cr.evaluate( m.points() );
    BOOST_CHECK_EQUAL( crEval.size1(), cr.nLocalDof );
    BOOST_CHECK_EQUAL( crEval.size2(), m.points().size2() );

    fem::CrouzeixRaviart<2,2,Vectorial,double,Hypercube> crv;
    auto crvEval = crv.evaluate( m.points() );
    BOOST_CHECK( crvEval.size1() > 0 );
    BOOST_CHECK_EQUAL( crvEval.size2(), m.points().size2() );

}

BOOST_AUTO_TEST_CASE( crouzeix_raviart_simplex_2d_dofs_are_edge_midpoints )
{
    using namespace Feel;
    using fe_type = fem::CrouzeixRaviart<2, 2, Scalar, double, Simplex>;
    checkCrouzeixRaviartSimplexElement<fe_type>( 1e-12 );
}

BOOST_AUTO_TEST_CASE( crouzeix_raviart_simplex_3d_dofs_are_face_barycenters )
{
    using namespace Feel;
    using fe_type = fem::CrouzeixRaviart<3, 3, Scalar, double, Simplex>;
    checkCrouzeixRaviartSimplexElement<fe_type>( 1e-12 );
}

BOOST_AUTO_TEST_SUITE_END()
