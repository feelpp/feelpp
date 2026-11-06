/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#pragma once
#include "natural_second_operator_checks.hpp"
#include <feel/feeldiscr/bdmh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/ned2h.hpp>
#include <feel/feeldiscr/pchv.hpp>

namespace Feel::natural_second_checks
{
/** All Hessian components versus physical differences of the gradient API. */
template <typename SpacePtr>
void checkVectorHessian( SpacePtr const& space, bool polynomial )
{
    using element_type = typename SpacePtr::element_type::element_type;
    constexpr int Dim = element_type::nDim;
    auto mesh = space->mesh();
    auto field = space->element();
    field.on( _range = elements( mesh ), _expr = quadraticField<Dim>(), _geomap = GeomapStrategyType::GEOMAP_HO );
    auto points = localPoints( mesh );
    typename element_type::matrix_node_type empty;
    typename element_type::hess_array_type hessians( boost::extents[points.size2()] );
    typename element_type::id_array_type laplacians( boost::extents[points.size2()] );
    typename element_type::id_array_type gradDiv( boost::extents[points.size2()] );
    typename element_type::id_array_type curlCurl( boost::extents[points.size2()] );
    for ( auto& value : hessians )
        value.setConstant( std::numeric_limits<double>::quiet_NaN() );
    field.hessInterpolate( points, hessians, false, empty );
    field.laplacianInterpolate( points, laplacians, false, empty );
    field.graddivInterpolate( points, gradDiv, false, empty );
    field.curlcurlInterpolate( points, curlCurl, false, empty );
    BOOST_TEST_CONTEXT( "full Hessian " << space->basisName() << " FE" << space->order() << " G" << mesh->order() )
    {
        constexpr double step = 1e-3;
        for ( int b = 0; b < Dim; ++b )
        {
            std::array<typename element_type::grad_array_type, 4> sampled;
            constexpr double offsets[4] = { 1., -1., .5, -.5 };
            for ( int sample = 0; sample < 4; ++sample )
            {
                auto shifted = points;
                for ( size_type q = 0; q < points.size2(); ++q )
                    shifted( b, q ) += step * offsets[sample];
                sampled[sample].resize( boost::extents[points.size2()] );
                field.gradInterpolate( shifted, sampled[sample], false, empty );
            }
            for ( size_type q = 0; q < points.size2(); ++q )
                for ( int c = 0; c < Dim; ++c )
                    for ( int a = 0; a < Dim; ++a )
                    {
                        double const coarse = ( sampled[0][q]( c, a ) - sampled[1][q]( c, a ) ) / ( 2 * step );
                        double const fine = ( sampled[2][q]( c, a ) - sampled[3][q]( c, a ) ) / step;
                        double expected = ( 4 * fine - coarse ) / 3;
                        if ( polynomial )
                        {
                            if ( c == 0 )
                                expected = 2 * ( a == 0 && b == 0 ) + ( a == 0 && b == 1 ) + ( a == 1 && b == 0 );
                            else if ( c == 1 )
                                expected = 2 * ( a == 1 && b == 1 ) + 2 * ( a == 0 && b == 1 ) + 2 * ( a == 1 && b == 0 );
                            else
                                expected = 2 * ( a == 2 && b == 2 ) + ( a == 0 && b == 2 ) + ( a == 2 && b == 0 ) + ( a == 1 && b == 2 ) + ( a == 2 && b == 1 );
                        }
                        BOOST_CHECK_SMALL( hessians[q]( c + Dim * a, b ) - expected, ( polynomial ? 2e-8 : 2e-5 ) * ( 1 + std::abs( expected ) ) );
                        BOOST_CHECK_SMALL( hessians[q]( c + Dim * a, b ) - hessians[q]( c + Dim * b, a ), 2e-9 * ( 1 + std::abs( expected ) ) );
                    }
        }
        for ( size_type q = 0; q < points.size2(); ++q )
            for ( int c = 0; c < Dim; ++c )
            {
                double trace = 0;
                for ( int a = 0; a < Dim; ++a )
                    trace += hessians[q]( c + Dim * a, a );
                BOOST_CHECK_SMALL( laplacians[q]( c, 0 ) - trace, 2e-8 * ( 1 + std::abs( trace ) ) );
                BOOST_CHECK_SMALL( gradDiv[q]( c, 0 ) - curlCurl[q]( c, 0 ) - trace, 2e-8 * ( 1 + std::abs( trace ) ) );
            }
        auto residual = graddivv( field ) - curlcurlv( field ) - laplacianv( field );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ), _expr = residual, _quad = 8,
                                   _geomap = GeomapStrategyType::GEOMAP_HO ),
                           2e-8 );
        // Exercise the full flattened VF shape on both cell and face contexts.
        double const energy = integrate( _range = elements( mesh ), _expr = inner( hessv( field ), hessv( field ) ),
                                         _quad = 8, _geomap = GeomapStrategyType::GEOMAP_HO )
                                  .evaluate()( 0, 0 );
        double const boundaryEnergy = integrate( _range = boundaryfaces( mesh ), _expr = inner( hessv( field ), hessv( field ) ),
                                                 _quad = 8, _geomap = GeomapStrategyType::GEOMAP_HO )
                                          .evaluate()( 0, 0 );
        BOOST_CHECK( std::isfinite( energy ) && energy >= 0 );
        BOOST_CHECK( std::isfinite( boundaryEnergy ) && boundaryEnergy >= 0 );
        if ( polynomial ) BOOST_CHECK_SMALL( energy - ( Dim == 2 ? 18. : 26. ), 2e-8 );
    }
}

/** Excite the highest-degree Piola modes and let VF choose quadrature. */
template <typename SpacePtr>
void checkAutomaticPiolaQuadrature( SpacePtr const& space )
{
    using element_type = typename SpacePtr::element_type::element_type;
    constexpr int Dim = element_type::nDim;
    auto mesh = space->mesh();
    auto field = space->element();
    if constexpr ( element_type::is_hdiv_conforming )
    {
        if constexpr ( Dim == 2 )
            field.on( _range = elements( mesh ), _expr = vec( Px() * Px() * Px(), Px() * Px() * Py() ) );
        else
            field.on( _range = elements( mesh ), _expr = vec( Px() * Px() * Px(), Px() * Px() * Py(), Px() * Px() * Pz() ) );
    }
    else
    {
        if constexpr ( Dim == 2 )
            field.on( _range = elements( mesh ), _expr = vec( -Px() * Px() * Py(), Px() * Px() * Px() ) );
        else
            field.on( _range = elements( mesh ), _expr = vec( -Px() * Px() * Py(), Px() * Px() * Px(), cst( 0. ) ) );
    }
    BOOST_TEST_CONTEXT( "automatic quadrature " << space->basisName() << " FE" << space->order() )
    {
        double const exactMass = 1. / 7. + ( element_type::is_hdiv_conforming ? Dim - 1 : 1 ) / 15.;
        double const mass = normL2( _range = elements( mesh ), _expr = idv( field ) );
        BOOST_CHECK_SMALL( mass * mass - exactMass, 2e-9 );
        if constexpr ( element_type::is_hdiv_conforming )
        {
            double const divergence = normL2( _range = elements( mesh ), _expr = divv( field ) );
            double const gradDiv = normL2( _range = elements( mesh ), _expr = graddivv( field ) );
            BOOST_CHECK_SMALL( divergence * divergence - ( Dim + 2. ) * ( Dim + 2. ) / 5., 2e-9 );
            BOOST_CHECK_SMALL( gradDiv * gradDiv - 4. * ( Dim + 2. ) * ( Dim + 2. ) / 3., 2e-9 );
        }
        else
        {
            double const curl = normL2( _range = elements( mesh ), _expr = curlv( field ) );
            double const curlCurl = normL2( _range = elements( mesh ), _expr = curlcurlv( field ) );
            BOOST_CHECK_SMALL( curl * curl - 16. / 5., 2e-9 );
            BOOST_CHECK_SMALL( curlCurl * curlCurl - 64. / 3., 2e-9 );
        }
    }
}

template <bool Assemble = true, typename SpacePtr>
void checkPiolaPolynomial( SpacePtr const& space )
{
    checkPolynomial<Assemble>( space );
    checkVectorHessian( space, true );
    if ( space->mesh()->order() == 1 ) checkAutomaticPiolaQuadrature( space );
    if constexpr ( Assemble )
    {
        constexpr int Dim = SpacePtr::element_type::nDim;
        auto field = space->element();
        field.on( _range = elements( space->mesh() ), _expr = quadraticField<Dim>(), _geomap = GeomapStrategyType::GEOMAP_HO );
        auto trial = space->element(), test = space->element();
        auto a = form2( _trial = space, _test = space );
        a = integrate( _range = elements( space->mesh() ),
                       _expr = inner( hesst( trial ), hess( test ) ) + inner( laplaciant( trial ), laplacian( test ) ),
                       _quad = 8, _geomap = GeomapStrategyType::GEOMAP_HO );
        BOOST_CHECK_SMALL( a( field, field ) - ( Dim == 2 ? 26. : 38. ), 2e-8 );
    }
}

template <int Dim, int GeometryOrder>
auto staticCurved()
{
    return createGMSHMesh(
        _mesh = new Mesh<Simplex<Dim, GeometryOrder>>,
        _desc = domain( _name = "piola-hessian-static-curved-" + std::to_string( Dim ) + "-g" + std::to_string( GeometryOrder ),
                        _shape = "ellipsoid", _dim = Dim, _xmin = -1., _ymin = -1., _zmin = -1., _h = .9 ) );
}

template <int Dim>
void checkPiolaDimension()
{
    for ( uint16_type geometryOrder : { 1, 2, 3 } )
    {
        auto mesh = box<Dim>( geometryOrder );
        checkVectorHessian( Pchv<Dynamic>( mesh, RuntimeOrder{ 0 } ), false );
        checkVectorHessian( Pchv<Dynamic>( mesh, RuntimeOrder{ 1 } ), false );
        checkVectorHessian( Pchv<Dynamic>( mesh, RuntimeOrder{ 2 } ), true );
        checkPiolaPolynomial<false>( RTh<2>( mesh ) );
        checkPiolaPolynomial<false>( BDMh<2>( mesh ) );
        checkPiolaPolynomial<false>( Ned1h<2>( mesh ) );
        checkPiolaPolynomial<false>( Ned2h<2>( mesh ) );
        checkPiolaPolynomial( RTh<Dynamic>( mesh, RuntimeOrder{ 2 } ) );
        checkPiolaPolynomial( BDMh<Dynamic>( mesh, RuntimeOrder{ 2 } ) );
        checkPiolaPolynomial( Ned1h<Dynamic>( mesh, RuntimeOrder{ 2 } ) );
        checkPiolaPolynomial( Ned2h<Dynamic>( mesh, RuntimeOrder{ 2 } ) );
    }
    for ( uint16_type geometryOrder : { 2, 3 } )
    {
        auto mesh = curved<Dim>( geometryOrder );
        for ( uint16_type order : { 0, 1 } )
        {
            auto check = [&]( auto const& space )
            {
                checkCurved( space );
                checkVectorHessian( space, false );
            };
            check( RTh<Dynamic>( mesh, RuntimeOrder{ order } ) );
            check( BDMh<Dynamic>( mesh, RuntimeOrder{ order } ) );
            check( Ned1h<Dynamic>( mesh, RuntimeOrder{ order } ) );
            check( Ned2h<Dynamic>( mesh, RuntimeOrder{ order } ) );
            check( Pchv<Dynamic>( mesh, RuntimeOrder{ order } ) );
        }
    }
    auto mesh = staticBox<Dim>();
    checkVectorHessian( Pchv<2>( mesh ), true );
    checkPiolaPolynomial<false>( RTh<2>( mesh ) );
    checkPiolaPolynomial<false>( Ned1h<2>( mesh ) );
    auto staticG2 = staticCurved<Dim, 2>();
    checkVectorHessian( RTh<1>( staticG2 ), false );
    checkVectorHessian( Ned1h<Dynamic>( staticG2, RuntimeOrder{ 1 } ), false );
    auto staticG3 = staticCurved<Dim, 3>();
    checkVectorHessian( RTh<Dynamic>( staticG3, RuntimeOrder{ 1 } ), false );
    checkVectorHessian( Ned1h<1>( staticG3 ), false );
}
} // namespace Feel::natural_second_checks
