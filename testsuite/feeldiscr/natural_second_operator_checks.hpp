/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#pragma once

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelvf/vf.hpp>
#include <limits>

namespace Feel::natural_second_checks
{
template <int Dim>
auto quadraticField()
{
    if constexpr ( Dim == 2 )
        return vec( Px() * Px() + Px() * Py(), Py() * Py() + 2 * Px() * Py() );
    else
        return vec( Px() * Px() + Px() * Py(), Py() * Py() + 2 * Px() * Py(),
                    Px() * Pz() + Py() * Pz() + Pz() * Pz() );
}

template <int Dim>
auto exactDivergence()
{
    if constexpr ( Dim == 2 )
        return 4 * Px() + 3 * Py();
    else
        return 5 * Px() + 4 * Py() + 2 * Pz();
}

template <int Dim>
auto exactCurl()
{
    if constexpr ( Dim == 2 )
        return 2 * Py() - Px();
    else
        return vec( Pz(), -Pz(), 2 * Py() - Px() );
}

template <int Dim>
auto exactGradDiv()
{
    if constexpr ( Dim == 2 )
        return vec( cst( 4. ), cst( 3. ) );
    else
        return vec( cst( 5. ), cst( 4. ), cst( 2. ) );
}

template <int Dim>
auto exactCurlCurl()
{
    if constexpr ( Dim == 2 )
        return vec( cst( 2. ), cst( 1. ) );
    else
        return vec( cst( 3. ), cst( 2. ), cst( 0. ) );
}

template <typename MeshPtr>
auto localPoints( MeshPtr const& mesh )
{
    using mesh_type = typename MeshPtr::element_type;
    std::vector<typename mesh_type::node_type> centers;
    for ( auto const& wrapped : elements( mesh ) )
    {
        centers.push_back( unwrap_ref( wrapped ).barycenter() );
        if ( centers.size() == 3 ) break;
    }
    typename matrix_node<double>::type points( mesh_type::nDim, centers.size() );
    for ( size_type q = 0; q < centers.size(); ++q )
        ublas::column( points, q ) = centers[centers.size() - q - 1];
    return points;
}

template <int Dim>
auto box( uint16_type order )
{
    using mesh_type = Mesh<Simplex<Dim, Dynamic>>;
    return createGMSHMesh(
        _mesh = new mesh_type( RuntimeOrder{ order } ),
        _desc = domain( _name = "natural-box-" + std::to_string( Dim ) + "-g" + std::to_string( order ),
                        _shape = "hypercube", _dim = Dim, _xmin = 0., _xmax = 1.,
                        _ymin = 0., _ymax = 1., _zmin = 0., _zmax = 1., _h = .7 ) );
}

template <int Dim>
auto curved( uint16_type order )
{
    using mesh_type = Mesh<Simplex<Dim, Dynamic>>;
    return createGMSHMesh(
        _mesh = new mesh_type( RuntimeOrder{ order } ),
        _desc = domain( _name = "natural-curved-" + std::to_string( Dim ) + "-g" + std::to_string( order ),
                        _shape = "ellipsoid", _dim = Dim, _xmin = -1., _ymin = -1., _zmin = -1., _h = .9 ) );
}

template <int Dim>
auto staticBox()
{
    return createGMSHMesh(
        _mesh = new Mesh<Simplex<Dim, 1>>,
        _desc = domain( _name = "natural-static-box-" + std::to_string( Dim ),
                        _shape = "hypercube", _dim = Dim, _xmin = 0., _xmax = 1.,
                        _ymin = 0., _ymax = 1., _zmin = 0., _zmax = 1., _h = .7 ) );
}

template <bool Assemble = true, typename SpacePtr>
void checkPolynomial( SpacePtr const& space )
{
    using element_type = typename SpacePtr::element_type::element_type;
    constexpr int Dim = element_type::nDim;
    auto mesh = space->mesh();
    auto field = space->element();
    auto exact = quadraticField<Dim>();
    field.on( _range = elements( mesh ), _expr = exact, _geomap = GeomapStrategyType::GEOMAP_HO );
    auto points = localPoints( mesh );
    typename element_type::matrix_node_type empty;
    typename element_type::id_array_type values( boost::extents[points.size2()] );
    auto checkValues = [&]( bool gradDiv )
    {
        for ( auto& value : values )
            value.setConstant( std::numeric_limits<double>::quiet_NaN() );
        if ( gradDiv )
            field.graddivInterpolate( points, values, false, empty );
        else
            field.curlcurlInterpolate( points, values, false, empty );
        constexpr double gd2[2] = { 4., 3. }, gd3[3] = { 5., 4., 2. };
        constexpr double cc2[2] = { 2., 1. }, cc3[3] = { 3., 2., 0. };
        for ( size_type q = 0; q < points.size2(); ++q )
            for ( int c = 0; c < Dim; ++c )
            {
                double expected = gradDiv ? ( Dim == 2 ? gd2[c] : gd3[c] ) : ( Dim == 2 ? cc2[c] : cc3[c] );
                BOOST_CHECK_SMALL( values[q]( c, 0 ) - expected, 2e-8 );
            }
    };
    BOOST_TEST_CONTEXT( space->basisName() << " FE" << space->order() << " G" << mesh->order() )
    {
        auto divError = divv( field ) - exactDivergence<Dim>();
        auto curlError = curlv( field ) - exactCurl<Dim>();
        double firstDerivativeError = integrate(
                                          _range = elements( mesh ),
                                          _expr = divError * divError + inner( curlError, curlError ),
                                          _quad = 12, _geomap = GeomapStrategyType::GEOMAP_HO )
                                          .evaluate()( 0, 0 );
        BOOST_CHECK_SMALL( std::sqrt( std::abs( firstDerivativeError ) ), 2e-8 );
        {
            checkValues( true );
            BOOST_CHECK_SMALL( normL2( _range = elements( mesh ), _expr = graddivv( field ) - exactGradDiv<Dim>(),
                                       _quad = 12, _geomap = GeomapStrategyType::GEOMAP_HO ),
                               2e-8 );
            BOOST_CHECK_SMALL( normL2( _range = boundaryfaces( mesh ), _expr = graddivv( field ) - exactGradDiv<Dim>(),
                                       _quad = 12, _geomap = GeomapStrategyType::GEOMAP_HO ),
                               2e-8 );
        }
        {
            checkValues( false );
            BOOST_CHECK_SMALL( normL2( _range = elements( mesh ), _expr = curlcurlv( field ) - exactCurlCurl<Dim>(),
                                       _quad = 12, _geomap = GeomapStrategyType::GEOMAP_HO ),
                               2e-8 );
            BOOST_CHECK_SMALL( normL2( _range = boundaryfaces( mesh ), _expr = curlcurlv( field ) - exactCurlCurl<Dim>(),
                                       _quad = 12, _geomap = GeomapStrategyType::GEOMAP_HO ),
                               2e-8 );
        }
        if constexpr ( Assemble )
        {
            // These are ordinary space elements: the trial/value suffixes
            // explicitly select their role in the element API.
            auto trial = space->element(), test = space->element();
            auto a = form2( _trial = space, _test = space );
            a = integrate( _range = elements( mesh ), _expr = inner( idt( trial ), id( test ) ), _quad = 12,
                           _geomap = GeomapStrategyType::GEOMAP_HO );
            double coefficient = 0;
            {
                a += integrate( _range = elements( mesh ), _expr = inner( graddivt( trial ), graddiv( test ) ), _quad = 12,
                                _geomap = GeomapStrategyType::GEOMAP_HO );
                coefficient += Dim == 2 ? 25. : 45.;
            }
            {
                a += integrate( _range = elements( mesh ), _expr = inner( curlcurlt( trial ), curlcurl( test ) ), _quad = 12,
                                _geomap = GeomapStrategyType::GEOMAP_HO );
                coefficient += Dim == 2 ? 5. : 13.;
            }
            double expected = integrate( _range = elements( mesh ), _expr = inner( exact, exact ) + cst( coefficient ),
                                         _quad = 12, _geomap = GeomapStrategyType::GEOMAP_HO )
                                  .evaluate()( 0, 0 );
            BOOST_CHECK_SMALL( a( field, field ) - expected, 2e-8 * ( 1 + std::abs( expected ) ) );
        }
    }
}

/** Independent physical finite differences of the existing div/curl API. */
template <typename SpacePtr>
void checkCurved( SpacePtr const& space )
{
    using element_type = typename SpacePtr::element_type::element_type;
    constexpr int Dim = element_type::nDim;
    auto mesh = space->mesh();
    auto field = space->element();
    field.on( _range = elements( mesh ), _expr = quadraticField<Dim>(), _geomap = GeomapStrategyType::GEOMAP_HO );
    auto points = localPoints( mesh );
    typename element_type::matrix_node_type empty;
    typename element_type::id_array_type actual( boost::extents[points.size2()] );
    for ( auto& value : actual )
        value.setConstant( std::numeric_limits<double>::quiet_NaN() );
    auto check = [&]<bool GradDiv>()
    {
        if constexpr ( GradDiv )
            field.graddivInterpolate( points, actual, false, empty );
        else
            field.curlcurlInterpolate( points, actual, false, empty );
        using first_array = std::conditional_t<GradDiv, typename element_type::div_array_type, typename element_type::curl_array_type>;
        std::vector<Eigen::Matrix<double, Dim, Dim>> derivatives( points.size2(), Eigen::Matrix<double, Dim, Dim>::Zero() );
        constexpr double step = 1e-3;
        for ( int d = 0; d < Dim; ++d )
        {
            std::array<first_array, 4> sampled;
            constexpr double offsets[4] = { 1., -1., .5, -.5 };
            for ( int sample = 0; sample < 4; ++sample )
            {
                auto shifted = points;
                for ( size_type q = 0; q < points.size2(); ++q )
                    shifted( d, q ) += step * offsets[sample];
                sampled[sample].resize( boost::extents[points.size2()] );
                for ( auto& value : sampled[sample] )
                    value.setConstant( std::numeric_limits<double>::quiet_NaN() );
                if constexpr ( GradDiv )
                    field.divInterpolate( shifted, sampled[sample], false, empty );
                else
                    field.curlInterpolate( shifted, sampled[sample], false, empty );
            }
            constexpr int components = GradDiv || Dim == 2 ? 1 : Dim;
            for ( size_type q = 0; q < points.size2(); ++q )
                for ( int c = 0; c < components; ++c )
                {
                    double coarse = ( sampled[0][q]( c, 0 ) - sampled[1][q]( c, 0 ) ) / ( 2 * step );
                    double fine = ( sampled[2][q]( c, 0 ) - sampled[3][q]( c, 0 ) ) / step;
                    derivatives[q]( c, d ) = ( 4 * fine - coarse ) / 3;
                }
        }
        for ( size_type q = 0; q < points.size2(); ++q )
        {
            Eigen::Matrix<double, Dim, 1> expected;
            if constexpr ( GradDiv )
                expected = derivatives[q].row( 0 ).transpose();
            else if constexpr ( Dim == 2 )
                expected << derivatives[q]( 0, 1 ), -derivatives[q]( 0, 0 );
            else
                expected << derivatives[q]( 2, 1 ) - derivatives[q]( 1, 2 ),
                    derivatives[q]( 0, 2 ) - derivatives[q]( 2, 0 ),
                    derivatives[q]( 1, 0 ) - derivatives[q]( 0, 1 );
            for ( int c = 0; c < Dim; ++c )
                BOOST_CHECK_SMALL( actual[q]( c, 0 ) - expected( c ), 2e-5 * ( 1 + std::abs( expected( c ) ) ) );
        }
    };
    BOOST_TEST_CONTEXT( space->basisName() << " curved G" << mesh->order() )
    {
        check.template operator()<true>();
        check.template operator()<false>();
    }
}
} // namespace Feel::natural_second_checks
