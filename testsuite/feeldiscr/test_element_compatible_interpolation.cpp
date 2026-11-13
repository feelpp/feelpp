/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define BOOST_TEST_MODULE compatible element interpolation
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/bdmh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/ned2h.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelvf/vf.hpp>
#include <limits>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template <typename Array, typename Points, typename Evaluate, typename Exact>
void checkValues( Points const& points, Evaluate&& evaluate, Exact&& exact )
{
    Array values( boost::extents[points.size2()] );
    for ( auto& value : values )
        value.setConstant( std::numeric_limits<double>::quiet_NaN() );
    evaluate( values );
    for ( size_type q = 0; q < points.size2(); ++q )
        for ( int i = 0; i < values[q].dimension( 0 ); ++i )
            BOOST_CHECK_SMALL( values[q]( i, 0 ) - exact( i, q ), 2e-8 );
}

template <typename Space, typename Points>
void checkSpace( Space const& space, Points const& points )
{
    using element_type = typename Space::element_type::element_type;
    constexpr int dim = element_type::nDim;
    constexpr double gradient[3][3] = { { 1, 2, 0 }, { 3, -1, 0 }, { 4, 1, 1 } };
    auto field = space->element();
    if constexpr ( dim == 2 )
        field.on( _range = elements( space->mesh() ), _expr = vec( cst( 1. ) + Px() + 2 * Py(), cst( 2. ) + 3 * Px() - Py() ) );
    else
        field.on( _range = elements( space->mesh() ), _expr = vec( cst( 1. ) + Px() + 2 * Py(), cst( 2. ) + 3 * Px() - Py(), cst( 3. ) + 4 * Px() + Py() + Pz() ) );
    typename element_type::matrix_node_type empty;
    BOOST_TEST_CONTEXT( space->basisName() << " order " << space->order() << " G" << space->mesh()->order() )
    {
        checkValues<typename element_type::id_array_type>( points, [&]( auto& out )
                                                           { field.idInterpolate( points, out, false, empty ); }, [&]( int i, size_type q )
                                                           {double value=i+1;for(int d=0;d<dim;++d)value+=gradient[i][d]*points(d,q);return value; } );
        typename element_type::grad_array_type gradients( boost::extents[points.size2()] );
        for ( auto& value : gradients )
            value.setConstant( std::numeric_limits<double>::quiet_NaN() );
        field.gradInterpolate( points, gradients, false, empty );
        for ( size_type q = 0; q < points.size2(); ++q )
            for ( int i = 0; i < dim; ++i )
                for ( int j = 0; j < dim; ++j )
                    BOOST_CHECK_SMALL( gradients[q]( i, j ) - gradient[i][j], 2e-8 );
        checkValues<typename element_type::id_array_type>( points, [&]( auto& out )
                                                           { field.dxInterpolate( points, out, false, empty ); }, [&]( int i, size_type )
                                                           { return gradient[i][0]; } );
        checkValues<typename element_type::id_array_type>( points, [&]( auto& out )
                                                           { field.dyInterpolate( points, out, false, empty ); }, [&]( int i, size_type )
                                                           { return gradient[i][1]; } );
        if constexpr ( dim == 3 )
            checkValues<typename element_type::id_array_type>( points, [&]( auto& out )
                                                               { field.dzInterpolate( points, out, false, empty ); }, [&]( int i, size_type )
                                                               { return gradient[i][2]; } );
        if constexpr ( element_type::is_hdiv_conforming )
            checkValues<typename element_type::div_array_type>( points, [&]( auto& out )
                                                                { field.divInterpolate( points, out, false, empty ); }, []( int, size_type )
                                                                { return dim == 2 ? 0. : 1.; } );
        if constexpr ( element_type::is_hcurl_conforming )
        {
            constexpr double curl[3] = { 1, -4, 1 };
            checkValues<typename element_type::curl_array_type>( points, [&]( auto& out )
                                                                 { field.curlInterpolate( points, out, false, empty ); }, [&]( int i, size_type )
                                                                 { return curl[i]; } );
            checkValues<typename element_type::comp_curl_array_type>( points, [&]( auto& out )
                                                                      { field.curlxInterpolate( points, out, false, empty ); }, []( int, size_type )
                                                                      { return 1.; } );
            checkValues<typename element_type::comp_curl_array_type>( points, [&]( auto& out )
                                                                      { field.curlyInterpolate( points, out, false, empty ); }, []( int, size_type )
                                                                      { return dim == 2 ? 1. : -4.; } );
            checkValues<typename element_type::comp_curl_array_type>( points, [&]( auto& out )
                                                                      { field.curlzInterpolate( points, out, false, empty ); }, []( int, size_type )
                                                                      { return 1.; } );
        }
    }
}

template <int Dim>
void checkDimension()
{
    using mesh_type = Mesh<Simplex<Dim, Dynamic>>;
    for ( uint16_type geometryOrder : { 1, 2, 3 } )
    {
        auto mesh = createGMSHMesh( _mesh = new mesh_type( RuntimeOrder{ geometryOrder } ),
                                    _desc = domain( _name = "compatible-points-" + std::to_string( Dim ) + "-G" + std::to_string( geometryOrder ),
                                                    _shape = "hypercube", _dim = Dim, _xmin = 0., _xmax = 1., _ymin = 0., _ymax = 1.,
                                                    _zmin = 0., _zmax = 1., _h = 0.7 ) );
        std::vector<typename mesh_type::node_type> centers;
        for ( auto const& element : elements( mesh ) )
            centers.push_back( unwrap_ref( element ).barycenter() );
        typename matrix_node<double>::type points( Dim, centers.size() );
        for ( size_type q = 0; q < centers.size(); ++q )
            ublas::column( points, q ) = centers[centers.size() - q - 1];
        for ( uint16_type order : { 1, 2 } )
        {
            checkSpace( RTh<Dynamic>( mesh, RuntimeOrder{ order } ), points );
            checkSpace( BDMh<Dynamic>( mesh, RuntimeOrder{ order } ), points );
            checkSpace( Ned1h<Dynamic>( mesh, RuntimeOrder{ order } ), points );
            checkSpace( Ned2h<Dynamic>( mesh, RuntimeOrder{ order } ), points );
        }
    }
}
template <typename Space>
void checkCurvedSpace( Space const& space )
{
    using element_type = typename Space::element_type::element_type;
    constexpr size_type geometryFlags = vm::POINT | vm::JACOBIAN | vm::KB | vm::HESSIAN;
    constexpr size_type naturalFlags = element_type::is_hdiv_conforming ? vm::DIV : vm::CURL;
    constexpr size_type basisFlags = vm::POINT | vm::JACOBIAN | vm::KB | vm::GRAD | vm::FIRST_DERIVATIVE | naturalFlags;
    auto field = space->element();
    field.on( _range = elements( space->mesh() ), _expr = vec( cst( 1. ) + Px() + 2 * Py(), cst( 2. ) + 3 * Px() - Py() ), _geomap = GeomapStrategyType::GEOMAP_HO );
    typename element_type::matrix_node_type reference( 2, 1 );
    reference( 0, 0 ) = -1. / 3.;
    reference( 1, 0 ) = -1. / 3.;
    auto gm = space->gm();
    auto geometricPrecompute = std::make_shared<typename element_type::geopc_type>( gm, reference );
    auto basisPrecompute = std::make_shared<typename element_type::pc_type>( space->fe(), reference );
    std::vector<Eigen::Vector2d> coordinates, expectedValues;
    std::vector<Eigen::Matrix2d> expectedGradients;
    std::vector<double> expectedNatural;
    for ( auto const& wrapped : elements( space->mesh() ) )
    {
        auto const& element = unwrap_ref( wrapped );
        auto geometry = gm->template context<geometryFlags>( element, geometricPrecompute );
        auto context = std::make_shared<typename element_type::template fec_t<basisFlags, geometryFlags>>(
            space->fe(), geometry, basisPrecompute );
        context->update( geometry, basisPrecompute );
        auto value = field.id( *context );
        auto gradient = field.grad( *context );
        coordinates.emplace_back( geometry->xReal( 0 )[0], geometry->xReal( 0 )[1] );
        expectedValues.emplace_back( value( 0, 0, 0 ), value( 1, 0, 0 ) );
        Eigen::Matrix2d expected;
        for ( int i = 0; i < 2; ++i )
            for ( int j = 0; j < 2; ++j )
                expected( i, j ) = gradient( i, j, 0 );
        expectedGradients.push_back( expected );
        if constexpr ( element_type::is_hdiv_conforming )
            expectedNatural.push_back( field.div( *context )( 0, 0, 0 ) );
        else
            expectedNatural.push_back( field.curl( *context )( 0, 0, 0 ) );
    }
    typename element_type::matrix_node_type points( 2, coordinates.size() ), empty;
    for ( size_type q = 0; q < coordinates.size(); ++q )
        for ( int d = 0; d < 2; ++d )
            points( d, q ) = coordinates[q][d];
    checkValues<typename element_type::id_array_type>( points, [&]( auto& out )
                                                       { field.idInterpolate( points, out, false, empty ); }, [&]( int i, size_type q )
                                                       { return expectedValues[q][i]; } );
    typename element_type::grad_array_type gradients( boost::extents[points.size2()] );
    for ( auto& value : gradients )
        value.setConstant( std::numeric_limits<double>::quiet_NaN() );
    field.gradInterpolate( points, gradients, false, empty );
    for ( size_type q = 0; q < coordinates.size(); ++q )
        for ( int i = 0; i < 2; ++i )
            for ( int j = 0; j < 2; ++j )
                BOOST_CHECK_SMALL( gradients[q]( i, j ) - expectedGradients[q]( i, j ), 2e-8 );
    if constexpr ( element_type::is_hdiv_conforming )
        checkValues<typename element_type::div_array_type>( points, [&]( auto& out )
                                                            { field.divInterpolate( points, out, false, empty ); }, [&]( int, size_type q )
                                                            { return expectedNatural[q]; } );
    else
        checkValues<typename element_type::curl_array_type>( points, [&]( auto& out )
                                                             { field.curlInterpolate( points, out, false, empty ); }, [&]( int, size_type q )
                                                             { return expectedNatural[q]; } );
}
} // namespace

BOOST_AUTO_TEST_CASE( compatible_curved_point_contexts )
{
    for ( uint16_type geometryOrder : { 2, 3 } )
    {
        auto mesh = unitCircle( RuntimeOrder{ geometryOrder }, 0.8 );
        checkCurvedSpace( RTh<Dynamic>( mesh, RuntimeOrder{ 1 } ) );
        checkCurvedSpace( BDMh<Dynamic>( mesh, RuntimeOrder{ 1 } ) );
        checkCurvedSpace( Ned1h<Dynamic>( mesh, RuntimeOrder{ 1 } ) );
        checkCurvedSpace( Ned2h<Dynamic>( mesh, RuntimeOrder{ 1 } ) );
    }
}

BOOST_AUTO_TEST_CASE( compatible_operators_2d )
{
    checkDimension<2>();
}
BOOST_AUTO_TEST_CASE( compatible_operators_3d )
{
    checkDimension<3>();
}
