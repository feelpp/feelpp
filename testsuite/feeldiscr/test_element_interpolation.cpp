/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define BOOST_TEST_MODULE element interpolation
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelvf/vf.hpp>

#include <limits>
#include <source_location>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template <typename Array, typename Points, typename Evaluate, typename Exact>
void checkPoints( Points const& points, Evaluate&& evaluate, Exact&& exact,
                  std::source_location const& caller = std::source_location::current() )
{
    BOOST_TEST_CONTEXT( "interpolation called at line " << caller.line() )
    {
        Array values( boost::extents[points.size2()] );
        for ( auto& value : values )
            value.setConstant( std::numeric_limits<double>::quiet_NaN() );
        evaluate( values );
        for ( size_type q = 0; q < points.size2(); ++q )
            for ( int i = 0; i < values[q].dimension( 0 ); ++i )
                for ( int j = 0; j < values[q].dimension( 1 ); ++j )
                    BOOST_CHECK_SMALL( values[q]( i, j ) - exact( i, j, q ), 2e-8 );
    }
}

template <typename Space, typename Points>
void checkScalar( Space const& space, Points const& points )
{
    using element_type = typename Space::element_type::element_type;
    constexpr int dim = element_type::nDim;
    auto field = space->element();
    if constexpr ( dim == 2 )
        field.on( _range = elements( space->mesh() ), _expr = cst( 1. ) + Px() * Px() + Px() * Py() + 3 * Py() * Py() );
    else
        field.on( _range = elements( space->mesh() ), _expr = cst( 1. ) + Px() * Px() + Px() * Py() + 3 * Py() * Py() + 2 * Pz() * Pz() + Px() * Pz() );

    auto gradient = [&]( int d, size_type q )
    {
        double x = points( 0, q ), y = points( 1, q ), z = dim == 3 ? points( 2, q ) : 0;
        return d == 0 ? 2 * x + y + z : ( d == 1 ? x + 6 * y : 4 * z + x );
    };
    typename element_type::matrix_node_type empty;
    checkPoints<typename element_type::id_array_type>( points, [&]( auto& out )
                                                       { field.idInterpolate( points, out, false, empty ); }, [&]( int, int, size_type q )
                                                       {
            double x=points(0,q), y=points(1,q), z=dim == 3 ? points(2,q) : 0;
            return 1+x*x+x*y+3*y*y+2*z*z+x*z; } );
    checkPoints<typename element_type::grad_array_type>( points, [&]( auto& out )
                                                         { field.gradInterpolate( points, out, false, empty ); }, [&]( int, int j, size_type q )
                                                         { return gradient( j, q ); } );
    checkPoints<typename element_type::id_array_type>( points, [&]( auto& out )
                                                       { field.dxInterpolate( points, out, false, empty ); }, [&]( int, int, size_type q )
                                                       { return gradient( 0, q ); } );
    checkPoints<typename element_type::id_array_type>( points, [&]( auto& out )
                                                       { field.dyInterpolate( points, out, false, empty ); }, [&]( int, int, size_type q )
                                                       { return gradient( 1, q ); } );
    if constexpr ( dim == 3 )
        checkPoints<typename element_type::id_array_type>( points, [&]( auto& out )
                                                           { field.dzInterpolate( points, out, false, empty ); }, [&]( int, int, size_type q )
                                                           { return gradient( 2, q ); } );
    checkPoints<typename element_type::hess_array_type>( points, [&]( auto& out )
                                                         { field.hessInterpolate( points, out, false, empty ); }, []( int i, int j, size_type )
                                                         {
            if ( i == j ) return i == 0 ? 2. : ( i == 1 ? 6. : 4. );
            return i == 0 || j == 0 ? 1. : 0.; } );
    checkPoints<typename element_type::id_array_type>( points, [&]( auto& out )
                                                       { field.laplacianInterpolate( points, out, false, empty ); }, []( int, int, size_type )
                                                       { return dim == 2 ? 8. : 12.; } );
    for ( auto const& wrapped : elements( space->mesh() ) )
    {
        auto const& element = unwrap_ref( wrapped );
        typename element_type::grad_array_type gradients;
        typename element_type::hess_array_type hessians;
        typename element_type::id_array_type laplacians;
        field.grad( element, gradients );
        field.hess( element, hessians );
        field.laplacian( element, laplacians );
        auto geometry = field.geomapPtr( element );
        for ( size_type q = 0; q < geometry->nPoints(); ++q )
        {
            double x = geometry->xReal( q )[0], y = geometry->xReal( q )[1];
            double z = dim == 3 ? geometry->xReal( q )[2] : 0;
            BOOST_CHECK_SMALL( gradients[q]( 0, 0 ) - ( 2 * x + y + z ), 2e-8 );
            BOOST_CHECK_SMALL( gradients[q]( 0, 1 ) - ( x + 6 * y ), 2e-8 );
            BOOST_CHECK_SMALL( hessians[q]( 0, 0 ) - 2., 2e-8 );
            BOOST_CHECK_SMALL( hessians[q]( 0, 1 ) - 1., 2e-8 );
            BOOST_CHECK_SMALL( hessians[q]( 1, 1 ) - 6., 2e-8 );
            BOOST_CHECK_SMALL( laplacians[q]( 0, 0 ) - ( dim == 2 ? 8. : 12. ), 2e-8 );
        }
        break;
    }
    BOOST_CHECK_SMALL( normL2( _range = boundaryfaces( space->mesh() ),
                               _expr = dnv( field ) - gradv( field ) * N() ),
                       2e-8 );
    BOOST_CHECK_SMALL( normL2( _range = elements( space->mesh() ),
                               _expr = laplacianv( field ) - trace( hessv( field ) ) ),
                       2e-8 );

    // All ranks evaluate the same point, including ranks not owning it.
    typename element_type::mesh_type::node_type point( dim );
    for ( int d = 0; d < dim; ++d )
        point[d] = 0.5;
    auto globalGradient = field.grad( point );
    BOOST_CHECK_SMALL( globalGradient( 0, 0, 0 ) - ( dim == 2 ? 1.5 : 2. ), 2e-8 );
    BOOST_CHECK_SMALL( globalGradient( 0, 1, 0 ) - 3.5, 2e-8 );

    // No local points must still allow participation in ghost exchange.
    auto localPoints = points;
    if ( Environment::worldComm().globalRank() != 0 ) localPoints.resize( dim, 0, false );
    checkPoints<typename element_type::id_array_type>( localPoints, [&]( auto& out )
                                                       { field.dxInterpolate( localPoints, out, false, empty ); }, [&]( int, int, size_type q )
                                                       { return gradient( 0, q ); } );
    typename element_type::dn_array_type normalValues( boost::extents[points.size2()] );
    BOOST_CHECK_THROW( field.dnInterpolate( points, normalValues, false, empty ), std::logic_error );
}

template <typename Space, typename Points>
void checkVector( Space const& space, Points const& points )
{
    using element_type = typename Space::element_type::element_type;
    constexpr int dim = element_type::nDim;
    auto field = space->element();
    if constexpr ( dim == 2 )
        field.on( _range = elements( space->mesh() ), _expr = vec( Px() * Px() + Px() * Py(), Py() * Py() + 2 * Px() * Py() ) );
    else
        field.on( _range = elements( space->mesh() ), _expr = vec( Px() * Px() + Px() * Py(), Py() * Py() + 2 * Px() * Py(), Px() * Pz() + Py() * Pz() + Pz() * Pz() ) );
    auto gradient = [&]( int c, int d, size_type q )
    {
        double x = points( 0, q ), y = points( 1, q ), z = dim == 3 ? points( 2, q ) : 0;
        if ( c == 0 ) return d == 0 ? 2 * x + y : ( d == 1 ? x : 0. );
        if ( c == 1 ) return d == 0 ? 2 * y : ( d == 1 ? 2 * x + 2 * y : 0. );
        return d == 2 ? x + y + 2 * z : z;
    };
    auto curl = [&]( int c, size_type q )
    {
        if constexpr ( dim == 3 )
            if ( c < 2 ) return ( c == 0 ? 1. : -1. ) * points( 2, q );
        return 2 * points( 1, q ) - points( 0, q );
    };
    typename element_type::matrix_node_type empty;
    checkPoints<typename element_type::grad_array_type>( points, [&]( auto& out )
                                                         { field.gradInterpolate( points, out, false, empty ); }, gradient );
    checkPoints<typename element_type::grad_array_type>( points, [&]( auto& out )
                                                         { field.symm_gradInterpolate( points, out, false, empty ); }, [&]( int i, int j, size_type q )
                                                         { return 0.5 * ( gradient( i, j, q ) + gradient( j, i, q ) ); } );
    checkPoints<typename element_type::id_array_type>( points, [&]( auto& out )
                                                       { field.dxInterpolate( points, out, false, empty ); }, [&]( int i, int, size_type q )
                                                       { return gradient( i, 0, q ); } );
    checkPoints<typename element_type::id_array_type>( points, [&]( auto& out )
                                                       { field.dyInterpolate( points, out, false, empty ); }, [&]( int i, int, size_type q )
                                                       { return gradient( i, 1, q ); } );
    if constexpr ( dim == 3 )
        checkPoints<typename element_type::id_array_type>( points, [&]( auto& out )
                                                           { field.dzInterpolate( points, out, false, empty ); }, [&]( int i, int, size_type q )
                                                           { return gradient( i, 2, q ); } );
    checkPoints<typename element_type::div_array_type>( points, [&]( auto& out )
                                                        { field.divInterpolate( points, out, false, empty ); }, [&]( int, int, size_type q )
                                                        { double sum=0; for(int d=0;d<dim;++d) sum+=gradient(d,d,q); return sum; } );
    checkPoints<typename element_type::curl_array_type>( points, [&]( auto& out )
                                                         { field.curlInterpolate( points, out, false, empty ); }, [&]( int i, int, size_type q )
                                                         { return curl( i, q ); } );
    checkPoints<typename element_type::comp_curl_array_type>( points, [&]( auto& out )
                                                              { field.curlxInterpolate( points, out, false, empty ); }, [&]( int, int, size_type q )
                                                              { return curl( 0, q ); } );
    checkPoints<typename element_type::comp_curl_array_type>( points, [&]( auto& out )
                                                              { field.curlyInterpolate( points, out, false, empty ); }, [&]( int, int, size_type q )
                                                              { return curl( dim == 3 ? 1 : 0, q ); } );
    checkPoints<typename element_type::comp_curl_array_type>( points, [&]( auto& out )
                                                              { field.curlzInterpolate( points, out, false, empty ); }, [&]( int, int, size_type q )
                                                              { return curl( dim == 3 ? 2 : 0, q ); } );
    checkPoints<typename element_type::id_array_type>( points, [&]( auto& out )
                                                       { field.laplacianInterpolate( points, out, false, empty ); }, []( int, int, size_type )
                                                       { return 2.; } );
    BOOST_CHECK_SMALL( normL2( _range = elements( space->mesh() ),
                               _expr = symm_gradv( field ) - 0.5 * ( gradv( field ) + trans( gradv( field ) ) ) ),
                       2e-8 );
}

template <typename Space, typename Points>
void checkTensor( Space const& space, Points const& points )
{
    using element_type = typename Space::element_type::element_type;
    constexpr bool symmetric = element_type::is_tensor2symm;
    auto field = space->element();
    if constexpr ( symmetric )
        field.on( _range = elements( space->mesh() ), _expr = mat<2, 2>( Px() * Px(), 2 * Px() * Py(), 2 * Px() * Py(), Py() * Py() ) );
    else
        field.on( _range = elements( space->mesh() ), _expr = mat<2, 2>( Px() * Px(), Px() * Py(), 2 * Px() * Py(), Py() * Py() ) );
    typename element_type::matrix_node_type empty;
    checkPoints<typename element_type::id_array_type>( points, [&]( auto& out )
                                                       { field.idInterpolate( points, out, false, empty ); }, [&]( int i, int j, size_type q )
                                                       {
            if ( i == j ) return points(i,q)*points(i,q);
            return (symmetric || i == 1 ? 2. : 1.)*points(0,q)*points(1,q); } );
    checkPoints<typename element_type::div_array_type>( points, [&]( auto& out )
                                                        { field.divInterpolate( points, out, false, empty ); }, [&]( int i, int, size_type q )
                                                        { return ( i == 0 && !symmetric ? 3. : 4. ) * points( i, q ); } );
    auto image = [&]()
    {
        if constexpr ( symmetric )
            return Pdhms<Dynamic>( space->mesh(), RuntimeOrder{ 3 } );
        else
            return Pdhm<Dynamic>( space->mesh(), RuntimeOrder{ 3 } );
    }();
    auto interpolation = opInterpolation( _domainSpace = space, _imageSpace = image );
    auto transferred = image->element();
    interpolation->apply( field, transferred );
    BOOST_CHECK_SMALL( normL2( _range = elements( space->mesh() ), _expr = idv( transferred ) - idv( field ) ), 2e-8 );
}

template <int Dim>
void checkGeometry( int geometryOrder )
{
    using mesh_type = Mesh<Simplex<Dim, Dynamic>>;
    auto mesh = createGMSHMesh(
        _mesh = new mesh_type( RuntimeOrder::checked( geometryOrder ) ),
        _desc = domain( _name = "interpolation-" + std::to_string( Dim ) + "-G" + std::to_string( geometryOrder ),
                        _shape = "hypercube", _dim = Dim, _xmin = 0., _xmax = 1., _ymin = 0., _ymax = 1.,
                        _zmin = 0., _zmax = 1., _h = 0.6 ) );
    std::vector<typename mesh_type::node_type> centers;
    for ( auto const& element : elements( mesh ) )
        centers.push_back( unwrap_ref( element ).barycenter() );
    typename matrix_node<double>::type points( Dim, centers.size() );
    // Reverse the order so scattering cannot rely on element traversal order.
    for ( size_type q = 0; q < centers.size(); ++q )
        ublas::column( points, q ) = centers[centers.size() - 1 - q];
    checkScalar( Pch<3>( mesh ), points );
    checkVector( Pchv<3>( mesh ), points );
    for ( uint16_type order : { 2, 3 } )
    {
        checkScalar( Pch<Dynamic>( mesh, RuntimeOrder{ order } ), points );
        checkVector( Pchv<Dynamic>( mesh, RuntimeOrder{ order } ), points );
    }
    if constexpr ( Dim == 2 )
    {
        checkTensor( Pdhm<2>( mesh ), points );
        checkTensor( Pdhm<Dynamic>( mesh, RuntimeOrder{ 2 } ), points );
        checkTensor( Pdhms<2>( mesh ), points );
        checkTensor( Pdhms<Dynamic>( mesh, RuntimeOrder{ 2 } ), points );
    }
}
} // namespace

BOOST_AUTO_TEST_CASE( polynomial_operators_2d )
{
    for ( int geometryOrder : { 1, 2, 3 } )
        checkGeometry<2>( geometryOrder );
}

BOOST_AUTO_TEST_CASE( polynomial_operators_3d )
{
    for ( int geometryOrder : { 1, 2, 3 } )
        checkGeometry<3>( geometryOrder );
}
