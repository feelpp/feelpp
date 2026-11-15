/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file test_eigen_polynomial_storage.cpp
 * @brief Order, logical extent, allocation and view contracts for polynomial storage.
 */
#define BOOST_TEST_MODULE eigen_polynomial_storage
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpoly/eigenextent.hpp>
#include <feel/feelpoly/policy.hpp>
#include <feel/feelalg/eigenmap.hpp>
#include <boost/mp11/list.hpp>
#include <array>
#include <limits>
#include <type_traits>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
using detail::NodalCellKind;

/** @brief A reference cell/dimension dataset sharing every storage check. */
template<int Dim, NodalCellKind Kind>
struct StorageCase
{
    /** Compile-time policy at a selected degree. */
    template<int P> using Static = StorageEigen<double>::ForBasis<Dim, P, Kind>;
    /** Runtime-degree policy for the same cell and dimension. */
    using Runtime = Static<Dynamic>;
    /** @brief Independent low-degree dimension formula used by these datasets. */
    static constexpr Eigen::Index expectedCount( int p )
    {
        if constexpr ( Dim == 0 ) return 1;
        if constexpr ( Kind == NodalCellKind::Hypercube )
        {
            if constexpr ( Dim == 1 ) return p + 1;
            if constexpr ( Dim == 2 ) return ( p + 1 ) * ( p + 1 );
            return ( p + 1 ) * ( p + 1 ) * ( p + 1 );
        }
        else
        {
            if constexpr ( Dim == 1 ) return p + 1;
            if constexpr ( Dim == 2 ) return ( p + 1 ) * ( p + 2 ) / 2;
            return ( p + 1 ) * ( p + 2 ) * ( p + 3 ) / 6;
        }
    }
};

/** Point entities and higher-dimensional cells share the same factory checks. */
using StorageCases = boost::mp11::mp_list<
    StorageCase<0, NodalCellKind::Simplex>, StorageCase<1, NodalCellKind::Simplex>,
    StorageCase<2, NodalCellKind::Simplex>, StorageCase<3, NodalCellKind::Simplex>,
    StorageCase<0, NodalCellKind::Hypercube>, StorageCase<1, NodalCellKind::Hypercube>,
    StorageCase<2, NodalCellKind::Hypercube>, StorageCase<3, NodalCellKind::Hypercube>>;

/** Inline, heap-backed and independently dynamic axes share one shape contract. */
using BufferShapes = boost::mp11::mp_list<
    detail::EigenMatrixBuffer<double, 2, 3>,
    detail::EigenMatrixBuffer<double, 2, 3, 1>,
    detail::EigenMatrixBuffer<double, Dynamic, 3>,
    detail::EigenMatrixBuffer<double, 2, Dynamic>,
    detail::EigenMatrixBuffer<double, Dynamic, Dynamic>>;

/** @brief Check all factories, independent axes and static/runtime conversion. */
template<typename Case, int P>
void checkFactories()
{
    using Static = typename Case::template Static<P>;
    using Runtime = typename Case::Runtime;
    constexpr auto count = Case::expectedCount( P );
    static_assert( Static::basisExtent == count );
    static_assert( Runtime::basisExtent == Eigen::Dynamic );
    static_assert( Static::Values::rowsAtCompileTime == count );
    static_assert( Runtime::Values::rowsAtCompileTime == Eigen::Dynamic );
    static_assert( Static::Points::rowsAtCompileTime == Static::dimension );
    static_assert( std::is_empty_v<Static> );
    static_assert( Static{}.basisCount() == count );
    Static fixed;
    Runtime runtime{ RuntimeOrder{ P } };
    BOOST_CHECK_EQUAL( fixed.basisCount(), count );
    BOOST_CHECK_EQUAL( runtime.basisCount(), count );
    for ( Eigen::Index samples : std::array<Eigen::Index, 3>{ 0, 1, 5 } )
    {
        auto values = fixed.makeValues( samples );
        auto dynamicValues = runtime.makeValues( samples );
        BOOST_CHECK_EQUAL( values.rows(), count );
        BOOST_CHECK_EQUAL( values.cols(), samples );
        BOOST_CHECK_SMALL( values.matrix().squaredNorm(), 1e-30 );
        BOOST_CHECK_SMALL( dynamicValues.matrix().squaredNorm(), 1e-30 );
        for ( Eigen::Index row = 0; row < count; ++row )
            for ( Eigen::Index col = 0; col < samples; ++col )
                values.matrix()( row, col ) = 100. * row + col;
        dynamicValues = runtime.convertValues( fixed, values );
        auto restored = fixed.convertValues( runtime, dynamicValues );
        BOOST_CHECK_SMALL( ( restored.matrix() - values.matrix() ).norm(), 1e-30 );
        BOOST_CHECK_THROW( dynamicValues.resize( count + 1, samples ), std::invalid_argument );
        auto points = fixed.makePoints( samples );
        BOOST_CHECK_EQUAL( points.rows(), Static::dimension );
        BOOST_CHECK_EQUAL( points.cols(), samples );
        auto gradient = runtime.makeGradient( samples );
        auto hessian = runtime.makeHessian( samples );
        static_assert( std::tuple_size_v<decltype( gradient )> == Static::dimension );
        static_assert( std::tuple_size_v<decltype( hessian )> == Static::dimension * Static::dimension );
        for ( auto const& direction : gradient )
        {
            BOOST_CHECK_EQUAL( direction.rows(), count );
            BOOST_CHECK_EQUAL( direction.cols(), samples );
            BOOST_CHECK_SMALL( direction.matrix().squaredNorm(), 1e-30 );
        }
        for ( auto const& direction : hessian )
            BOOST_CHECK_SMALL( direction.matrix().squaredNorm(), 1e-30 );
    }
    for ( auto polynomials : std::array<Eigen::Index, 4>{ 0, count - 1, count, count + 2 } )
    {
        auto coefficients = runtime.makeCoefficients( polynomials );
        BOOST_CHECK_EQUAL( coefficients.rows(), polynomials );
        BOOST_CHECK_EQUAL( coefficients.cols(), count );
        BOOST_CHECK_SMALL( coefficients.matrix().squaredNorm(), 1e-30 );
        BOOST_CHECK_THROW( coefficients.resize( polynomials, count + 1 ), std::invalid_argument );
    }
    auto nodes = runtime.makeNodes();
    auto vector = fixed.makeVector();
    auto op = runtime.makeOperator();
    BOOST_CHECK_EQUAL( nodes.rows(), Static::dimension );
    BOOST_CHECK_EQUAL( nodes.cols(), count );
    BOOST_CHECK_EQUAL( vector.rows(), count );
    BOOST_CHECK_EQUAL( vector.cols(), 1 );
    BOOST_CHECK_EQUAL( op.rows(), count );
    BOOST_CHECK_EQUAL( op.cols(), count );
    BOOST_CHECK_SMALL( nodes.matrix().squaredNorm() + vector.matrix().squaredNorm() + op.matrix().squaredNorm(), 1e-30 );
    Runtime wrongOrder{ RuntimeOrder{ P + 1 } };
    BOOST_CHECK_THROW( (void)wrongOrder.convertValues( runtime, runtime.makeValues( 1 ) ), std::invalid_argument );
    auto destination = runtime.makeValues( 1 );
    auto wrongValues = wrongOrder.makeValues( 1 );
    if constexpr ( Case::expectedCount( P + 1 ) != count )
        BOOST_CHECK_THROW( destination = wrongValues, std::invalid_argument );
    BOOST_CHECK_THROW( (void)runtime.makeValues( -1 ), std::invalid_argument );
    BOOST_CHECK_THROW( (void)runtime.makeCoefficients( -1 ), std::invalid_argument );
    BOOST_CHECK_THROW( (void)fixed.valuesFrom( Eigen::MatrixXd::Zero( count + 1, 1 ) ), std::invalid_argument );
}

/** @brief Detect unsafe borrowed views of temporary owners by substitution. */
template<typename Buffer>
concept HasTemporaryMatrixView = requires( Buffer&& buffer ) { std::move( buffer ).matrix(); };
/** @brief Detect lifetime-unsafe adapters without instantiating a deleted function. */
template<typename Owner>
concept HasTemporaryEigenMap = requires( Owner&& owner ) { eigenMap( std::move( owner ) ); };
} // namespace

/** @test Cardinalities, shapes and initialization agree for both order modes. */
BOOST_AUTO_TEST_CASE_TEMPLATE( polynomialStorageFactories, Case, StorageCases )
{
    checkFactories<Case, 0>();
    checkFactories<Case, 1>();
    checkFactories<Case, 2>();
    checkFactories<Case, 3>();
    checkFactories<Case, 4>();
    checkFactories<Case, 8>();
}

/** @test A retained shape preserves storage; construction and changes still validate bindings. */
BOOST_AUTO_TEST_CASE_TEMPLATE( retainedBufferShape, Buffer, BufferShapes )
{
    BOOST_CHECK_THROW( (void)Buffer( 0, 0, 2, 3 ), std::invalid_argument );
    BOOST_CHECK_THROW( (void)Buffer::forOverwrite( 0, 0, 2, 3 ), std::invalid_argument );
    BOOST_CHECK_THROW( (void)Buffer( 2, 2, 2, 3 ), std::invalid_argument );
    Buffer buffer( 2, 3, 2, 3 );
    buffer.matrix().setConstant( 7. );
    auto const* data = buffer.matrix().data();
    buffer.resizeForOverwrite( 2, 3 );
    BOOST_CHECK_EQUAL( buffer.matrix().data(), data );
    BOOST_CHECK_EQUAL( buffer.matrix().sum(), 42. );
    BOOST_CHECK_THROW( buffer.resizeForOverwrite( 0, 0 ), std::invalid_argument );
    BOOST_CHECK_THROW( buffer.resizeForOverwrite( 2, 4 ), std::invalid_argument );
    BOOST_CHECK_EQUAL( buffer.matrix().sum(), 42. );
    buffer.resize( 2, 3 );
    BOOST_CHECK_SMALL( buffer.matrix().squaredNorm(), 1e-30 );
}

/** @test Heap fallback retains logical extents and checks assignment under NDEBUG. */
BOOST_AUTO_TEST_CASE( boundedStaticStorage )
{
    using Storage = detail::PolynomialStorage<double, 3, 4, NodalCellKind::Simplex, 64>;
    static_assert( Storage::basisExtent == 35 );
    static_assert( !Storage::Nodes::usesInlineStorage );
    static_assert( !Storage::Vector::usesInlineStorage );
    static_assert( !Storage::Operator::usesInlineStorage );
    static_assert( Storage::Operator::rowsAtCompileTime == 35 );
    static_assert( Storage::Operator::matrix_type::RowsAtCompileTime == 35 );
    static_assert( Storage::Operator::matrix_type::ColsAtCompileTime == Dynamic );
    static_assert( !HasTemporaryMatrixView<Storage::Operator> );
    Storage storage;
    auto op = storage.makeOperator();
    op.matrix().setIdentity();
    auto copy = op;
    auto moved = std::move( copy );
    BOOST_CHECK_SMALL( ( moved.matrix() - Eigen::MatrixXd::Identity( 35, 35 ) ).norm(), 1e-30 );
    BOOST_CHECK_THROW( moved.resize( 34, 35 ), std::invalid_argument );
    BOOST_CHECK_THROW( (void)Storage::Operator::forOverwrite( 34, 35 ), std::invalid_argument );
    BOOST_CHECK_THROW( moved.assign( Eigen::MatrixXd::Zero( 34, 35 ) ), std::invalid_argument );
    auto batch = storage.makeValues<32>( 32 );
    static_assert( decltype( batch )::colsAtCompileTime == 32 );
    static_assert( decltype( batch )::matrix_type::ColsAtCompileTime == Dynamic );
    detail::PolynomialStorage<double, 3, 4, NodalCellKind::Simplex> largerBudget;
    auto converted = largerBudget.convertValues( storage, batch );
    static_assert( decltype( converted )::colsAtCompileTime == 32 );
    static_assert( decltype( converted )::usesInlineStorage );
    BOOST_CHECK_SMALL( converted.matrix().squaredNorm(), 1e-30 );
    using Large = detail::PolynomialStorage<double, 2, 400, NodalCellKind::Simplex>;
    static_assert( Large::basisExtent == 80601 );
    Large large;
    auto nodes = large.makeNodes();
    auto vector = large.makeVector();
    BOOST_CHECK_EQUAL( nodes.cols(), 80601 );
    BOOST_CHECK_EQUAL( vector.rows(), 80601 );
    BOOST_CHECK_THROW( (void)detail::validatedNodalOrder<detail::SimplexCardinality>( 2, RuntimeOrder{ 400 } ), std::length_error );
}

/** @test Aliasing, strided inputs, layout and const views preserve coefficients. */
BOOST_AUTO_TEST_CASE( matrixViewsAndAliasing )
{
    detail::EigenMatrixBuffer<double, Dynamic, Dynamic> owner( 2, 3 );
    owner.matrix() << 1., 2., 3., 4., 5., 6.;
    auto expected = owner.matrix().transpose().eval();
    owner.assign( owner.matrix().transpose() );
    BOOST_CHECK_SMALL( ( owner.matrix() - expected ).norm(), 1e-30 );
    Eigen::Matrix<double, 6, 6, Eigen::RowMajor> source;
    for ( Eigen::Index i = 0; i < 6; ++i )
        for ( Eigen::Index j = 0; j < 6; ++j ) source( i, j ) = 10. * i + j;
    owner.assign( source.block( 1, 1, 3, 2 ) );
    BOOST_CHECK_SMALL( ( owner.matrix() - source.block( 1, 1, 3, 2 ) ).norm(), 1e-30 );
    auto const& constOwner = owner;
    static_assert( std::is_const_v<std::remove_pointer_t<decltype( constOwner.matrix().data() )>> );
    static_assert( detail::EigenMatrixBuffer<double, 1, Dynamic>::storageOrder == Eigen::RowMajor );
    static_assert( detail::EigenMatrixBuffer<double, Dynamic, 1>::storageOrder == Eigen::ColMajor );
    boost::numeric::ublas::matrix<double> legacy( 3, 2 );
    eigenMap( legacy ) = owner.matrix();
    auto borrowed = eigenMap( legacy );
    borrowed( 1, 1 ) = 91.;
    BOOST_CHECK_EQUAL( legacy( 1, 1 ), 91. );
    auto fixedRows = eigenMap<3, Dynamic>( legacy );
    static_assert( decltype( fixedRows )::RowsAtCompileTime == 3 );
    BOOST_CHECK_EQUAL( fixedRows( 1, 1 ), 91. );
    BOOST_CHECK_THROW( (void)( eigenMap<2, Dynamic>( legacy ) ), std::invalid_argument );
    static_assert( !HasTemporaryEigenMap<decltype( legacy )> );
    static_assert( !HasTemporaryEigenMap<ublas::vector<double>> );
}

/** @test Arithmetic fails before overflow or allocation; division does not overflow. */
BOOST_AUTO_TEST_CASE( checkedPolynomialAndMatrixExtents )
{
    auto const maxIndex = std::numeric_limits<Eigen::Index>::max();
    BOOST_CHECK_THROW( (void)detail::matrixElementCount( maxIndex, 2, sizeof( double ) ), std::length_error );
    BOOST_CHECK_THROW( (void)detail::matrixElementCount( maxIndex, 1, sizeof( double ) ), std::length_error );
    BOOST_CHECK_THROW( (void)detail::matrixElementCount( -1, 0, sizeof( double ) ), std::invalid_argument );
    BOOST_CHECK_THROW( (void)detail::staticEigenExtent( static_cast<std::size_t>( std::numeric_limits<int>::max() ) + 1 ), std::length_error );
    BOOST_CHECK_THROW( (void)detail::SimplexCardinality::total( 3, std::numeric_limits<int>::max() ), std::length_error );
    BOOST_CHECK_THROW( (void)detail::HypercubeCardinality::total( 3, std::numeric_limits<int>::max() ), std::length_error );
    BOOST_CHECK_THROW( (void)detail::validatedNodalOrder<detail::SimplexCardinality>( 2, RuntimeOrder{ 1 }, std::numeric_limits<size_type>::max() ), std::length_error );
    BOOST_CHECK_EQUAL( detail::checkedLocalIndexValue( 65535 ), 65535 );
    BOOST_CHECK_THROW( (void)detail::checkedLocalIndexValue( 65536 ), std::length_error );
    if constexpr ( sizeof( size_type ) == 8 )
    {
        BOOST_CHECK_EQUAL( detail::binomial( 67, 33 ), 14226520737620288370ULL );
        BOOST_CHECK_THROW( (void)detail::binomial( 68, 34 ), std::length_error );
    }
}
