/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file testreferenceoperators.cpp
 * @brief Analytic differentiation, immutable ownership, legacy construction and failure recovery.
 */
#define BOOST_TEST_MODULE referenceOperators
#include <atomic>
#include <boost/mp11/list.hpp>
#include <boost/test/data/test_case.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <future>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Fixed dimensions share one tensor differentiation and ownership contract. */
using Dimensions = boost::mp11::mp_list_c<int, 1, 2, 3>;

/** @brief Compare all coefficient directions against independently sampled native derivatives. */
template <int Dim, bool Normalized, template <class> class Storage = StorageUBlas>
void checkLegendre( uint16_type degree )
{
    using Modal = Legendre<Dim, Dim, 1, Feel::Normalized<Normalized>, double, Storage>;
    using OtherSeed = Legendre<Dim, Dim, 5, Feel::Normalized<Normalized>, double, Storage>;
    auto data = Modal::referenceOperators( RuntimeOrder{ degree } );
    BOOST_CHECK( data == OtherSeed::referenceOperators( RuntimeOrder{ degree } ) );
    BOOST_CHECK_EQUAL( data->degree, degree );

    // Interior and boundary samples are independent of the coefficient construction.
    Eigen::Matrix<double, Dim, 4> points;
    for ( int d = 0; d < Dim; ++d )
        for ( int q = 0; q < 4; ++q )
            points( d, q ) = -.9 + .5 * q + .02 * d;
    points.col( 0 ).setConstant( -1 );
    points.col( 3 ).setConstant( 1 );
    auto values = Modal::values( points, RuntimeOrder{ degree } );
    auto gradients = Modal::gradients( points, RuntimeOrder{ degree } );
    auto hessians = Modal::hessians( points, RuntimeOrder{ degree } );
    for ( int i = 0; i < Dim; ++i )
    {
        decltype( auto ) matrix = ::Feel::detail::referenceOperatorView( data->matrix( i ) );
        BOOST_CHECK_SMALL( ( matrix * values.matrix() - gradients[i].matrix() ).cwiseAbs().maxCoeff(), 2e-10 );
        auto copy = Modal::d( i, degree );
        BOOST_CHECK_SMALL( ( ::Feel::detail::referenceOperatorView( copy ) - matrix ).norm(), 1e-14 );
        for ( int j = 0; j < Dim; ++j )
        {
            decltype( auto ) other = ::Feel::detail::referenceOperatorView( data->matrix( j ) );
            BOOST_CHECK_SMALL( ( matrix * other * values.matrix() - hessians[i * Dim + j].matrix() ).cwiseAbs().maxCoeff(), 2e-8 );
        }
    }
    BOOST_CHECK_THROW( (void)data->matrix( Dim ), std::out_of_range );
    BOOST_CHECK_THROW( (void)Modal::d( Dim, degree ), std::out_of_range );
}

/** @brief Probe the reachable legacy base without reading derived state during base construction. */
template <typename Tag>
class LegacyProbe : public Basis<Tag, double>
{
    /** Legacy provider registration and static direction access. */
    using Base = Basis<Tag, double>;

  public:
    /** Stable scalar modal primal supplying the historical interpolation contract. */
    using Modal = Legendre<2, 2, 3>;
    /** Traits include the static interior differentiation point policy. */
    using traits_type = typename Modal::traits_type;
    /** Count sampled evaluations to distinguish registration from initialization. */
    inline static std::atomic<int> evaluations = 0;
    /** Inject one deterministic initialization failure. */
    inline static std::atomic<bool> fail = false;
    /** @brief Initialize derived members only after registering the provider. */
    LegacyProbe()
        : Base( *this ), M_ready( true ) {}
    /** @brief Evaluate a fully initialized primal or fail before data publication. */
    auto evaluate( typename Modal::points_type const& points ) const
    {
        if ( !M_ready ) throw std::logic_error( "legacy provider evaluated before derived initialization" );
        ++evaluations;
        if ( fail.exchange( false ) ) throw std::runtime_error( "injected reference initialization failure" );
        return Modal::evaluate( points );
    }
    /** @brief Supply sampled derivatives from the same independent primal. */
    auto derivate( typename Modal::points_type const& points ) const { return Modal::derivate( points ); }

  private:
    /** Sentinel detects evaluation from an unfinished base/derived constructor chain. */
    bool M_ready;
};

/** Independent identities keep first-use and failure tests isolated. */
struct LegacyConcurrentTag
{
};
/** Failure configuration cannot contaminate the successful first-use test. */
struct LegacyFailureTag
{
};
/** Cache test configuration has no dependency on modal-family caches. */
struct SmallCacheFamily
{
};
/** Failure publication test has its own bounded index. */
struct FailedCacheFamily
{
};
} // namespace

/** @test Normalized and raw tensor coefficients use the same exact formula in every dimension. */
BOOST_AUTO_TEST_CASE_TEMPLATE( tensorDerivativeIdentity, Dimension, Dimensions )
{
    for ( uint16_type degree : { 0, 1, 2, 3, 4, 5, 6, 8 } )
        BOOST_TEST_CONTEXT( "dimension=" << Dimension::value << " degree=" << degree )
        {
            checkLegendre<Dimension::value, true>( degree );
            checkLegendre<Dimension::value, false>( degree );
        }
}

/** @test Eigen compatibility storage uses identical normalized and raw coefficients. */
BOOST_DATA_TEST_CASE( nativeStorageParity, boost::unit_test::data::make( { 0, 1, 2, 3, 5 } ), degree )
{
    checkLegendre<2, true, StorageEigen>( degree );
    checkLegendre<2, false, StorageEigen>( degree );
}

/** @test Independent polynomial instances share operators and copies survive order changes. */
BOOST_AUTO_TEST_CASE( polynomialOwnershipAndOrders )
{
    using Modal = Legendre<2, 2, 1>;
    using Polynomial = PolynomialSet<Modal, Scalar, Dynamic>;
    Polynomial first( 3 ), second( 3 );
    auto const* original = &first.d( 0 );
    BOOST_CHECK_EQUAL( original, &second.d( 0 ) );
    auto copy = first;
    first.setOrder( 4 );
    BOOST_CHECK_NE( original, &first.d( 0 ) );
    BOOST_CHECK_EQUAL( original, &copy.d( 0 ) );
    BOOST_CHECK_EQUAL( original->size1(), 16 );
    BOOST_CHECK_EQUAL( first.d( 1 ).size1(), 25 );

    // Normalization is a mathematical cache input, not just a degree or point-count key.
    auto normalized = Modal::referenceOperators( RuntimeOrder{ 3 } );
    auto raw = Legendre<2, 2, 1, Normalized<false>>::referenceOperators( RuntimeOrder{ 3 } );
    BOOST_CHECK( normalized != raw );
    auto independent = Modal::d( 0, 3 );
    eigenMap( independent ).setZero();
    BOOST_CHECK_GT( eigenMap( normalized->matrix( 0 ) ).norm(), 0. );
}

/** @test Every simplex direction shares one runtime interpolation result and owned lifetime. */
BOOST_DATA_TEST_CASE( simplexRuntimeOwnership, boost::unit_test::data::make( { 0, 2, 3, 5 } ), degree )
{
    using Modal = Dubiner<2, 2, 1>;
    using OtherSeed = Dubiner<2, 2, 6>;
    auto data = Modal::referenceOperators( RuntimeOrder{ degree } );
    BOOST_CHECK( data == OtherSeed::referenceOperators( RuntimeOrder{ degree } ) );
    Eigen::Matrix<double, 2, 3> points;
    points << -.8, -.7, -.6, -.9, -.8, -.7;
    auto values = Modal::values( points, RuntimeOrder{ degree } );
    auto gradients = Modal::gradients( points, RuntimeOrder{ degree } );
    for ( int direction = 0; direction < 2; ++direction )
        BOOST_CHECK_SMALL( ( eigenMap( data->matrix( direction ) ) * values.matrix() - gradients[direction].matrix() ).norm(), 1e-10 );
}

/** @test Concurrent legacy first access publishes complete data only, and failures retry safely. */
BOOST_AUTO_TEST_CASE( legacyInitializationAndRetry )
{
    using Probe = LegacyProbe<LegacyConcurrentTag>;
    BOOST_CHECK_THROW( (void)Probe::d( 0 ), std::logic_error );
    {
        Probe instance;
        BOOST_CHECK_EQUAL( Probe::evaluations.load(), 0 );
    }

    // Registration retains no pointer to the now-destroyed descriptor.
    std::array<std::future<typename Probe::matrix_type const*>, 8> tasks;
    for ( auto& task : tasks )
        task = std::async( std::launch::async, []
                           { return &Probe::d( 0 ); } );
    auto const* first = tasks[0].get();
    for ( std::size_t i = 1; i < tasks.size(); ++i )
        BOOST_CHECK_EQUAL( first, tasks[i].get() );
    BOOST_CHECK_EQUAL( Probe::evaluations.load(), 1 );
    BOOST_CHECK_THROW( (void)Probe::d( 2 ), std::out_of_range );

    // A throwing local-static initializer must leave no partial matrices or initialized flag.
    using Failing = LegacyProbe<LegacyFailureTag>;
    Failing registered;
    Failing::fail = true;
    BOOST_CHECK_THROW( (void)Failing::d( 0 ), std::runtime_error );
    BOOST_CHECK_EQUAL( Failing::evaluations.load(), 1 );
    BOOST_CHECK_EQUAL( Failing::d( 1 ).size1(), 16 );
    BOOST_CHECK_EQUAL( Failing::evaluations.load(), 2 );
}

/** @test Concurrent requested-order first access retains one complete immutable owner. */
BOOST_AUTO_TEST_CASE( concurrentRuntimeOwnership )
{
    using Modal = Legendre<1, 1, 1>;
    using Handle = decltype( Modal::referenceOperators( RuntimeOrder{ 11 } ) );
    std::array<std::future<Handle>, 8> tasks;
    for ( auto& task : tasks )
        task = std::async( std::launch::async, []
                           { return Modal::referenceOperators( RuntimeOrder{ 11 } ); } );
    auto owner = tasks[0].get();
    for ( std::size_t i = 1; i < tasks.size(); ++i )
        BOOST_CHECK( owner == tasks[i].get() );
    BOOST_CHECK_EQUAL( owner->matrix( 0 ).size1(), 12 );
    BOOST_CHECK( owner == Modal::referenceOperators( RuntimeOrder{ 11 } ) );
}

/** @test Reachable simplex descriptors initialize completely and support independent sample counts. */
BOOST_AUTO_TEST_CASE_TEMPLATE( reachableLegacyDescriptor, Dimension, Dimensions )
{
    BoundaryAdapted<Dimension::value, 3> basis;
    for ( int count : { 1, 3, 7 } )
        BOOST_TEST_CONTEXT( "dimension=" << Dimension::value << " points=" << count )
        {
            typename decltype( basis )::points_type points( Dimension::value, count );
            for ( int direction = 0; direction < Dimension::value; ++direction )
                for ( int q = 0; q < count; ++q )
                    points( direction, q ) = -.95 + .04 * ( q + 1 ) + .01 * direction;
            auto values = basis.evaluate( points );
            auto gradients = basis.derivate( points );
            for ( int direction = 0; direction < Dimension::value; ++direction )
            {
                BOOST_CHECK_SMALL( ( eigenMap( basis.d( direction ) ) * eigenMap( values ) - eigenMap( gradients[direction] ) ).norm(), 1e-9 );
                BOOST_CHECK_EQUAL( &basis.d( direction ), &basis.derivate( uint16_type( direction ) ) );

                // Independent differences check the primal formulas, not just the interpolated operators.
                auto plus = points, minus = points;
                for ( int q = 0; q < count; ++q )
                {
                    plus( direction, q ) += 2e-6;
                    minus( direction, q ) -= 2e-6;
                }
                auto upper = basis.evaluate( plus ), lower = basis.evaluate( minus );
                BOOST_CHECK_SMALL( ( ( eigenMap( upper ) - eigenMap( lower ) ) / 4e-6 - eigenMap( gradients[direction] ) ).norm(), 2e-8 );
            }
        }
}

/** @test The collapsed-coordinate product rule agrees with independent value differences. */
BOOST_DATA_TEST_CASE( principalProductRule, boost::unit_test::data::make( { 1, 2, 3, 5, 6 } ), degree )
{
    Principal<double> principal( degree );
    typename Principal<double>::vector_type points( 3 );
    points[0] = -.7;
    points[1] = -.1;
    points[2] = .4;
    auto plus = points, minus = points;
    for ( int q = 0; q < 3; ++q )
    {
        plus[q] += 2e-6;
        minus[q] -= 2e-6;
    }
    auto upper = principal.evaluate_2( plus ), lower = principal.evaluate_2( minus );
    auto gradients = principal.derivate_2( points );
    for ( std::size_t i = 0; i < gradients.size(); ++i )
        BOOST_CHECK_SMALL( ( ( eigenMap( upper[i] ) - eigenMap( lower[i] ) ) / 4e-6 - eigenMap( gradients[i] ) ).norm(), 5e-7 );
}

/** @test Byte/entry eviction preserves retained views and failed builders cannot poison a key. */
BOOST_AUTO_TEST_CASE( boundedRetentionAndFailure )
{
    using Matrix = Eigen::MatrixXd;
    using Data = ::Feel::detail::ReferenceDerivativeOperators<Matrix, 1>;
    using Cache = ::Feel::detail::ReferenceOperatorCache<SmallCacheFamily, Data, 2, 128>;
    int builds = 0;
    auto get = [&]( int degree )
    {
        return Cache::get( RuntimeOrder{ degree }, [&]
                           { ++builds; return Data{ static_cast<uint16_type>( degree ), { Matrix::Constant( 4, 4, degree ) } }; } );
    };
    auto first = get( 2 );
    BOOST_CHECK( first == get( 2 ) );
    BOOST_CHECK_EQUAL( builds, 1 );
    std::weak_ptr<Data const> weak = first;
    auto other = get( 3 );
    auto last = get( 4 );
    BOOST_CHECK_EQUAL( first->matrix( 0 )( 0, 0 ), 2 );
    first.reset();
    BOOST_CHECK( weak.expired() );
    BOOST_CHECK( last == get( 4 ) );

    // Construction failure does not publish a sentinel or evict the preceding configuration.
    using FailingCache = ::Feel::detail::ReferenceOperatorCache<FailedCacheFamily, Data>;
    auto preceding = FailingCache::get( RuntimeOrder{ 2 }, []
                                        { return Data{ 2, { Matrix::Identity( 2, 2 ) } }; } );
    BOOST_CHECK_THROW( (void)FailingCache::get( RuntimeOrder{ 1 }, []() -> Data
                                                { throw std::runtime_error( "injected failure" ); } ),
                       std::runtime_error );
    BOOST_CHECK( preceding == FailingCache::get( RuntimeOrder{ 2 }, []() -> Data
                                                 { throw std::logic_error( "successful entry was evicted" ); } ) );
    auto recovered = FailingCache::get( RuntimeOrder{ 1 }, []
                                        { return Data{ 1, { Matrix::Identity( 2, 2 ) } }; } );
    BOOST_CHECK_EQUAL( recovered->matrix( 0 ).trace(), 2 );
}

/** @test High-order value-only construction never requests an interpolation solve or dense D. */
BOOST_AUTO_TEST_CASE( valueOnlyConstruction )
{
    Legendre<2, 2, 32> basis;
    Eigen::Matrix<double, 2, 1> point;
    point << -.7, .4;
    auto values = basis.values( point );
    BOOST_CHECK_EQUAL( values.rows(), 1089 );
    BOOST_CHECK( values.matrix().allFinite() );
}
