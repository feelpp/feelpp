/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/** @file test_reference_tabulation.cpp
 * @brief Requested modal/nodal tables, retained storage, components and FE precompute compatibility.
 */
#define BOOST_TEST_MODULE reference_tabulation
#include <array>
#include <boost/mp11/list.hpp>
#include <boost/test/data/test_case.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>
#include <tuple>
#include <type_traits>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{

/** Cells share identical tabulation and coefficient-transform contracts. */
using Cells = boost::mp11::mp_list<Simplex<2>, Hypercube<2>>;

/** Fixed reference dimensions use the same modal test body. */
using Dimensions = boost::mp11::mp_list_c<int, 1, 2, 3>;

/** Whether a temporary owner would expose a dangling value reference. */
template <typename Tables>
concept TemporaryValues = requires( Tables&& tables ) { std::move( tables ).values(); };

/** @brief Column-oriented interior points, independent of basis degree and cardinality. */
template <int Dim>
auto referencePoints()
{
    // Keep the same interior samples for every order and cell-family comparison.
    Eigen::Matrix<double, Dim, 3> points;
    for ( int d = 0; d < Dim; ++d )
        for ( int q = 0; q < 3; ++q )
            points( d, q ) = -.9 + .03 * q + .01 * d;

    return points;
}

/** @brief Combined and selected requests agree with individual modal operations. */
template <int Dim, int Degree, template <int, int, typename, bool> class Family>
void checkModalTables()
{
    // Establish independent single-operation results before checking combined requests.
    Family<Dim, Degree, double, true> basis{ RuntimeOrder{ 2 } };
    auto points = referencePoints<Dim>();
    auto values = basis.values( points );
    auto gradients = basis.gradients( points );
    auto hessians = basis.hessians( points );

    // Exercise every request mask while retaining compile-time degree and sample extents.
    for ( unsigned mask = 1; mask <= 7; ++mask )
        withEvaluationRequest( static_cast<Evaluation>( mask ), [&]( auto selection )
                               {
            constexpr auto request = decltype( selection )::value;
            auto tables = tabulate<request>( basis, points );
            static_assert( decltype( tables )::order_spec == Degree );
            static_assert( decltype( tables )::rowsAtCompileTime == Family<Dim, Degree, double, true>::basisExtent );
            static_assert( decltype( tables )::samplesAtCompileTime == 3 );
            static_assert( !TemporaryValues<decltype( tables )> );

            // Compare requested groups and reject access to an absent value table.
            if constexpr ( includes( request, Evaluation::Values ) )
                BOOST_CHECK_SMALL( ( tables.values().matrix() - values.matrix() ).norm(), 1e-12 );
            else BOOST_CHECK_THROW( (void)tables.values(), std::logic_error );
            for ( int i = 0; i < Dim; ++i )
            {
                if constexpr ( includes( request, Evaluation::Gradients ) )
                    BOOST_CHECK_SMALL( ( tables.gradient( i ).matrix() - gradients[i].matrix() ).norm(), 1e-12 );
                if constexpr ( includes( request, Evaluation::Hessians ) )
                    for ( int j = 0; j < Dim; ++j )
                        BOOST_CHECK_SMALL( ( tables.hessian( i, j ).matrix() - hessians[i * Dim + j].matrix() ).norm(), 1e-12 );
            }

            // Reuse existing output storage with an Eigen coordinate expression.
            TabulationWorkspace<decltype( basis ), request, 3> workspace{ RuntimeOrder{ 2 } };
            auto const* retained = tables.table( 0 ).matrix().data();
            tabulateInto( basis, points + 0. * points, tables, workspace );
            BOOST_CHECK_EQUAL( retained, tables.table( 0 ).matrix().data() );

            // Empty batches retain their row shape and invalid directions are rejected.
            BOOST_CHECK_THROW( (void)tables.gradient( Dim ), std::out_of_range );
            BOOST_CHECK_THROW( (void)tables.hessian( -1, 0 ), std::out_of_range );
            Eigen::Matrix<double, Dim, Dynamic> empty( Dim, 0 );
            auto emptyTables = tabulate<request>( basis, empty );
            BOOST_CHECK_EQUAL( emptyTables.cols(), 0 );
            BOOST_CHECK_EQUAL( emptyTables.rows(), basis.basisCount() ); } );

    // Coordinates stored inside an output table must survive alias-safe evaluation.
    auto aliased = tabulate( basis, points );
    aliased.values().matrix().topRows( Dim ) = points;
    TabulationWorkspace<decltype( basis ), Evaluation::All, 3> workspace{ RuntimeOrder{ 2 } };
    tabulateInto( basis, aliased.values().matrix().topRows( Dim ), aliased, workspace );
    BOOST_CHECK_SMALL( ( aliased.values().matrix() - values.matrix() ).norm(), 1e-12 );
    for ( int i = 0; i < Dim; ++i )
        BOOST_CHECK_SMALL( ( aliased.gradient( i ).matrix() - gradients[i].matrix() ).norm(), 1e-12 );

    // Owning copies retain their data when the original tables are overwritten.
    auto saved = aliased;
    aliased.values().matrix().setZero();
    BOOST_CHECK_SMALL( ( saved.values().matrix() - values.matrix() ).norm(), 1e-12 );

    // Check empty owning results and reject invalid runtime request masks.
    Eigen::Matrix<double, Dim, Dynamic> empty( Dim, 0 );
    auto emptyTables = tabulate( basis, empty );
    BOOST_CHECK_EQUAL( emptyTables.cols(), 0 );
    BOOST_CHECK_EQUAL( emptyTables.rows(), basis.basisCount() );
    BOOST_CHECK_THROW( withTabulation( basis, points, static_cast<Evaluation>( 0 ), []( auto const& ) {} ), std::invalid_argument );
    BOOST_CHECK_THROW( withTabulation( basis, points, static_cast<Evaluation>( 8 ), []( auto const& ) {} ), std::invalid_argument );
}

/** @brief Preserve component packing, API return types, request selection and retained precompute storage. */
template <typename Cell, int Degree, template <uint16_type> class Shape, int EffectiveDegree = 2>
void checkNodalPrecompute()
{
    // Build equivalent nodal and legacy-coordinate inputs for every field shape.
    using Element = typename Lagrange<Degree, Shape, Continuous, PointSetEquiSpaced>::
        template Apply<2, 2, double, Cell>::type;
    auto element = makeOrderAware<Element>( RuntimeOrder{ EffectiveDegree } );
    using Source = typename Element::polynomialset_type;
    using Precompute = typename Element::precompute_type;
    auto source = std::make_shared<Source>( element );
    auto points = referencePoints<2>();
    typename Source::points_type legacy( 2, points.cols() );
    eigenMap( legacy ) = points;

    // Degree, sample extent and component layout remain independent public contracts.
    auto tables = element.tabulate( points );
    static_assert( decltype( tables )::order_spec == Degree );
    static_assert( decltype( tables )::rowsAtCompileTime == Dynamic );
    static_assert( decltype( tables )::samplesAtCompileTime == 3 );
    static_assert( decltype( tables )::components == Source::nComponents );
    BOOST_CHECK_EQUAL( tables.rows(), source->coefficients().rows() );
    BOOST_CHECK_EQUAL( tables.functionCount(), tables.rows() / Source::nComponents );

    // Compare the shared coefficient transform with the existing separate operations.
    auto individualValues = source->values( legacy );
    auto individualGradients = source->gradients( legacy );
    auto individualHessians = source->hessians( legacy );
    BOOST_CHECK_SMALL( ( tables.values().matrix() - individualValues ).norm(), 1e-11 );

    // Preserve exact native accessor return types and the legacy Hessian tensor type.
    Precompute precompute( source, legacy );
    using Matrix = eigen_matrix_xx_col_type<double>;
    static_assert( std::is_same_v<decltype( precompute.phiEigen() ), Matrix const&> );
    static_assert( std::is_same_v<decltype( precompute.gradientEigen( 0 ) ), Matrix const&> );
    static_assert( std::is_same_v<decltype( precompute.hessianEigen( 0, 0 ) ), Matrix const&> );
    using LegacyHessian = Eigen::TensorFixedSize<double, Eigen::Sizes<Precompute::h_comp1, Precompute::h_comp2, Precompute::h_comp3>>;
    static_assert( std::is_same_v<typename Precompute::h_type, LegacyHessian> );

    // Same-shape updates retain native and compatibility allocations.
    auto const* valueStorage = precompute.phiEigen().data();
    auto const* gradientStorage = precompute.gradientEigen( 0 ).data();
    auto const* legacyStorage = precompute.phi().data();
    precompute.update( legacy );
    BOOST_CHECK_EQUAL( valueStorage, precompute.phiEigen().data() );
    BOOST_CHECK_EQUAL( gradientStorage, precompute.gradientEigen( 0 ).data() );
    BOOST_CHECK_EQUAL( legacyStorage, precompute.phi().data() );

    // Native reference tables agree with the existing scalar/vector/tensor evaluators.
    BOOST_CHECK_SMALL( ( precompute.phiEigen() - individualValues ).norm(), 1e-11 );
    for ( int i = 0; i < 2; ++i )
    {
        BOOST_CHECK_SMALL( ( precompute.gradientEigen( i ) - individualGradients[i] ).norm(), 1e-11 );
        for ( int j = 0; j < 2; ++j )
            BOOST_CHECK_SMALL( ( precompute.hessianEigen( i, j ) - individualHessians[i * 2 + j] ).norm(), 1e-11 );
    }

    // Verify every interleaved component against its legacy tensor accessor.
    for ( Eigen::Index function = 0; function < tables.functionCount(); ++function )
        for ( int c1 = 0; c1 < Source::nComponents1; ++c1 )
            for ( int c2 = 0; c2 < Source::nComponents2; ++c2 )
                for ( int q = 0; q < points.cols(); ++q )
                {
                    auto const row = function * Source::nComponents + c1 * Source::nComponents2 + c2;
                    BOOST_CHECK_SMALL( precompute.phi( function, c1, c2, q ) - tables.values().matrix()( row, q ), 1e-12 );
                    for ( int i = 0; i < 2; ++i )
                    {
                        double gradient;
                        if constexpr ( Source::rank < 2 )
                            gradient = precompute.grad( function, c1, i, q );
                        else
                            gradient = precompute.grad( function, c1, c2, i, q );
                        BOOST_CHECK_SMALL( gradient - tables.gradient( i ).matrix()( row, q ), 1e-12 );
                        for ( int j = 0; j < 2; ++j )
                        {
                            double hessian;
                            if constexpr ( Source::rank == 0 ) hessian = precompute.hessian( function, i, j, q );
                            if constexpr ( Source::rank == 1 ) hessian = precompute.hessian( function, c1, i, j, q );
                            if constexpr ( Source::rank == 2 ) hessian = precompute.hessian( function, c1, c2, i, j, q );
                            BOOST_CHECK_SMALL( hessian - tables.hessian( i, j ).matrix()( row, q ), 1e-12 );
                        }
                    }
                }

    // Request changes expose selected groups and reject unrequested native access.
    for ( unsigned mask = 1; mask <= 7; ++mask )
    {
        auto const request = static_cast<Evaluation>( mask );
        precompute.update( legacy, request );
        BOOST_CHECK_EQUAL( static_cast<unsigned>( precompute.requestedTables() ), mask );
        if ( includes( request, Evaluation::Values ) )
            BOOST_CHECK_SMALL( ( precompute.phiEigen() - individualValues ).norm(), 1e-11 );
        else
            BOOST_CHECK_THROW( (void)precompute.phiEigen(), std::logic_error );
        if ( includes( request, Evaluation::Gradients ) )
            BOOST_CHECK_SMALL( ( precompute.gradientEigen( 0 ) - individualGradients[0] ).norm(), 1e-11 );
        else
            BOOST_CHECK_THROW( (void)precompute.gradientEigen( 0 ), std::logic_error );
        if ( includes( request, Evaluation::Hessians ) )
            BOOST_CHECK_SMALL( ( precompute.hessianEigen( 0, 1 ) - individualHessians[1] ).norm(), 1e-11 );
        else
        {
            BOOST_CHECK_THROW( (void)precompute.hessianEigen( 0, 1 ), std::logic_error );
            if constexpr ( Source::rank == 2 )
                BOOST_CHECK_THROW( (void)precompute.hessian( 0, 0, 0, 0, 1, 0 ), std::logic_error );
        }
    }

    // Copies own their data independently through empty updates and invalid requests.
    precompute.update( legacy, Evaluation::All );
    auto copy = precompute;
    BOOST_CHECK( copy.phiEigen().data() != precompute.phiEigen().data() );
    typename Source::points_type empty( 2, 0 );
    precompute.update( empty );
    BOOST_CHECK_EQUAL( precompute.phiEigen().cols(), 0 );
    BOOST_CHECK_EQUAL( precompute.gradientEigen( 0 ).cols(), 0 );
    BOOST_CHECK_EQUAL( precompute.hessianEigen( 0, 1 ).cols(), 0 );
    BOOST_CHECK_SMALL( ( copy.phiEigen() - individualValues ).norm(), 1e-11 );
    BOOST_CHECK_THROW( copy.update( legacy, static_cast<Evaluation>( 8 ) ), std::invalid_argument );
    BOOST_CHECK_THROW( (void)copy.gradientEigen( 2 ), std::out_of_range );
    BOOST_CHECK_THROW( (void)copy.hessianEigen( 0, 2 ), std::out_of_range );
}
} // namespace

/** @test Every requested group uses the same simplex/tensor numerical kernel in all dimensions. */
BOOST_AUTO_TEST_CASE_TEMPLATE( modalRequestedTables, Dimension, Dimensions )
{
    constexpr int dim = Dimension::value;
    checkModalTables<dim, 2, SimplexDubinerBasis>();
    checkModalTables<dim, Dynamic, SimplexDubinerBasis>();
    checkModalTables<dim, 2, TensorLegendreBasis>();
    checkModalTables<dim, Dynamic, TensorLegendreBasis>();
}

/** @test Static/runtime scalar, vector and tensor fields preserve their precompute API and components. */
BOOST_AUTO_TEST_CASE_TEMPLATE( nodalPrecomputeCompatibility, Cell, Cells )
{
    // P1 sampled Hessians and fused tables must agree across every field shape.
    checkNodalPrecompute<Cell, 1, Scalar, 1>();
    checkNodalPrecompute<Cell, Dynamic, Scalar, 1>();
    checkNodalPrecompute<Cell, 1, Vectorial, 1>();
    checkNodalPrecompute<Cell, Dynamic, Vectorial, 1>();
    checkNodalPrecompute<Cell, 1, Tensor2, 1>();
    checkNodalPrecompute<Cell, Dynamic, Tensor2, 1>();
    checkNodalPrecompute<Cell, 1, Tensor2Symm, 1>();
    checkNodalPrecompute<Cell, Dynamic, Tensor2Symm, 1>();

    // P2 exercises direct sampled Hessians through the same compatibility checks.
    checkNodalPrecompute<Cell, 2, Scalar>();
    checkNodalPrecompute<Cell, Dynamic, Scalar>();
    checkNodalPrecompute<Cell, 2, Vectorial>();
    checkNodalPrecompute<Cell, Dynamic, Vectorial>();
    checkNodalPrecompute<Cell, 2, Tensor2>();
    checkNodalPrecompute<Cell, Dynamic, Tensor2>();
    checkNodalPrecompute<Cell, 2, Tensor2Symm>();
    checkNodalPrecompute<Cell, Dynamic, Tensor2Symm>();
}

/** @test Every legacy accessor rejects a table group that was not requested. */
BOOST_DATA_TEST_CASE( unrequestedLegacyAccessors, boost::unit_test::data::xrange( 1u, 8u ), mask )
{
    // Use one real precompute object for each of the seven runtime request masks.
    using Element = typename Lagrange<2, Scalar, Continuous, PointSetEquiSpaced>::template Apply<2, 2, double, Simplex<2>>::type;
    using Source = typename Element::polynomialset_type;
    Element element;
    auto source = std::make_shared<Source>( element );
    typename Source::points_type points( 2, 1 );
    eigenMap( points ) = referencePoints<2>().template leftCols<1>();
    auto const request = static_cast<Evaluation>( mask );
    typename Element::precompute_type precompute( source, points, request );
    auto const& immutable = precompute;
    uint16_type const zero = 0;
    size_type const dof = 0;

    // Describe the accessors as data so every request uses the same exception check.
    auto accessors = std::tuple{
        std::tuple{ Evaluation::Values, "phi collection", [&]
                    { (void)immutable.phi(); } },
        std::tuple{ Evaluation::Values, "phi pointer", [&]
                    { (void)immutable.phiPtr(); } },
        std::tuple{ Evaluation::Values, "phi indexed", [&]
                    { (void)immutable.phi( zero, zero, zero, zero ); } },
        std::tuple{ Evaluation::Values, "phi scalar", [&]
                    { (void)immutable.phi( zero, zero, zero, zero, rank_t<0>{} ); } },
        std::tuple{ Evaluation::Values, "phi vector", [&]
                    { (void)immutable.phi( zero, zero, zero, zero, rank_t<1>{} ); } },
        std::tuple{ Evaluation::Values, "phi matrix", [&]
                    { (void)immutable.phi( zero, zero, zero, zero, rank_t<2>{} ); } },
        std::tuple{ Evaluation::Values, "phi const tensor", [&]
                    { (void)immutable.phi( zero, zero ); } },
        std::tuple{ Evaluation::Values, "phi mutable tensor", [&]
                    { (void)precompute.phi( zero, zero ); } },
        std::tuple{ Evaluation::Values, "phi logging", [&]
                    { precompute.print(); } },
        std::tuple{ Evaluation::Gradients, "gradient collection", [&]
                    { (void)immutable.grad(); } },
        std::tuple{ Evaluation::Gradients, "gradient pointer", [&]
                    { (void)immutable.gradPtr(); } },
        std::tuple{ Evaluation::Gradients, "gradient scalar/vector", [&]
                    { (void)immutable.grad( dof, zero, zero, zero ); } },
        std::tuple{ Evaluation::Gradients, "gradient matrix", [&]
                    { (void)immutable.grad( dof, zero, zero, zero, zero ); } },
        std::tuple{ Evaluation::Hessians, "Hessian collection", [&]
                    { (void)immutable.hessian(); } },
        std::tuple{ Evaluation::Hessians, "Hessian short index", [&]
                    { (void)immutable.hessian( zero, zero ); } },
        std::tuple{ Evaluation::Hessians, "Hessian dof index", [&]
                    { (void)immutable.hessian( dof, zero ); } },
        std::tuple{ Evaluation::Hessians, "Hessian scalar", [&]
                    { (void)immutable.hessian( dof, zero, zero, zero ); } },
        std::tuple{ Evaluation::Hessians, "Hessian vector", [&]
                    { (void)immutable.hessian( dof, zero, zero, zero, zero ); } } };

    // Missing groups must throw before any empty compatibility array is indexed.
    auto check = [&]( auto const& accessor )
    {
        if ( !includes( request, std::get<0>( accessor ) ) )
        {
            BOOST_TEST_CONTEXT( std::get<1>( accessor ) )
            {
                BOOST_CHECK_THROW( std::get<2>( accessor )(), std::logic_error );
            }
        }
    };
    std::apply( [&]( auto const&... accessor )
                { ( check( accessor ), ... ); }, accessors );
}

/** @test A constant modal space retains value normalization and exactly zero derivatives. */
BOOST_AUTO_TEST_CASE_TEMPLATE( constantReferenceTables, Dimension, Dimensions )
{
    constexpr int dim = Dimension::value;
    auto check = [&]( auto const& basis )
    {
        // Constant spaces preserve normalization and produce exact zero derivatives.
        auto points = referencePoints<dim>();
        auto tables = tabulate( basis, points );
        auto expected = basis.values( points );
        BOOST_CHECK_EQUAL( tables.rows(), 1 );
        BOOST_CHECK_SMALL( ( tables.values().matrix() - expected.matrix() ).norm(), 1e-15 );
        for ( int i = 0; i < dim; ++i )
        {
            BOOST_CHECK_EQUAL( tables.gradient( i ).matrix().norm(), 0. );
            for ( int j = 0; j < dim; ++j )
                BOOST_CHECK_EQUAL( tables.hessian( i, j ).matrix().norm(), 0. );
        }
    };

    // Apply the same invariant to both modal families and both order policies.
    check( SimplexDubinerBasis<dim, 0>{} );
    check( SimplexDubinerBasis<dim, Dynamic>{ RuntimeOrder{ 0 } } );
    check( TensorLegendreBasis<dim, 0>{} );
    check( TensorLegendreBasis<dim, Dynamic>{ RuntimeOrder{ 0 } } );
}

/** @test An independently specified quadratic and its derivatives survive the nodal transform. */
BOOST_AUTO_TEST_CASE_TEMPLATE( analyticNodalDerivatives, Cell, Cells )
{
    // Interpolate an explicitly specified quadratic without using derivative evaluators.
    using Element = typename Lagrange<Dynamic, Scalar, Continuous, PointSetEquiSpaced>::template Apply<2, 2, double, Cell>::type;
    Element element{ RuntimeOrder{ 2 } };
    Eigen::Matrix<double, 1, Dynamic> weights( 1, element.points().size2() );
    for ( Eigen::Index q = 0; q < weights.cols(); ++q )
    {
        double const x = element.points()( 0, q ), y = element.points()( 1, q );
        weights( q ) = 1 + x + x * y + y * y;
    }

    // Evaluate the polynomial and compare against analytic first and second partials.
    typename Element::polynomial_type polynomial( element, weights );
    auto points = referencePoints<2>();
    auto tables = polynomial.tabulate( points );
    for ( int q = 0; q < points.cols(); ++q )
    {
        double const x = points( 0, q ), y = points( 1, q );
        BOOST_CHECK_SMALL( tables.values().matrix()( 0, q ) - ( 1 + x + x * y + y * y ), 1e-11 );
        BOOST_CHECK_SMALL( tables.gradient( 0 ).matrix()( 0, q ) - ( 1 + y ), 1e-11 );
        BOOST_CHECK_SMALL( tables.gradient( 1 ).matrix()( 0, q ) - ( x + 2 * y ), 1e-11 );
        BOOST_CHECK_SMALL( tables.hessian( 0, 0 ).matrix()( 0, q ), 1e-11 );
        BOOST_CHECK_SMALL( tables.hessian( 0, 1 ).matrix()( 0, q ) - 1, 1e-11 );
        BOOST_CHECK_SMALL( tables.hessian( 1, 1 ).matrix()( 0, q ) - 2, 1e-11 );
    }
}

/** @test Existing face maps and additive oriented generators select the correct reference coordinates. */
BOOST_AUTO_TEST_CASE_TEMPLATE( orientedReferenceRestrictions, Cell, Cells )
{
    // Use the established quadrature orientation maps as reference coordinates.
    using Element = typename Lagrange<2, Scalar, Continuous, PointSetEquiSpaced>::template Apply<2, 2, double, Cell>::type;
    using Source = typename Element::polynomialset_type;
    using Permutation = typename Source::convex_type::permutation_type;
    using EdgePermutation = typename Source::convex_type::edge_permutation_type;
    using Quadrature = std::conditional_t<Cell::is_simplex, IMGeneral<2, double, Simplex>, IMGeneral<2, double, Hypercube>>;
    Element element;
    auto source = std::make_shared<Source>( element );
    Quadrature quadrature( 4 );
    auto const& points = quadrature.allfpoints();
    auto maps = source->preComputeOnFaces( source, points );

    // Explicit face generators must preserve each permutation and derivative table.
    auto generator = [&]( uint16_type face, Permutation permutation ) -> typename Source::points_type const&
    { return points[face].at( permutation.value() ); };
    auto generated = source->preComputeOnFaces( source, generator );
    for ( std::size_t face = 0; face < maps.size(); ++face )
        for ( auto const& [permutation, precompute] : maps[face] )
        {
            auto const& expected = points[face].at( permutation.value() );
            auto direct = source->tabulate( eigenMap( expected ) );
            BOOST_CHECK_SMALL( ( eigenMap( precompute->nodes() ) - eigenMap( expected ) ).norm(), 1e-15 );
            BOOST_CHECK_SMALL( ( precompute->phiEigen() - direct.values().matrix() ).norm(), 1e-12 );
            BOOST_CHECK_SMALL( ( generated[face].at( permutation )->gradientEigen( 0 ) - direct.gradient( 0 ).matrix() ).norm(), 1e-12 );
            BOOST_CHECK_SMALL( ( generated[face].at( permutation )->hessianEigen( 0, 1 ) - direct.hessian( 0, 1 ).matrix() ).norm(), 1e-12 );
        }

    // Compare explicit and legacy edge generators with the same reversal convention.
    auto edgeGenerator = [&]( uint16_type edge, EdgePermutation permutation ) -> typename Source::points_type const&
    { return points[edge].at( permutation.value() ); };
    auto edges = source->preComputeOnEdges( source, edgeGenerator );
    auto legacyEdges = source->preComputeOnEdges( source, [&]( uint16_type edge ) -> typename Source::points_type const&
                                                  { return points[edge].at( EdgePermutation::IDENTITY ); } );
    for ( std::size_t edge = 0; edge < edges.size(); ++edge )
        for ( unsigned orientation : { EdgePermutation::IDENTITY, EdgePermutation::REVERSE_PERMUTATION } )
        {
            auto const key = Permutation{ static_cast<uint16_type>( orientation ) };
            auto const& expected = points[edge].at( orientation );
            BOOST_CHECK_SMALL( ( eigenMap( edges[edge].at( key )->nodes() ) - eigenMap( expected ) ).norm(), 1e-14 );
            BOOST_CHECK_SMALL( ( edges[edge].at( key )->phiEigen() - source->values( expected ) ).norm(), 1e-12 );
            BOOST_CHECK_SMALL( ( legacyEdges[edge].at( key )->phiEigen() - source->values( expected ) ).norm(), 1e-12 );
        }
}

/** @test Existing physical evaluators accept value-only precomputation and reject missing gradients. */
BOOST_AUTO_TEST_CASE( requestedTablesInPhysicalContext )
{
    // Keep mesh creation serial while the remaining reference tests run under MPI.
    if ( Environment::worldComm().localSize() > 1 )
    {
        BOOST_TEST_MESSAGE( "physical value-only context fixture is serial-only" );
        return;
    }

    // Build a real geometry context backed by value-only reference precomputation.
    using MeshType = Mesh<Simplex<2, 1>>;
    using Element = typename Lagrange<2, Scalar, Continuous, PointSetEquiSpaced>::template Apply<2, 2, double, Simplex<2>>::type;
    auto mesh = unitSquare( .5 );
    auto geometry = mesh->gm();
    auto basis = std::make_shared<Element>();
    IMGeneral<2, double, Simplex> quadrature( 3 );
    auto geometryPrecompute = geometry->preCompute( geometry, quadrature.points() );
    auto basisPrecompute = basis->preCompute( quadrature.points(), Evaluation::Values );

    using ValueContext = typename Element::template Context<vm::POINT, Element, typename MeshType::gm_type, typename MeshType::element_type>;
    using GradientContext = typename Element::template Context<vm::GRAD, Element, typename MeshType::gm_type, typename MeshType::element_type>;

    // Physical value evaluation consumes the requested tables on every element.
    int checkedElements = 0;
    for ( auto const& wrapped : elements( mesh ) )
    {
        auto const& element = unwrap_ref( wrapped );
        auto mapping = geometry->template context<vm::POINT | vm::JACOBIAN | vm::KB>( element, geometryPrecompute );
        ValueContext context( basis, mapping, basisPrecompute );
        auto const reference = basis->values( quadrature.points() );

        // Point evaluation must agree with the canonical reference values.
        for ( int q = 0; q < quadrature.points().size2(); ++q )
        {
            for ( int function = 0; function < reference.rows(); ++function )
            {
                BOOST_CHECK_SMALL( context.id( function, 0, 0, q ) - reference( function, q ), 1e-12 );
            }
        }

        // A gradient evaluator must reject an insufficient precomputation request.
        BOOST_CHECK_THROW( GradientContext( basis, mapping, basisPrecompute ), std::invalid_argument );
        ++checkedElements;
    }

    // Ensure that the mesh fixture actually exercised physical evaluation.
    BOOST_CHECK_GT( checkedElements, 0 );
}

BOOST_AUTO_TEST_CASE( owned_tables_outlive_input_coordinates )
{
    std::array<Eigen::MatrixXd, 7> expected;
    auto tables = [&]
    {
        SimplexDubinerBasis<2, Dynamic> basis{ RuntimeOrder{ 2 } };
        Eigen::Matrix<double, 2, Dynamic> points( 2, 3 );
        points << -.8, -.7, -.6, -.9, -.8, -.7;
        auto result = tabulate( basis, points );
        for ( std::size_t slot = 0; slot < expected.size(); ++slot )
            expected[slot] = result.table( slot ).matrix();
        points.setConstant( 23.0 );
        for ( std::size_t slot = 0; slot < expected.size(); ++slot )
            BOOST_CHECK_SMALL( ( result.table( slot ).matrix() - expected[slot] ).norm(), 1e-14 );
        return result;
    }();
    BOOST_CHECK_EQUAL( tables.order(), 2 );
    BOOST_CHECK_EQUAL( tables.rows(), 6 );
    BOOST_CHECK_EQUAL( tables.cols(), 3 );
    for ( std::size_t slot = 0; slot < expected.size(); ++slot )
        BOOST_CHECK_SMALL( ( tables.table( slot ).matrix() - expected[slot] ).norm(), 1e-14 );
    auto saved = tables;
    tables.values().matrix().setZero();
    BOOST_CHECK_SMALL( ( saved.values().matrix() - expected[0] ).norm(), 1e-14 );
}
