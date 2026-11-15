/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file testlagrangetransform.cpp
 * @brief Nodal unisolvence, reproduction, mass, component layout and transform ownership.
 */
#define BOOST_TEST_MODULE lagrangeTransform
#include <boost/mp11/list.hpp>
#include <boost/test/data/test_case.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>
#include <future>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Reference cells share every scalar reproduction and mass check. */
using Cells = boost::mp11::mp_list<Simplex<1>, Simplex<2>, Simplex<3>, Hypercube<2>>;

/** Floating-point precisions exercise the same solver contract. */
using Scalars = boost::mp11::mp_list<float, double>;

/** @brief Expose coordinate replacement to exercise arbitrary user-supplied interpolation nodes. */
template <typename Cell, int Degree, typename T = double>
struct AdjustableNodes : PointSetEquiSpaced<Cell, Degree, T>
{
    /** Underlying family retains its static/runtime degree contract. */
    using Base = PointSetEquiSpaced<Cell, Degree, T>;
    using Base::Base;
    using Base::setPoints;
};

/** @brief Check monomials, their derivatives and independently integrated modal/nodal mass. */
template <typename Element>
void checkScalarElement( Element const& element )
{
    constexpr int dim = Element::nDim;
    auto const degree = element.order();
    auto const& nodes = element.points();
    auto const cardinal = element.values( nodes );
    BOOST_CHECK_SMALL( ( cardinal - Eigen::MatrixXd::Identity( cardinal.rows(), cardinal.cols() ) ).norm(), 2e-9 );

    // Interior samples deliberately differ from both interpolation and quadrature nodes.
    Eigen::Matrix<double, dim, 3> points;
    for ( int d = 0; d < dim; ++d )
        for ( int q = 0; q < 3; ++q )
            points( d, q ) = -.87 + .02 * d + .03 * q;
    auto const tables = element.tabulate( points );

    // Coordinate monomials through the effective degree supply analytic first/second derivatives.
    for ( int direction = 0; direction < dim; ++direction )
        for ( int power = 0; power <= degree; ++power )
        {
            Eigen::VectorXd samples( nodes.size2() );
            for ( size_type node = 0; node < nodes.size2(); ++node )
                samples[node] = std::pow( nodes( direction, node ), power );
            for ( int q = 0; q < points.cols(); ++q )
            {
                double const x = points( direction, q );
                BOOST_CHECK_SMALL( samples.dot( tables.values().matrix().col( q ) ) - std::pow( x, power ), 2e-9 );
                for ( int i = 0; i < dim; ++i )
                {
                    double const first = i == direction && power > 0 ? power * std::pow( x, power - 1 ) : 0.;
                    BOOST_CHECK_SMALL( samples.dot( tables.gradient( i ).matrix().col( q ) ) - first, 2e-8 );
                    for ( int j = 0; j < dim; ++j )
                    {
                        double const second = i == direction && j == direction && power > 1 ? power * ( power - 1 ) * std::pow( x, power - 2 ) : 0.;
                        BOOST_CHECK_SMALL( samples.dot( tables.hessian( i, j ).matrix().col( q ) ) - second, 2e-7 );
                    }
                }
            }
        }

    // Sufficient quadrature checks normalization independently of the nodal solve.
    using Cell = std::conditional_t<Element::convex_type::is_simplex, Simplex<dim>, Hypercube<dim>>;
    PointSetQuadrature<Cell, double> quadrature( 2 * degree + 2 );
    auto const weights = eigenMap( quadrature.weights() );
    auto const modal = element.basis().values( quadrature.points(), RuntimeOrder{ degree } );
    Eigen::MatrixXd const modalMass = modal.matrix() * weights.asDiagonal() * modal.matrix().transpose();
    BOOST_CHECK_SMALL( ( modalMass - Eigen::MatrixXd::Identity( modalMass.rows(), modalMass.cols() ) ).norm(), 2e-9 );
    auto const nodal = element.values( quadrature.points() );
    Eigen::MatrixXd const nodalMass = nodal * weights.asDiagonal() * nodal.transpose();
    auto const coefficients = element.coefficients();
    BOOST_CHECK_SMALL( ( nodalMass - coefficients * coefficients.transpose() ).norm(), 2e-9 );
}

/** @brief Exercise equivalent static/runtime and continuous/discontinuous scalar elements. */
template <typename Cell, template <class, int, class> class Nodes>
void checkDegree( uint16_type degree )
{
    constexpr int dim = Cell::nDim;
    using DynamicElement = typename Lagrange<Dynamic, Scalar, Continuous, Nodes>::template Apply<dim, dim, double, Cell>::type;
    using DiscontinuousElement = typename Lagrange<Dynamic, Scalar, Discontinuous, Nodes>::template Apply<dim, dim, double, Cell>::type;
    DynamicElement dynamic{ RuntimeOrder{ degree } };
    DiscontinuousElement discontinuous{ RuntimeOrder{ degree } };
    checkScalarElement( dynamic );
    BOOST_CHECK_SMALL( ( dynamic.coefficients() - discontinuous.coefficients() ).norm(), 1e-13 );
    BOOST_CHECK( dynamic.dual().nodalTransform() == discontinuous.dual().nodalTransform() );

    // A static facade retains its modal column extent while sharing the same immutable data.
    if ( degree == 3 )
    {
        using StaticElement = typename Lagrange<3, Scalar, Continuous, Nodes>::template Apply<dim, dim, double, Cell>::type;
        StaticElement fixed;
        static_assert( decltype( fixed.dual().nodalCoefficients() )::ColsAtCompileTime != Dynamic );
        BOOST_CHECK( fixed.dual().nodalTransform() == dynamic.dual().nodalTransform() );
        BOOST_CHECK_SMALL( ( fixed.coefficients() - dynamic.coefficients() ).norm(), 1e-12 );
        checkScalarElement( fixed );
    }
}

/** @brief Verify component-major functions and interleaved output components against a scalar oracle. */
template <template <uint16_type> class Field>
void checkComponents()
{
    using Element = typename Lagrange<Dynamic, Field, Continuous, PointSetRecursiveLGL>::template Apply<2, 2, double, Simplex<2>>::type;
    using ScalarElement = typename Lagrange<Dynamic, Scalar, Continuous, PointSetRecursiveLGL>::template Apply<2, 2, double, Simplex<2>>::type;
    Element element{ RuntimeOrder{ 3 } };
    ScalarElement scalar{ RuntimeOrder{ 3 } };
    BOOST_CHECK( element.dual().nodalTransform() == scalar.dual().nodalTransform() );
    auto const expected = scalar.tabulate( scalar.points() );
    auto const actual = element.tabulate( scalar.points() );
    auto const nodes = scalar.points().size2();
    constexpr int components = Element::nComponents;
    for ( int functionComponent = 0; functionComponent < components; ++functionComponent )
        for ( size_type node = 0; node < nodes; ++node )
            for ( int component = 0; component < components; ++component )
            {
                auto const row = ( functionComponent * nodes + node ) * components + component;
                double const active = component == functionComponent ? 1. : 0.;
                BOOST_CHECK_SMALL( ( actual.values().matrix().row( row ) - active * expected.values().matrix().row( node ) ).norm(), 1e-12 );
                for ( int i = 0; i < 2; ++i )
                {
                    BOOST_CHECK_SMALL( ( actual.gradient( i ).matrix().row( row ) - active * expected.gradient( i ).matrix().row( node ) ).norm(), 1e-11 );
                    for ( int j = 0; j < 2; ++j )
                        BOOST_CHECK_SMALL( ( actual.hessian( i, j ).matrix().row( row ) - active * expected.hessian( i, j ).matrix().row( node ) ).norm(), 1e-10 );
                }
            }
}

/** @brief Compare cell traces with independently constructed lower-dimensional nodal bases. */
template <int Dim>
void checkConformingTrace( uint16_type degree )
{
    using Element = typename Lagrange<Dynamic, Scalar, Continuous, PointSetRecursiveLGL>::template Apply<Dim, Dim, double, Simplex<Dim>>::type;
    using Facet = typename Lagrange<Dynamic, Scalar, Continuous, PointSetRecursiveLGL>::template Apply<Dim - 1, Dim - 1, double, Simplex<Dim - 1>>::type;
    Element element{ RuntimeOrder{ degree } };
    Facet facet{ RuntimeOrder{ degree } };
    auto const nodeData = element.dual().pointSet().referenceData();

    // Embed arbitrary facet samples in the face opposite the last reference vertex.
    Eigen::Matrix<double, Dim - 1, 3> local;
    for ( int d = 0; d < Dim - 1; ++d )
        for ( int q = 0; q < 3; ++q )
            local( d, q ) = -.82 + .02 * d + .03 * q;
    Eigen::Matrix<double, Dim, 3> points;
    points.template topRows<Dim - 1>() = local;
    points.row( Dim - 1 ).setConstant( -1 );
    auto const trace = element.tabulate( points );
    auto const expected = facet.tabulate( local );

    // Boundary functions restrict to the matching facet basis; other functions vanish.
    for ( size_type node = 0; node < element.points().size2(); ++node )
    {
        if ( std::abs( nodeData->barycentric( Dim, node ) ) > 1e-13 )
        {
            BOOST_CHECK_SMALL( trace.values().matrix().row( node ).norm(), 1e-10 );
            continue;
        }
        Eigen::Index match;
        auto const coordinate = eigenMap( element.points() ).col( node ).template head<Dim - 1>().eval();
        double const distance = ( eigenMap( facet.points() ).colwise() - coordinate ).colwise().norm().minCoeff( &match );
        BOOST_REQUIRE_SMALL( distance, 1e-12 );
        BOOST_CHECK_SMALL( ( trace.values().matrix().row( node ) - expected.values().matrix().row( match ) ).norm(), 1e-10 );
        for ( int i = 0; i < Dim - 1; ++i )
        {
            BOOST_CHECK_SMALL( ( trace.gradient( i ).matrix().row( node ) - expected.gradient( i ).matrix().row( match ) ).norm(), 1e-9 );
            for ( int j = 0; j < Dim - 1; ++j )
                BOOST_CHECK_SMALL( ( trace.hessian( i, j ).matrix().row( node ) - expected.hessian( i, j ).matrix().row( match ) ).norm(), 1e-8 );
        }
    }

    // Reverse the local vertex order on the same physical face and compare oriented traces.
    auto reversed = points;
    reversed.row( 0 ) = ( 1 - ( points.array() + 1 ).colwise().sum() ).matrix();
    auto const other = element.values( reversed );
    std::array<uint16_type, Dim + 1> vertices;
    std::iota( vertices.begin(), vertices.end(), 0 );
    std::swap( vertices[0], vertices[1] );
    auto const permutation = element.dual().pointSet().nodePermutation( vertices );
    for ( size_type node = 0; node < element.points().size2(); ++node )
        BOOST_CHECK_SMALL( ( trace.values().matrix().row( node ) - other.row( permutation[node] ) ).norm(), 1e-10 );
}
} // namespace

/** @test One dataset covers P0 through P8 on simplex and tensor-product cells. */
BOOST_AUTO_TEST_CASE_TEMPLATE( equispacedReproduction, Cell, Cells )
{
    for ( uint16_type degree : { 0, 1, 2, 3, 5, 8 } )
        BOOST_TEST_CONTEXT( "dimension=" << Cell::nDim << " degree=" << degree )
        {
            checkDegree<Cell, PointSetEquiSpaced>( degree );
        }
}

/** @test Recursive nodes preserve reproduction and mass in every supported reference dimension. */
BOOST_DATA_TEST_CASE( recursiveReproduction, boost::unit_test::data::make( { 0, 1, 2, 3, 5, 8 } ), degree )
{
    checkDegree<Simplex<1>, PointSetRecursiveLGL>( degree );
    checkDegree<Simplex<2>, PointSetRecursiveLGL>( degree );
    checkDegree<Simplex<3>, PointSetRecursiveLGL>( degree );
}

/** @test Vector, tensor and symmetric-tensor fields use the same scalar transform and packing. */
BOOST_AUTO_TEST_CASE( componentOrdering )
{
    checkComponents<Vectorial>();
    checkComponents<Tensor2>();
    checkComponents<Tensor2Symm>();
}

/** @test Cell traces, tangential derivatives and reversed orientations agree with facet elements. */
BOOST_DATA_TEST_CASE( conformingNodalTraces, boost::unit_test::data::make( { 1, 2, 3, 5, 8 } ), degree )
{
    checkConformingTrace<2>( degree );
    checkConformingTrace<3>( degree );
}

/** @test The existing default point policies retain their numerical nodal representation. */
BOOST_DATA_TEST_CASE( defaultNodePolicies, boost::unit_test::data::make( { 1, 2, 3, 5, 8 } ), degree )
{
    checkDegree<Simplex<2>, PointSetFekete>( degree );
    checkDegree<Simplex<3>, PointSetFekete>( degree );
    checkDegree<Hypercube<2>, PointSetFekete>( degree );
}

/** @test LU and QR agree on a well-conditioned problem; conditioning alone does not reject a solve. */
BOOST_AUTO_TEST_CASE_TEMPLATE( checkedSolverStrategies, T, Scalars )
{
    using Transform = ::Feel::detail::NodalPolynomialTransform<T>;
    using Matrix = typename Transform::Matrix;
    using Method = ::Feel::detail::NodalSolveMethod;
    Matrix matrix( 2, 2 );
    matrix << 2, 1, 1, 3;
    Transform lu( matrix, Method::PartialPivLU );
    Transform qr( matrix, Method::ColumnPivQR );
    T const tolerance = T( 100 ) * Eigen::NumTraits<T>::epsilon();
    BOOST_CHECK_SMALL( ( lu.coefficients() - qr.coefficients() ).norm(), tolerance );
    BOOST_CHECK_EQUAL( qr.diagnostics().rank, 2 );

    // A diagonal, poorly scaled but full-rank system triggers QR and retains its diagnostic.
    matrix.setIdentity();
    matrix( 1, 1 ) = std::pow( Eigen::NumTraits<T>::epsilon(), T( .65 ) );
    Transform difficult( matrix );
    BOOST_CHECK( difficult.diagnostics().method == Method::ColumnPivQR );
    BOOST_CHECK( difficult.diagnostics().illConditioned );
    BOOST_CHECK_SMALL( ( difficult.coefficients() * matrix.transpose() - Matrix::Identity( 2, 2 ) ).norm(), tolerance );

    // Global scaling does not change rank, conditioning or the cardinal representation.
    for ( T scale : { T( 1e-20 ), T( 1e20 ) } )
    {
        Transform scaled( scale * matrix );
        BOOST_CHECK_SMALL( ( scaled.coefficients() * ( scale * matrix ).transpose() - Matrix::Identity( 2, 2 ) ).norm(), tolerance );
    }
    matrix.row( 1 ) = matrix.row( 0 );
    BOOST_CHECK_THROW( Transform{ matrix }, std::runtime_error );
    BOOST_CHECK_THROW( Transform{ Matrix::Zero( 2, 2 ) }, std::runtime_error );
    BOOST_CHECK_THROW( Transform{ Matrix::Zero( 2, 3 ) }, std::invalid_argument );
    BOOST_CHECK_THROW( Transform{ Matrix( 0, 0 ) }, std::invalid_argument );
    matrix( 0, 0 ) = std::numeric_limits<T>::infinity();
    BOOST_CHECK_THROW( Transform{ matrix }, std::invalid_argument );
}

/** @test Exact custom coordinates distinguish configurations, and failed construction cannot poison reuse. */
BOOST_AUTO_TEST_CASE( configurationIdentityAndLifetime )
{
    using Element = typename Lagrange<Dynamic, Scalar, Continuous, AdjustableNodes>::template Apply<1, 1, double, Simplex<1>>::type;
    using Primal = typename Element::primal_space_type;
    using Dual = typename Element::dual_space_type;
    Primal primal{ RuntimeOrder{ 3 } };
    Dual original( primal );
    auto retained = original.nodalTransform();
    auto const expected = retained->transform.coefficients().eval();
    std::vector<decltype( retained )> owners;

    // More live configurations than weak index slots exercise eviction without invalidating owners.
    for ( int i = 0; i < 20; ++i )
    {
        typename Dual::pointset_type nodes{ RuntimeOrder{ 3 } };
        auto points = nodes.points();
        points( 0, 2 ) += .001 * ( i + 1 );
        nodes.setPoints( points );
        Dual custom( primal, nodes );
        BOOST_CHECK( custom.nodalTransform() != retained );
        auto const values = custom.nodalTransform()->transform.coefficients() * custom.nodalTransform()->transform.vandermonde().transpose();
        BOOST_CHECK_SMALL( ( values - Eigen::MatrixXd::Identity( 4, 4 ) ).norm(), 1e-12 );
        owners.push_back( custom.nodalTransform() );
    }
    BOOST_CHECK_SMALL( ( retained->transform.coefficients() - expected ).norm(), 1e-14 );
    Dual copy( original );
    Dual moved( std::move( copy ) );
    BOOST_CHECK( moved.nodalTransform() == retained );

    // Duplicate coordinates fail before publication; subsequent valid construction still succeeds.
    typename Dual::pointset_type duplicate{ RuntimeOrder{ 3 } };
    auto points = duplicate.points();
    ublas::column( points, 2 ) = ublas::column( points, 3 );
    duplicate.setPoints( points );
    BOOST_CHECK_THROW( ( Dual{ primal, duplicate } ), std::runtime_error );
    BOOST_CHECK_NO_THROW( ( Dual{ primal } ) );

    // Altered full-rank primals retain the generic solve; singular primals still fail.
    /** @brief Concrete shell exposing the generic construction path for a modified primal. */
    struct TestElement : Feel::FiniteElement<Primal, fem::details::LagrangeDual, AdjustableNodes>
    {
        /** Common Ciarlet constructor supplied by the production base. */
        using Base = Feel::FiniteElement<Primal, fem::details::LagrangeDual, AdjustableNodes>;
        using Base::Base;
        /** @brief Name the test shell without altering its numerical representation. */
        std::string familyName() const override { return "Lagrange"; }
    };
    primal.setCoefficient( 2. * primal.coefficients(), true );
    Dual modified( primal );
    BOOST_CHECK( !modified.hasCanonicalPrimal() );
    TestElement finiteElement( modified );
    BOOST_CHECK_SMALL( ( finiteElement.values( finiteElement.points() ) - Eigen::MatrixXd::Identity( 4, 4 ) ).norm(), 1e-12 );
    primal.setCoefficient( Eigen::MatrixXd::Zero( 4, 4 ), true );
    BOOST_CHECK_THROW( ( TestElement{ Dual{ primal } } ), std::runtime_error );
}

/** @test Concurrent first construction publishes one immutable transform for every caller. */
BOOST_AUTO_TEST_CASE( concurrentReferenceConstruction )
{
    using Element = typename Lagrange<Dynamic, Scalar, Continuous, PointSetRecursiveLGL>::template Apply<2, 2, double, Simplex<2>>::type;
    using Handle = decltype( std::declval<Element const&>().dual().nodalTransform() );
    std::array<std::future<Handle>, 8> tasks;
    for ( auto& task : tasks )
        task = std::async( std::launch::async, []
                           {
            Element element{ RuntimeOrder{ 7 } };
            return element.dual().nodalTransform(); } );
    auto first = tasks[0].get();
    for ( size_type i = 1; i < tasks.size(); ++i )
        BOOST_CHECK( tasks[i].get() == first );
}
