/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file test_polynomial_coefficients.cpp
 * @brief Native coefficient ownership, legacy snapshots and analytic FE tabulation.
 */
#define BOOST_TEST_MODULE polynomial_coefficients
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/functionalset.hpp>
#include <boost/mp11/list.hpp>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Both reference-cell families exercise the same coefficient API and mathematics. */
using Cells = boost::mp11::mp_list<Simplex<2>, Hypercube<2>>;

/** @brief Independent quadratic fixture shared by both order modes and cells. */
double polynomialValue( double x, double y ) { return 1. + x + x * y + y * y; }

/** @brief Check coefficient composition and all sampled derivatives against the fixture. */
template<typename Cell, int OrderSpec>
void checkPolynomialCoefficients()
{
    using Element = typename Lagrange<OrderSpec, Scalar, Continuous, PointSetEquiSpaced>::
        template Apply<2, 2, double, Cell>::type;
    auto element = makeOrderAware<Element>( RuntimeOrder{ 2 } );
    auto source = element.coefficients();
    static_assert( decltype( source )::ColsAtCompileTime ==
        ( OrderSpec == Dynamic ? Dynamic : ( Cell::is_simplex ? 6 : 9 ) ) );
    Eigen::Matrix<double, 1, Eigen::Dynamic> weights( 1, element.points().size2() );
    for ( Eigen::Index q = 0; q < weights.cols(); ++q )
        weights( q ) = polynomialValue( element.points()( 0, q ), element.points()( 1, q ) );
    typename Element::polynomial_type polynomial( element, weights );
    auto set = polynomial.toSet();
    typename Element::matrix_type points( 2, 3 );
    eigenMap( points ) << -.7, -.5, -.2, -.8, -.7, -.9;
    auto modalCoefficients = polynomial.coefficients().eval();
    using Set = typename Element::polynomialset_type;
    for ( bool asIs : { false, true } )
    {
        BOOST_TEST_CONTEXT( "asIs=" << asIs )
        {
            Eigen::MatrixXd input = weights;
            if ( asIs ) input = modalCoefficients;
            typename Element::matrix_type legacyInput( input.rows(), input.cols() );
            eigenMap( legacyInput ) = input;
            typename Element::polynomial_type nativePolynomial( element, input, asIs );
            typename Element::polynomial_type legacyPolynomial( element, legacyInput, asIs );
            Set nativeSet( element, input, asIs );
            Set legacySet( element, legacyInput, asIs );
            input.setZero();
            eigenMap( legacyInput ).setZero();
            auto checkValues = [&]( auto const& candidate )
            {
                BOOST_CHECK_EQUAL( candidate.order(), 2 );
                auto evaluated = candidate.values( points );
                for ( Eigen::Index q = 0; q < evaluated.cols(); ++q )
                    BOOST_CHECK_SMALL( evaluated( 0, q ) - polynomialValue( points( 0, q ), points( 1, q ) ), 1e-11 );
            };
            checkValues( nativePolynomial );
            checkValues( legacyPolynomial );
            checkValues( nativeSet );
            checkValues( legacySet );
        }
    }
    auto values = set.values( points );
    auto gradients = set.gradients( points );
    auto hessians = set.hessians( points );
    if constexpr ( Cell::is_hypercube )
    {
        auto coordinates = eigenMap<2, 3>( points );
        auto modal = set.basis().template values<OrderSpec>( coordinates, RuntimeOrder{ 2 } );
        static_assert( decltype( modal )::rowsAtCompileTime == ( OrderSpec == Dynamic ? Dynamic : 9 ) );
        static_assert( decltype( modal )::colsAtCompileTime == 3 );
        BOOST_CHECK_SMALL( ( set.coefficients() * modal.matrix() - values ).norm(), 1e-11 );
        auto modalGradients = set.basis().template gradients<OrderSpec>( coordinates, RuntimeOrder{ 2 } );
        static_assert( std::tuple_size_v<decltype( modalGradients )> == 2 );
        auto runtime = [&]
        {
            if constexpr ( OrderSpec == Dynamic ) return set;
            else return set.toDynamic();
        }();
        auto converted = runtime.template toStatic<2>();
        BOOST_REQUIRE( converted );
        BOOST_CHECK_SMALL( ( converted->values( points ) - values ).norm(), 1e-11 );
        ::Feel::detail::PolynomialBasisState<typename Element::basis_type, 2> state;
        auto fixedModal = state.withValues( set.basis(), coordinates, []( auto const& table )
        {
            static_assert( std::remove_cvref_t<decltype( table )>::RowsAtCompileTime == 9 );
            return table.eval();
        } );
        BOOST_CHECK_SMALL( ( set.coefficients() * fixedModal - values ).norm(), 1e-11 );
    }
    static_assert( std::tuple_size_v<decltype( gradients )> == 2 );
    static_assert( std::tuple_size_v<decltype( hessians )> == 4 );
    for ( Eigen::Index q = 0; q < values.cols(); ++q )
    {
        double const x = points( 0, q ), y = points( 1, q );
        BOOST_CHECK_SMALL( values( 0, q ) - polynomialValue( x, y ), 1e-11 );
        BOOST_CHECK_SMALL( gradients[0]( 0, q ) - ( 1. + y ), 1e-11 );
        BOOST_CHECK_SMALL( gradients[1]( 0, q ) - ( x + 2. * y ), 1e-11 );
        BOOST_CHECK_SMALL( hessians[0]( 0, q ), 1e-11 );
        BOOST_CHECK_SMALL( hessians[1]( 0, q ) - 1., 1e-11 );
        BOOST_CHECK_SMALL( hessians[2]( 0, q ) - 1., 1e-11 );
        BOOST_CHECK_SMALL( hessians[3]( 0, q ) - 2., 1e-11 );
    }
    auto snapshot = set.coeff();
    auto expected = set.coefficients().eval();
    eigenMap( snapshot ).setZero();
    BOOST_CHECK_SMALL( ( set.coefficients() - expected ).norm(), 1e-30 );
    auto copy = set;
    copy.setCoefficient( copy.coefficients() * 2., true );
    BOOST_CHECK_SMALL( ( copy.coefficients() - 2. * expected ).norm(), 1e-30 );
    BOOST_CHECK_SMALL( ( set.coefficients() - expected ).norm(), 1e-30 );
    // Self-aliasing assignment completes before replacing the owner.
    set.setCoefficient( set.coefficients().topRows( 1 ), true );
    BOOST_CHECK_SMALL( ( set.coefficients() - expected ).norm(), 1e-30 );
    if constexpr ( OrderSpec != Dynamic )
    {
        auto runtime = set.toDynamic();
        auto restored = runtime.template toStatic<2>();
        BOOST_REQUIRE( restored );
        BOOST_CHECK_SMALL( ( restored->coefficients() - expected ).norm(), 1e-30 );
        BOOST_CHECK_THROW( set.setCoefficient( Eigen::Matrix<double, 1, 2>::Zero(), true ), std::invalid_argument );
    }
}

/** @brief Differentiate twice without rebuilding the effective space's operators.
 * @return A derivative owning its state after all source copies are destroyed.
 */
template<typename Set>
Set checkDerivativeCache( Set source )
{
    using Basis = typename Set::basis_type;
    constexpr bool cached = Set::order_spec == Dynamic || Set::order_spec != Basis::nOrder;
    auto const* dx = &source.d( 0 );
    auto const* dy = &source.d( 1 );
    auto expected = ( source.coefficients() * eigenMap( *dx ) * eigenMap( *dy ) ).eval();
    auto first = source.derivative( 0 );
    auto second = first.derivative( 1 );
    BOOST_CHECK_EQUAL( second.order(), source.order() );
    BOOST_CHECK_SMALL( ( second.coefficients() - expected ).norm(), 1e-11 );
    if constexpr ( cached )
    {
        BOOST_CHECK( &first.d( 0 ) == dx );
        BOOST_CHECK( &first.d( 1 ) == dy );
        BOOST_CHECK( &second.d( 0 ) == dx );
        BOOST_CHECK( &second.d( 1 ) == dy );
    }
    if constexpr ( Set::order_spec == Dynamic )
    {
        auto changed = source;
        changed.setOrder( source.order() );
        BOOST_CHECK( &changed.d( 0 ) == dx );
        changed.setOrder( source.order() + 1 );
        auto const& different = changed.d( 0 );
        BOOST_CHECK( &different != dx );
        BOOST_CHECK_GT( different.size1(), dx->size1() );
        BOOST_CHECK( &source.d( 0 ) == dx );
        BOOST_CHECK( &second.d( 1 ) == dy );
        changed.setOrder( source.order() );
        // Returning to a live effective order reuses its immutable reference matrices.
        BOOST_CHECK( &changed.d( 0 ) == dx );
        BOOST_CHECK_SMALL( ( eigenMap( changed.d( 0 ) ) - eigenMap( *dx ) ).norm(), 1e-11 );
    }
    return second;
}
} // namespace

/** @test Scalar FE coefficients preserve order metadata, values, gradients and Hessians. */
BOOST_AUTO_TEST_CASE_TEMPLATE( nativeCoefficientOwnership, Cell, Cells )
{
    checkPolynomialCoefficients<Cell, 2>();
    checkPolynomialCoefficients<Cell, Dynamic>();
}

/** @test Cached operators survive derivative chains, order detachment and source destruction. */
BOOST_AUTO_TEST_CASE_TEMPLATE( derivativeOperatorStateIsRetained, Cell, Cells )
{
    using Element = typename Lagrange<Dynamic, Scalar, Continuous, PointSetEquiSpaced>::
        template Apply<2, 2, double, Cell>::type;
    using Set = typename Element::polynomialset_type;
    for ( int degree : { 1, 2, 3 } )
    {
        BOOST_TEST_CONTEXT( "degree=" << degree )
        {
            auto makeSet = [&]
            {
                auto element = makeOrderAware<Element>( RuntimeOrder{ degree } );
                return Set( static_cast<Set const&>( element ) );
            };
            auto survivor = checkDerivativeCache( makeSet() );
            BOOST_CHECK_EQUAL( survivor.d( 0 ).size1(), survivor.coefficients().cols() );
            auto polynomial = checkDerivativeCache( makeSet().polynomial( 0 ) );
            BOOST_CHECK_EQUAL( polynomial.d( 1 ).size1(), polynomial.coefficients().cols() );
            if ( degree == 2 )
            {
                auto converted = makeSet().template toStatic<2>();
                BOOST_REQUIRE( converted );
                auto fixed = checkDerivativeCache( std::move( *converted ) );
                BOOST_CHECK_SMALL( ( fixed.coefficients() - survivor.coefficients() ).norm(), 1e-11 );
                BOOST_CHECK_EQUAL( fixed.d( 1 ).size1(), fixed.coefficients().cols() );
            }
        }
    }
}

/** Static and runtime column extents share operator-composition ownership checks. */
using CoefficientColumns = boost::mp11::mp_list_c<int, 2, Dynamic>;

/** @test C D owns borrowed operands and rejects invalid products without changing either source. */
BOOST_AUTO_TEST_CASE_TEMPLATE( coefficientRightProductPreservesOwnership, Columns, CoefficientColumns )
{
    Eigen::Matrix2d initial;
    initial << 1., 2., 3., 4.;
    ::Feel::detail::PolynomialCoefficientStorage<double, Columns::value> coefficients( initial );
    ublas::matrix<float> legacyWeights( 2, 2 );
    legacyWeights( 0, 0 ) = 1.f;
    legacyWeights( 0, 1 ) = 2.f;
    legacyWeights( 1, 0 ) = 0.f;
    legacyWeights( 1, 1 ) = 1.f;
    auto legacyComposition = coefficients;
    legacyComposition.template set<Scalar<2>>( legacyWeights, false );
    Eigen::Matrix2d expectedLeft;
    expectedLeft << 7., 10., 3., 4.;
    BOOST_CHECK_SMALL( ( legacyComposition.matrix() - expectedLeft ).norm(), 1e-30 );
    legacyComposition.assign( legacyWeights + legacyWeights );
    expectedLeft << 2., 4., 0., 2.;
    BOOST_CHECK_SMALL( ( legacyComposition.matrix() - expectedLeft ).norm(), 1e-30 );
    auto square = coefficients.rightProduct( coefficients.matrix() );
    Eigen::Matrix2d expected;
    expected << 7., 10., 15., 22.;
    BOOST_CHECK_SMALL( ( square.matrix() - expected ).norm(), 1e-30 );
    BOOST_CHECK_SMALL( ( coefficients.matrix() - initial ).norm(), 1e-30 );
    Eigen::Matrix2d right;
    right << 1., 2., 0., 1.;
    auto composed = square.rightProduct( right );
    expected << 7., 24., 15., 52.;
    BOOST_CHECK_SMALL( ( composed.matrix() - expected ).norm(), 1e-30 );
    Eigen::MatrixXd wrong = Eigen::Matrix3d::Identity();
    BOOST_CHECK_THROW( (void)composed.rightProduct( wrong ), std::invalid_argument );
    BOOST_CHECK_SMALL( ( composed.matrix() - expected ).norm(), 1e-30 );
    Eigen::Matrix<double, 2, 3> rectangular = Eigen::Matrix<double, 2, 3>::Ones();
    if constexpr ( Columns::value == Dynamic )
    {
        auto result = composed.rightProduct( rectangular );
        BOOST_CHECK_SMALL( ( result.matrix() - expected * rectangular ).norm(), 1e-30 );
    }
    else
    {
        BOOST_CHECK_THROW( (void)composed.rightProduct( rectangular ), std::invalid_argument );
        BOOST_CHECK_SMALL( ( composed.matrix() - expected ).norm(), 1e-30 );
    }
}

/** @test Single-sample products preserve declared axes and handle strided/expression inputs. */
BOOST_AUTO_TEST_CASE( coefficientProductsRespectSampleAxes )
{
    Eigen::Matrix2d initial;
    initial << 1., 2., 3., 4.;
    ::Feel::detail::PolynomialCoefficientStorage<double, Dynamic> coefficients( initial );
    Eigen::MatrixXd column( 2, 1 );
    column << 1., 2.;
    auto values = coefficients.multiply( column );
    static_assert( decltype( values )::ColsAtCompileTime == Dynamic );
    BOOST_CHECK_EQUAL( values( 0, 0 ), 5. );
    BOOST_CHECK_EQUAL( values( 1, 0 ), 11. );
    auto legacy = coefficients.template multiplyAs<ublas::matrix<double>>( column );
    BOOST_CHECK_SMALL( ( eigenMap( legacy ) - values ).norm(), 1e-30 );
    Eigen::Matrix<double, 2, 4, Eigen::RowMajor> storage = Eigen::Matrix<double, 2, 4, Eigen::RowMajor>::Zero();
    storage.col( 0 ) = column;
    auto strided = coefficients.multiply( storage.leftCols( 1 ) );
    BOOST_CHECK_SMALL( ( strided - values ).norm(), 1e-30 );
    auto expression = coefficients.multiply( storage.leftCols( 1 ) * 2. );
    BOOST_CHECK_SMALL( ( expression - 2. * values ).norm(), 1e-30 );
    auto empty = coefficients.multiply( storage.leftCols( 0 ) );
    BOOST_CHECK_EQUAL( empty.rows(), 2 );
    BOOST_CHECK_EQUAL( empty.cols(), 0 );
    Eigen::MatrixXd wrong = Eigen::MatrixXd::Zero( 3, 1 );
    BOOST_CHECK_THROW( (void)coefficients.multiply( wrong ), std::invalid_argument );
}

/** @test Ordered dual assembly owns its coefficients and rejects inconsistent input shapes. */
BOOST_AUTO_TEST_CASE( nativeFunctionalCoefficients )
{
    using Element = typename Lagrange<2, Scalar, Continuous, PointSetEquiSpaced>::
        Apply<2, 2, double, Simplex<2>>::type;
    using Space = Element::primal_space_type;
    Space primal;
    Functional<Space> functional( primal );
    Eigen::Matrix<double, 1, 6> coefficients = Eigen::Matrix<double, 1, 6>::Zero();
    coefficients( 0 ) = 1.;
    functional.setCoefficient( coefficients );
    static_assert( decltype( functional.rieszRepresentation() )::ColsAtCompileTime == 6 );
    auto snapshot = functional.coeff();
    functional += functional;
    BOOST_CHECK_EQUAL( functional.rieszRepresentation()( 0, 0 ), 2. );
    BOOST_CHECK_EQUAL( snapshot( 0, 0 ), 1. );
    FunctionalSet<Space> dual( primal, { functional } );
    static_assert( decltype( dual.eigenDualMatrix() )::ColsAtCompileTime == 6 );
    BOOST_CHECK_EQUAL( dual.eigenDualMatrix()( 0, 0 ), 2. );
    functional.setCoefficient( 3. * coefficients );
    BOOST_CHECK_EQUAL( dual.eigenDualMatrix()( 0, 0 ), 2. );
    Functional<Space> invalid( primal );
    invalid.setCoefficient( Eigen::Matrix<double, 2, 6>::Zero() );
    BOOST_CHECK_THROW( dual.setFunctionalSet( { invalid } ), std::invalid_argument );
    BOOST_CHECK_EQUAL( dual.eigenDualMatrix()( 0, 0 ), 2. );
    dual.setFunctionalSet( {} );
    BOOST_CHECK_EQUAL( dual.eigenDualMatrix().rows(), 0 );
}
