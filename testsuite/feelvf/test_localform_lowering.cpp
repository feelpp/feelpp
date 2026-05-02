#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_localform_lowering
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feeldiscr/pch.hpp>

#include <feel/feelvf/vf_eval.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/detail/localform.hpp>

using namespace Feel;
using namespace Feel::vf;

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( localform_lowering_suite )

BOOST_AUTO_TEST_CASE( scalar_mass_and_source )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();
    auto one = Xh->element();
    one.on( _range=elements( mesh ), _expr=cst( 1. ) );

    auto expr_mass = cst( 2.5 ) * idt( u ) * id( v );
    auto expr_source = cst( 3.0 ) * id( v );
    auto expr_unsupported = id( v ) + cst( 1.0 );

    auto a_generic = form2( _trial=Xh, _test=Xh );
    auto a_lowered = form2( _trial=Xh, _test=Xh );
    auto l_generic = form1( _test=Xh );
    auto l_lowered = form1( _test=Xh );

    using bilinear_form_type = std::decay_t<decltype( a_generic )>;
    using linear_form_type = std::decay_t<decltype( l_generic )>;
    constexpr bool canLowerMass =
        vf::detail::can_lower_scalar_bilinear_localform_for_form_v<decltype( expr_mass ), bilinear_form_type>;
    constexpr bool canLowerSource =
        vf::detail::can_lower_scalar_linear_localform_for_form_v<decltype( expr_source ), linear_form_type>;
    constexpr bool canLowerUnsupported =
        vf::detail::can_lower_scalar_linear_localform_v<decltype( expr_unsupported )>;

    BOOST_TEST( canLowerMass );
    BOOST_TEST( canLowerSource );
    BOOST_TEST( !canLowerUnsupported );

    auto lowered_mass = vf::detail::lower_scalar_bilinear_localform( expr_mass );
    auto lowered_source = vf::detail::lower_scalar_linear_localform( expr_source );

    a_generic = integrate( _range=elements( mesh ), _expr=expr_mass );
    a_lowered = integrate( _range=elements( mesh ), _expr=lowered_mass );
    l_generic = integrate( _range=elements( mesh ), _expr=expr_source );
    l_lowered = integrate( _range=elements( mesh ), _expr=lowered_source );

    auto measure = integrate( _range=elements( mesh ), _expr=cst( 1.0 ) ).evaluate()( 0, 0 );
    auto expected_mass = 2.5 * measure;
    auto expected_source = 3.0 * measure;

    BOOST_CHECK_CLOSE( a_generic( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_lowered( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_generic( one, one ), a_lowered( one, one ), 1e-10 );

    BOOST_CHECK_CLOSE( l_generic( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_lowered( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_generic( one ), l_lowered( one ), 1e-10 );
}

BOOST_AUTO_TEST_CASE( scalar_chi_weighted_mass_and_source )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();
    auto one = Xh->element();
    one.on( _range=elements( mesh ), _expr=cst( 1. ) );

    auto indicator = chi( Px() > cst( 0.5 ) );
    auto expr_mass = indicator * cst( 4.0 ) * idt( u ) * id( v );
    auto expr_source = indicator * cst( 1.25 ) * id( v );
    auto expr_nonlinear = chi( id( v ) > cst( 0.5 ) ) * id( v );

    static_assert( vf::detail::scalar_localform_coefficient_expr<decltype( indicator )>::supported );
    static_assert( vf::detail::scalar_localform_factors<decltype( indicator * cst( 4.0 ) )>::supported );
    static_assert( vf::detail::scalar_localform_factors<decltype( indicator * cst( 4.0 ) * idt( u ) )>::supported );

    auto a_generic = form2( _trial=Xh, _test=Xh );
    auto a_lowered = form2( _trial=Xh, _test=Xh );
    auto l_generic = form1( _test=Xh );
    auto l_lowered = form1( _test=Xh );

    using bilinear_form_type = std::decay_t<decltype( a_generic )>;
    using linear_form_type = std::decay_t<decltype( l_generic )>;
    constexpr bool canLowerMass =
        vf::detail::can_lower_scalar_bilinear_localform_for_form_v<decltype( expr_mass ), bilinear_form_type>;
    constexpr bool canLowerSource =
        vf::detail::can_lower_scalar_linear_localform_for_form_v<decltype( expr_source ), linear_form_type>;
    constexpr bool canLowerNonlinear =
        vf::detail::can_lower_scalar_linear_localform_v<decltype( expr_nonlinear )>;

    BOOST_TEST( canLowerMass );
    BOOST_TEST( canLowerSource );
    BOOST_TEST( !canLowerNonlinear );

    auto lowered_mass = vf::detail::lower_scalar_bilinear_localform( expr_mass );
    auto lowered_source = vf::detail::lower_scalar_linear_localform( expr_source );

    a_generic = integrate( _range=elements( mesh ), _expr=expr_mass );
    a_lowered = integrate( _range=elements( mesh ), _expr=lowered_mass );
    l_generic = integrate( _range=elements( mesh ), _expr=expr_source );
    l_lowered = integrate( _range=elements( mesh ), _expr=lowered_source );

    auto indicator_measure = integrate( _range=elements( mesh ), _expr=indicator ).evaluate()( 0, 0 );
    auto expected_source = 1.25 * indicator_measure;

    BOOST_TEST( indicator_measure > 0.45 );
    BOOST_TEST( indicator_measure < 0.55 );
    BOOST_TEST( a_generic( one, one ) > 0.0 );
    BOOST_TEST( a_lowered( one, one ) > 0.0 );
    BOOST_CHECK_CLOSE( a_generic( one, one ), a_lowered( one, one ), 1e-10 );

    BOOST_CHECK_CLOSE( l_generic( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_lowered( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_generic( one ), l_lowered( one ), 1e-10 );
}

BOOST_AUTO_TEST_CASE( scalar_field_coefficient_mass_and_source )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto Mh = Pch<1>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();
    auto one = Xh->element();
    auto beta = Mh->element();
    one.on( _range=elements( mesh ), _expr=cst( 1. ) );
    beta.on( _range=elements( mesh ), _expr=Px() + cst( 2. ) );

    auto expr_mass = idv( beta ) * idt( u ) * id( v );
    auto expr_source = ( idv( beta ) + cst( 0.5 ) ) * id( v );
    auto expr_unsupported = idv( beta ) * id( v ) * id( v );

    static_assert( vf::detail::scalar_localform_coefficient_expr<decltype( idv( beta ) )>::supported );
    static_assert( vf::detail::scalar_localform_coefficient_expr<decltype( idv( beta ) )>::kind ==
                   vf::detail::localform_scalar_expr_kind::scalar_value_leaf );

    auto a_generic = form2( _trial=Xh, _test=Xh );
    auto a_lowered = form2( _trial=Xh, _test=Xh );
    auto l_generic = form1( _test=Xh );
    auto l_lowered = form1( _test=Xh );

    using bilinear_form_type = std::decay_t<decltype( a_generic )>;
    using linear_form_type = std::decay_t<decltype( l_generic )>;
    constexpr bool canLowerMass =
        vf::detail::can_lower_scalar_bilinear_localform_for_form_v<decltype( expr_mass ), bilinear_form_type>;
    constexpr bool canLowerSource =
        vf::detail::can_lower_scalar_linear_localform_for_form_v<decltype( expr_source ), linear_form_type>;
    constexpr bool canLowerUnsupported =
        vf::detail::can_lower_scalar_bilinear_localform_v<decltype( expr_unsupported )>;

    BOOST_TEST( canLowerMass );
    BOOST_TEST( canLowerSource );
    BOOST_TEST( !canLowerUnsupported );

    auto lowered_mass = vf::detail::lower_scalar_bilinear_localform( expr_mass );
    auto lowered_source = vf::detail::lower_scalar_linear_localform( expr_source );

    a_generic = integrate( _range=elements( mesh ), _expr=expr_mass );
    a_lowered = integrate( _range=elements( mesh ), _expr=lowered_mass );
    l_generic = integrate( _range=elements( mesh ), _expr=expr_source );
    l_lowered = integrate( _range=elements( mesh ), _expr=lowered_source );

    auto expected_mass = integrate( _range=elements( mesh ), _expr=idv( beta ) ).evaluate()( 0, 0 );
    auto expected_source = integrate( _range=elements( mesh ), _expr=idv( beta ) + cst( 0.5 ) ).evaluate()( 0, 0 );

    BOOST_CHECK_CLOSE( a_generic( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_lowered( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_generic( one, one ), a_lowered( one, one ), 1e-10 );

    BOOST_CHECK_CLOSE( l_generic( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_lowered( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_generic( one ), l_lowered( one ), 1e-10 );
}

BOOST_AUTO_TEST_CASE( scalar_geometry_coefficient_mass_and_source )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();
    auto one = Xh->element();
    one.on( _range=elements( mesh ), _expr=cst( 1. ) );

    auto expr_mass = Px() * idt( u ) * id( v );
    auto expr_source = ( Px() + cst( 0.5 ) ) * id( v );
    auto expr_unsupported = Px() * id( v ) * id( v );

    static_assert( vf::detail::scalar_localform_coefficient_expr<decltype( Px() )>::supported );
    static_assert( vf::detail::scalar_localform_coefficient_expr<decltype( Px() )>::kind ==
                   vf::detail::localform_scalar_expr_kind::geometry_scalar_leaf );

    auto a_generic = form2( _trial=Xh, _test=Xh );
    auto a_lowered = form2( _trial=Xh, _test=Xh );
    auto l_generic = form1( _test=Xh );
    auto l_lowered = form1( _test=Xh );

    using bilinear_form_type = std::decay_t<decltype( a_generic )>;
    using linear_form_type = std::decay_t<decltype( l_generic )>;
    constexpr bool canLowerMass =
        vf::detail::can_lower_scalar_bilinear_localform_for_form_v<decltype( expr_mass ), bilinear_form_type>;
    constexpr bool canLowerSource =
        vf::detail::can_lower_scalar_linear_localform_for_form_v<decltype( expr_source ), linear_form_type>;
    constexpr bool canLowerUnsupported =
        vf::detail::can_lower_scalar_bilinear_localform_v<decltype( expr_unsupported )>;

    BOOST_TEST( canLowerMass );
    BOOST_TEST( canLowerSource );
    BOOST_TEST( !canLowerUnsupported );

    auto lowered_mass = vf::detail::lower_scalar_bilinear_localform( expr_mass );
    auto lowered_source = vf::detail::lower_scalar_linear_localform( expr_source );

    a_generic = integrate( _range=elements( mesh ), _expr=expr_mass );
    a_lowered = integrate( _range=elements( mesh ), _expr=lowered_mass );
    l_generic = integrate( _range=elements( mesh ), _expr=expr_source );
    l_lowered = integrate( _range=elements( mesh ), _expr=lowered_source );

    auto expected_mass = integrate( _range=elements( mesh ), _expr=Px() ).evaluate()( 0, 0 );
    auto expected_source = integrate( _range=elements( mesh ), _expr=Px() + cst( 0.5 ) ).evaluate()( 0, 0 );

    BOOST_CHECK_CLOSE( a_generic( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_lowered( one, one ), expected_mass, 1e-10 );
    BOOST_CHECK_CLOSE( a_generic( one, one ), a_lowered( one, one ), 1e-10 );

    BOOST_CHECK_CLOSE( l_generic( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_lowered( one ), expected_source, 1e-10 );
    BOOST_CHECK_CLOSE( l_generic( one ), l_lowered( one ), 1e-10 );
}

BOOST_AUTO_TEST_CASE( scalar_mixed_coefficient_mass_and_source )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto Mh = Pch<1>( mesh );

    auto u = Xh->element();
    auto v = Xh->element();
    auto one = Xh->element();
    auto beta = Mh->element();
    one.on( _range=elements( mesh ), _expr=cst( 1. ) );
    beta.on( _range=elements( mesh ), _expr=Px() + cst( 2. ) );

    auto indicator = chi( Px() > cst( 0.5 ) );
    auto mixed_coefficient = indicator * ( idv( beta ) + cst( 0.5 ) );
    auto expr_mass = mixed_coefficient * idt( u ) * id( v );
    auto expr_source = mixed_coefficient * id( v );
    auto expr_unsupported = indicator * ( id( v ) + cst( 0.5 ) ) * id( v );

    using coefficient_support = vf::detail::scalar_localform_coefficient_expr<decltype( mixed_coefficient )>;
    static_assert( coefficient_support::supported );
    static_assert( vf::detail::is_localform_coeff_ir_v<typename coefficient_support::coefficient_expr_type> );

    auto a_generic = form2( _trial=Xh, _test=Xh );
    auto a_lowered = form2( _trial=Xh, _test=Xh );
    auto l_generic = form1( _test=Xh );
    auto l_lowered = form1( _test=Xh );

    using bilinear_form_type = std::decay_t<decltype( a_generic )>;
    using linear_form_type = std::decay_t<decltype( l_generic )>;
    constexpr bool canLowerMass =
        vf::detail::can_lower_scalar_bilinear_localform_for_form_v<decltype( expr_mass ), bilinear_form_type>;
    constexpr bool canLowerSource =
        vf::detail::can_lower_scalar_linear_localform_for_form_v<decltype( expr_source ), linear_form_type>;
    constexpr bool canLowerUnsupported =
        vf::detail::can_lower_scalar_linear_localform_v<decltype( expr_unsupported )>;

    BOOST_TEST( canLowerMass );
    BOOST_TEST( canLowerSource );
    BOOST_TEST( !canLowerUnsupported );

    auto lowered_mass = vf::detail::lower_scalar_bilinear_localform( expr_mass );
    auto lowered_source = vf::detail::lower_scalar_linear_localform( expr_source );

    a_generic = integrate( _range=elements( mesh ), _expr=expr_mass );
    a_lowered = integrate( _range=elements( mesh ), _expr=lowered_mass );
    l_generic = integrate( _range=elements( mesh ), _expr=expr_source );
    l_lowered = integrate( _range=elements( mesh ), _expr=lowered_source );

    BOOST_TEST( a_generic( one, one ) > 0.0 );
    BOOST_TEST( a_lowered( one, one ) > 0.0 );
    BOOST_CHECK_CLOSE( a_generic( one, one ), a_lowered( one, one ), 1e-10 );

    BOOST_TEST( l_generic( one ) > 0.0 );
    BOOST_TEST( l_lowered( one ) > 0.0 );
    BOOST_CHECK_CLOSE( l_generic( one ), l_lowered( one ), 1e-10 );
}

BOOST_AUTO_TEST_SUITE_END()
