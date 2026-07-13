/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
       Date: 2026-06-07

  Copyright (C) 2026 Université de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#define BOOST_TEST_MODULE vf symbolic ad expression testsuite

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( symbolic_ad_expr )

namespace
{
template <typename ExprT>
double
scalarValue( ExprT const& e )
{
    return e.evaluate()( 0, 0 );
}

bool
contains( std::string const& text, std::string const& needle )
{
    return text.find( needle ) != std::string::npos;
}

template <typename ExprT>
void
checkDerivativeValue( std::string const& name, ExprT const& expr, std::string const& symbol, double expected, double tolerance = 1e-12 )
{
    BOOST_TEST_CONTEXT( name )
    {
        auto derivative = expr.template diff<1>( symbol );
        BOOST_CHECK_SMALL( std::abs( scalarValue( derivative ) - expected ), tolerance );

        auto report = Feel::vf::describeDifferentiability( expr, symbol );
        BOOST_CHECK( report.ok() );
        BOOST_CHECK( report.supportsSymbolicDiff );
    }
}

template <typename ExprT>
double
finiteDifference( ExprT const& expr, std::string const& symbol, double value,
                  std::map<std::string,double> const& fixedParameters, double eps = 1e-6 )
{
    auto plus = expr;
    auto minus = expr;
    auto plusParameters = fixedParameters;
    auto minusParameters = fixedParameters;
    plusParameters[symbol] = value + eps;
    minusParameters[symbol] = value - eps;
    plus.setParameterValues( plusParameters );
    minus.setParameterValues( minusParameters );
    return ( scalarValue( plus ) - scalarValue( minus ) )/( 2.0*eps );
}
} // namespace

BOOST_AUTO_TEST_CASE( ginac_parameter_derivatives_match_analytic_and_finite_difference )
{
    using namespace Feel;
    using namespace Feel::vf;

    double const mu0 = 2.0;
    double const u0 = 0.25;
    double const eps = 1e-6;

    auto law = expr( "mu*(exp(u)-1):mu:u" );
    law.setParameterValues( { { "mu", mu0 }, { "u", u0 } } );

    auto dlawDu = law.diff<1>( "u" );
    auto dlawDmu = law.diff<1>( "mu" );

    BOOST_CHECK_SMALL( std::abs( scalarValue( dlawDu ) - mu0*std::exp( u0 ) ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( scalarValue( dlawDmu ) - ( std::exp( u0 ) - 1.0 ) ), 1e-12 );

    auto finiteDifference = [&law,eps]( std::string const& symbol, double value, std::string const& otherSymbol, double otherValue )
    {
        auto plus = law;
        auto minus = law;
        plus.setParameterValues( { { symbol, value + eps }, { otherSymbol, otherValue } } );
        minus.setParameterValues( { { symbol, value - eps }, { otherSymbol, otherValue } } );
        return ( scalarValue( plus ) - scalarValue( minus ) )/( 2.0*eps );
    };

    BOOST_CHECK_SMALL( std::abs( scalarValue( dlawDu ) - finiteDifference( "u", u0, "mu", mu0 ) ), 1e-9 );
    BOOST_CHECK_SMALL( std::abs( scalarValue( dlawDmu ) - finiteDifference( "mu", mu0, "u", u0 ) ), 1e-9 );
}

BOOST_AUTO_TEST_CASE( nested_symbols_and_apply_symbols_expr_keep_dependency_contract )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto e1 = expr( "u*u:u" );
    auto base = expr( "sin(e1)+mu*e1:e1:mu" );
    auto expanded = base.applySymbolsExpr( symbolsExpr( symbolExpr( "e1", e1 ),
                                                        symbolExpr( "mu", cst( 3.0 ) ) ) );
    expanded.setParameterValues( { { "u", 2.0 } } );

    BOOST_CHECK( expanded.hasSymbolDependency( "u" ) );

    auto expandedDiff = expanded.diff<1>( "u" );
    double const expected = 4.0*( std::cos( 4.0 ) + 3.0 );

    BOOST_CHECK_SMALL( std::abs( scalarValue( expandedDiff ) - expected ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( vector_component_suffixes_are_differentiable_symbols )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto vectorSource = expr<2,1>( "{a,b}:a:b" );
    auto base = expr( "u_0*u_0 + 3*u_1:u_0:u_1" );
    auto expanded = expr( base, symbolExpr( "u", vectorSource, SymbolExprComponentSuffix( 2, 1 ) ) );

    expanded.setParameterValues( { { "a", 2.0 }, { "b", 4.0 } } );

    auto expandedDa = expanded.diff<1>( "a" );
    auto expandedDb = expanded.diff<1>( "b" );

    BOOST_CHECK_SMALL( std::abs( scalarValue( expandedDa ) - 4.0 ), 1e-12 );
    BOOST_CHECK_SMALL( std::abs( scalarValue( expandedDb ) - 3.0 ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( concepts_and_diagnostic_reports_describe_symbolic_ad_capabilities )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto law = expr( "mu*(u*u+1):mu:u" );
    law.setParameterValues( { { "mu", 2.0 }, { "u", 3.0 } } );

    static_assert( Feel::VfExpr<decltype( law )> );
    static_assert( Feel::HasSymbolDependency<decltype( law )> );
    static_assert( Feel::ParametricSymbolExpr<decltype( law )> );
    static_assert( Feel::SymbolicallyDifferentiableExpr<decltype( law )> );

    auto se = symbolsExpr( symbolExpr( "u", cst( 3.0 ) ) );
    static_assert( Feel::AppliesSymbolExpr<decltype( law ), decltype( se )> );

    auto dependencies = describeSymbolicDependencies( law, { "u", "mu", "unused" } );
    BOOST_CHECK( dependencies.supportsDependencyQuery );
    BOOST_CHECK( dependencies.hasDependency() );
    BOOST_CHECK( dependencies.dependsOn( "u" ) );
    BOOST_CHECK( dependencies.dependsOn( "mu" ) );
    BOOST_CHECK( !dependencies.dependsOn( "unused" ) );
    BOOST_CHECK( contains( dependencies.summary(), "dependentSymbols" ) );

    auto lawDiff = describeDifferentiability( law, "u" );
    BOOST_CHECK( lawDiff.ok() );
    BOOST_CHECK( lawDiff.supportsDependencyQuery );
    BOOST_CHECK( lawDiff.hasDependency );
    BOOST_CHECK( lawDiff.supportsSymbolicDiff );
    BOOST_CHECK( lawDiff.kind == SymbolicDifferentiationDiagnosticKind::None );

    auto nativeUnsupported = abs( law );
    auto nativeUnsupportedDiff = describeDifferentiability( nativeUnsupported, "u" );
    BOOST_CHECK( !nativeUnsupportedDiff.ok() );
    BOOST_CHECK( nativeUnsupportedDiff.supportsDependencyQuery );
    BOOST_CHECK( nativeUnsupportedDiff.hasDependency );
    BOOST_CHECK( nativeUnsupportedDiff.supportsSymbolicDiff );
    BOOST_CHECK( nativeUnsupportedDiff.kind == SymbolicDifferentiationDiagnosticKind::UnsupportedPiecewiseFunction );
    BOOST_CHECK( contains( nativeUnsupportedDiff.message, "native unary function" ) );
    BOOST_CHECK( contains( nativeUnsupportedDiff.summary(), "unsupported piecewise function" ) );
}

BOOST_AUTO_TEST_CASE( coordinate_derivatives_propagate_through_fe_field_substitutions )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = loadMesh( _mesh = new Mesh<Simplex<2,1>> );
    auto Xh = Pch<2>( mesh );
    auto uh = Xh->element( Px() + 2.0*Py() );

    auto base = expr( "u*u:u" );
    auto expanded = expr( base, symbolExpr( "u", idv( uh ) ) );

    auto diffX = expanded.diff<1>( "x" );
    auto diffY = expanded.diff<1>( "y" );
    auto gradExpanded = grad<2>( expanded );

    auto exactX = 2.0*( Px() + 2.0*Py() );
    auto exactY = 4.0*( Px() + 2.0*Py() );
    auto exactGrad = trans( vec( exactX, exactY ) );

    double const errorX = normL2( _range = elements( mesh ), _expr = diffX - exactX );
    double const errorY = normL2( _range = elements( mesh ), _expr = diffY - exactY );
    double const errorGrad = normL2( _range = elements( mesh ), _expr = gradExpanded - exactGrad );

    BOOST_CHECK_SMALL( errorX, 1e-10 );
    BOOST_CHECK_SMALL( errorY, 1e-10 );
    BOOST_CHECK_SMALL( errorGrad, 1e-10 );
}

BOOST_AUTO_TEST_CASE( native_math_functor_chain_rules_match_analytic_values )
{
    using namespace Feel;
    using namespace Feel::vf;

    double const x = 0.4;
    auto u = expr( "u:u" );
    u.setParameterValues( { { "u", x } } );

    checkDerivativeValue( "sin", sin( u ), "u", std::cos( x ) );
    checkDerivativeValue( "cos", cos( u ), "u", -std::sin( x ) );
    checkDerivativeValue( "tan", tan( u ), "u", 1.0/( std::cos( x )*std::cos( x ) ) );
    checkDerivativeValue( "acos", acos( u ), "u", -1.0/std::sqrt( 1.0 - x*x ) );
    checkDerivativeValue( "asin", asin( u ), "u", 1.0/std::sqrt( 1.0 - x*x ) );
    checkDerivativeValue( "atan", atan( u ), "u", 1.0/( 1.0 + x*x ) );
    checkDerivativeValue( "cosh", cosh( u ), "u", std::sinh( x ) );
    checkDerivativeValue( "sinh", sinh( u ), "u", std::cosh( x ) );
    checkDerivativeValue( "tanh", tanh( u ), "u", 1.0/( std::cosh( x )*std::cosh( x ) ) );
    checkDerivativeValue( "exp", exp( u ), "u", std::exp( x ) );
    checkDerivativeValue( "loge", loge( u ), "u", 1.0/x );
    checkDerivativeValue( "sqrt", sqrt( u ), "u", 1.0/( 2.0*std::sqrt( x ) ) );

    double const y = 0.8;
    auto v = expr( "v:v" );
    v.setParameterValues( { { "v", y } } );

    auto angle = atan2( u, v );
    double const denominator = x*x + y*y;
    checkDerivativeValue( "atan2 first argument", angle, "u", y/denominator );
    checkDerivativeValue( "atan2 second argument", angle, "v", -x/denominator );
}

BOOST_AUTO_TEST_CASE( power_chain_rules_and_domain_policy_are_explicit )
{
    using namespace Feel;
    using namespace Feel::vf;

    double const x = 1.4;
    double const y = 2.3;

    auto u = expr( "u:u" );
    auto g = expr( "g:g" );
    std::map<std::string,double> const parameters = { { "u", x }, { "g", y } };

    auto cubic = pow( u, 3.0 );
    cubic.setParameterValues( parameters );
    checkDerivativeValue( "pow constant exponent", cubic, "u", 3.0*x*x );
    BOOST_CHECK_SMALL( std::abs( scalarValue( cubic.diff<1>( "u" ) ) -
                                 finiteDifference( cubic, "u", x, parameters ) ), 1e-8 );

    auto variablePowerWithoutDomain = pow( u, g );
    variablePowerWithoutDomain.setParameterValues( parameters );
    checkDerivativeValue( "pow symbol-independent exponent", variablePowerWithoutDomain, "u",
                          y*std::pow( x, y - 1.0 ) );

    auto missingDomain = describeDifferentiability( variablePowerWithoutDomain, "g" );
    BOOST_CHECK( !missingDomain.ok() );
    BOOST_CHECK( missingDomain.kind == SymbolicDifferentiationDiagnosticKind::MissingDomainAssumption );
    BOOST_CHECK( contains( missingDomain.message, "positive" ) );

    auto positivePower = powPositiveBase( u, g );
    positivePower.setParameterValues( parameters );
    double const xPowerY = std::pow( x, y );
    checkDerivativeValue( "pow positive base first argument", positivePower, "u", y*std::pow( x, y - 1.0 ) );
    checkDerivativeValue( "pow positive base exponent", positivePower, "g", xPowerY*std::log( x ) );
    BOOST_CHECK_SMALL( std::abs( scalarValue( positivePower.diff<1>( "u" ) ) -
                                 finiteDifference( positivePower, "u", x, parameters ) ), 1e-8 );
    BOOST_CHECK_SMALL( std::abs( scalarValue( positivePower.diff<1>( "g" ) ) -
                                 finiteDifference( positivePower, "g", y, parameters ) ), 1e-8 );

    auto invalidPositivePower = powPositiveBase( u, g );
    invalidPositivePower.setParameterValues( { { "u", -0.5 }, { "g", 0.5 } } );
    BOOST_CHECK( std::isnan( scalarValue( invalidPositivePower ) ) );
}

BOOST_AUTO_TEST_CASE( unsupported_piecewise_functions_report_clear_diagnostics )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto u = expr( "u:u" );
    u.setParameterValues( { { "u", 0.25 } } );

    auto nativeAbs = abs( u );
    BOOST_CHECK( nativeAbs.hasSymbolDependency( "u" ) );

    try
    {
        auto ignored = nativeAbs.diff<1>( "u" );
        (void)ignored;
        BOOST_FAIL( "native abs differentiation should report an unsupported path" );
    }
    catch ( SymbolicDifferentiationError const& e )
    {
        BOOST_CHECK( e.kind() == SymbolicDifferentiationDiagnosticKind::UnsupportedPiecewiseFunction );
        std::string const message = e.what();
        BOOST_CHECK( contains( message, "native unary function" ) );
        BOOST_CHECK( contains( message, "symbolic differentiation" ) );
        BOOST_CHECK( contains( message, "source:" ) );
    }

    auto nativeFloor = floor( u );
    BOOST_CHECK( nativeFloor.hasSymbolDependency( "u" ) );

    try
    {
        auto ignored = nativeFloor.diff<1>( "u" );
        (void)ignored;
        BOOST_FAIL( "native floor differentiation should report an unsupported path" );
    }
    catch ( SymbolicDifferentiationError const& e )
    {
        BOOST_CHECK( e.kind() == SymbolicDifferentiationDiagnosticKind::UnsupportedPiecewiseFunction );
        std::string const message = e.what();
        BOOST_CHECK( contains( message, "native unary function" ) );
        BOOST_CHECK( contains( message, "symbolic differentiation" ) );
        BOOST_CHECK( contains( message, "source:" ) );
    }

    auto v = expr( "v:v" );
    auto minExpr = min( u, v );
    minExpr.setParameterValues( { { "u", 0.25 }, { "v", 0.5 } } );
    auto minReport = describeDifferentiability( minExpr, "u" );
    BOOST_CHECK( !minReport.ok() );
    BOOST_CHECK( minReport.kind == SymbolicDifferentiationDiagnosticKind::UnsupportedPiecewiseFunction );
    BOOST_CHECK( contains( minReport.message, "min" ) );

    auto maxExpr = max( u, v );
    maxExpr.setParameterValues( { { "u", 0.25 }, { "v", 0.5 } } );
    auto maxReport = describeDifferentiability( maxExpr, "u" );
    BOOST_CHECK( !maxReport.ok() );
    BOOST_CHECK( maxReport.kind == SymbolicDifferentiationDiagnosticKind::UnsupportedPiecewiseFunction );
    BOOST_CHECK( contains( maxReport.message, "max" ) );
}

BOOST_AUTO_TEST_SUITE_END()
