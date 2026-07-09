/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#define BOOST_TEST_MODULE test_tensorfield_evaluate
#include <feel/feelcore/testsuite.hpp>

#include <cmath>
#include <concepts>
#include <cstdio>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <string>

#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( tensorfield_evaluate )

template<typename MatrixT>
void
check_row( MatrixT const& values, int row, std::initializer_list<double> expected )
{
    BOOST_REQUIRE_EQUAL( values.cols(), expected.size() );
    int col = 0;
    for ( double v : expected )
    {
        auto const actual = values( row, col );
        auto const diff = actual - v;
        BOOST_CHECK_MESSAGE( std::abs( diff ) <= 1e-12,
                             "row=" << row << " col=" << col
                                    << " actual=" << actual
                                    << " expected=" << v
                                    << " diff=" << diff );
        ++col;
    }
}

template<typename MatrixT>
void
check_rows( MatrixT const& values, int npoints, std::initializer_list<double> expected )
{
    BOOST_REQUIRE_EQUAL( values.rows(), npoints );
    for ( int p = 0; p < npoints; ++p )
        check_row( values, p, expected );
}

std::string
read_file( std::string const& filename )
{
    std::ifstream input( filename );
    BOOST_REQUIRE_MESSAGE( input.good(), "cannot open " << filename );
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

void
check_contains( std::string const& text, std::string const& pattern )
{
    BOOST_CHECK_MESSAGE( text.find( pattern ) != std::string::npos,
                         "missing pattern: " << pattern << "\ncontent:\n" << text );
}

template<typename ElementT>
concept HasTensorComponentApi = requires( ElementT& element, ElementT const& constElement )
{
    element.tensorComponent( Component::X, Component::X );
    constElement.tensorComponent( Component::X, Component::X );
};

template<typename ElementT>
concept HasEvaluateTensorApi = requires( ElementT const& element,
                                         typename ElementT::functionspace_type::Context const& context )
{
    element.evaluateTensor( context );
};

template<typename ElementT>
concept HasEvaluateSymmetricApi = requires( ElementT const& element,
                                            typename ElementT::functionspace_type::Context const& context )
{
    element.evaluateSymmetric( context );
    element.evaluateSymmetricStorage( context );
};

template<typename ElementT>
concept HasTensorMatlabDiagnosticApi = requires( ElementT const& element,
                                                 typename ElementT::functionspace_type::Context const& context )
{
    element.printMatlab( std::string{}, context, SymmetricTensorFormat{} );
};

BOOST_AUTO_TEST_CASE( symmetric_tensor2d_formats_are_point_major )
{
    auto mesh = unitSquare();
    auto Xh = Pdhms<1>( mesh );
    auto u = Xh->element();
    static_assert( HasTensorComponentApi<decltype( u )> );
    static_assert( HasEvaluateTensorApi<decltype( u )> );
    static_assert( HasEvaluateSymmetricApi<decltype( u )> );
    static_assert( HasTensorMatlabDiagnosticApi<decltype( u )> );

    auto uxx = u.tensorComponent( Component::X, Component::X );
    auto uxy = u.tensorComponent( Component::X, Component::Y );
    auto uyy = u.tensorComponent( Component::Y, Component::Y );
    uxx.on( _range=elements( mesh ), _expr=cst( 11. ) );
    uxy.on( _range=elements( mesh ), _expr=cst( 12. ) );
    uyy.on( _range=elements( mesh ), _expr=cst( 22. ) );

    auto ctx = Xh->context();
    node_type p0( 2 ), p1( 2 );
    p0( 0 ) = 0.25; p0( 1 ) = 0.25;
    p1( 0 ) = 0.75; p1( 1 ) = 0.75;
    ctx.add( p0 );
    ctx.add( p1 );

    auto xyComponent = u.tensorComponent( Component::X, Component::Y );
    auto yxComponent = u.tensorComponent( Component::Y, Component::X );
    auto componentCtx = xyComponent.functionSpace()->context();
    componentCtx.add( p0 );
    componentCtx.add( p1 );
    auto xy = xyComponent.evaluate( componentCtx );
    auto yx = yxComponent.evaluate( componentCtx );
    BOOST_REQUIRE_EQUAL( xy.size(), ctx.nPoints() );
    for ( int p = 0; p < ctx.nPoints(); ++p )
        BOOST_CHECK_SMALL( xy( p ) - yx( p ), 1e-12 );

    auto storage = u.evaluateSymmetricStorage( ctx );
    check_rows( storage, ctx.nPoints(), { 11., 12., 22. } );

    auto diagonalFirst = u.evaluateSymmetric(
        ctx, SymmetricTensorFormat{ SymmetricTensorOrder::DiagonalFirst,
                                    SymmetricTensorScaling::Tensor } );
    check_rows( diagonalFirst, ctx.nPoints(), { 11., 22., 12. } );
    BOOST_CHECK_SMALL( diagonalFirst( 0, 2 ) - xy( 0 ), 1e-12 );

    auto engineering = u.evaluateSymmetric(
        ctx, SymmetricTensorFormat{ SymmetricTensorOrder::DiagonalFirst,
                                    SymmetricTensorScaling::EngineeringShear } );
    check_rows( engineering, ctx.nPoints(), { 11., 22., 24. } );

    auto mandel = u.evaluateSymmetric(
        ctx, SymmetricTensorFormat{ SymmetricTensorOrder::DiagonalFirst,
                                    SymmetricTensorScaling::Mandel } );
    check_rows( mandel, ctx.nPoints(), { 11., 22., std::sqrt( 2. )*12. } );

    auto vtk = u.evaluateSymmetric(
        ctx, SymmetricTensorFormat{ SymmetricTensorOrder::VtkTensor6,
                                    SymmetricTensorScaling::Tensor } );
    check_rows( vtk, ctx.nPoints(), { 11., 22., 0., 12., 0., 0. } );

    auto xdmf = u.evaluateSymmetric(
        ctx, SymmetricTensorFormat{ SymmetricTensorOrder::XdmfTensor6,
                                    SymmetricTensorScaling::Tensor } );
    check_rows( xdmf, ctx.nPoints(), { 11., 12., 0., 22., 0., 0. } );

    auto full = u.evaluateTensor( ctx );
    check_rows( full, ctx.nPoints(), { 11., 12., 12., 22. } );

    auto raw = u.evaluate( ctx );
    BOOST_CHECK_EQUAL( raw.size(), ctx.nPoints()*decltype( u )::nComponents );
}

BOOST_AUTO_TEST_CASE( matlab_tensor_diagnostic_is_self_describing )
{
    auto mesh = unitSquare();
    auto Xh = Pdhms<1>( mesh );
    auto u = Xh->element();

    auto uxx = u.tensorComponent( Component::X, Component::X );
    auto uxy = u.tensorComponent( Component::X, Component::Y );
    auto uyy = u.tensorComponent( Component::Y, Component::Y );
    uxx.on( _range=elements( mesh ), _expr=cst( 11. ) );
    uxy.on( _range=elements( mesh ), _expr=cst( 12. ) );
    uyy.on( _range=elements( mesh ), _expr=cst( 22. ) );

    auto ctx = Xh->context();
    node_type p( 2 );
    p( 0 ) = 0.5; p( 1 ) = 0.5;
    ctx.add( p );

    std::string const tensorBase = "tensorfield_components_phase4";
    std::string const tensorFile = tensorBase + ".m";
    std::string const rawFile = "tensorfield_raw_phase4.m";
    std::remove( tensorFile.c_str() );
    std::remove( rawFile.c_str() );

    u.printMatlab( tensorBase,
                   ctx,
                   SymmetricTensorFormat{ SymmetricTensorOrder::DiagonalFirst,
                                          SymmetricTensorScaling::EngineeringShear },
                   true,
                   "epsilon_components" );
    u.printMatlab( rawFile );

    if ( Environment::worldComm().isMasterRank() )
    {
        auto const tensorText = read_file( tensorFile );
        check_contains( tensorText, "% Feel++ tensor field diagnostic" );
        check_contains( tensorText, "epsilon_components_order = 'DiagonalFirst';" );
        check_contains( tensorText, "epsilon_components_scaling = 'EngineeringShear';" );
        check_contains( tensorText, "epsilon_components_point_major = true;" );
        check_contains( tensorText, "epsilon_components_npoints = 1;" );
        check_contains( tensorText, "epsilon_components_ncomponents = 3;" );
        check_contains( tensorText, "epsilon_components_labels = {'xx','yy','xy'};" );
        check_contains( tensorText, "1.1000000000000000e+01" );
        check_contains( tensorText, "2.2000000000000000e+01" );
        check_contains( tensorText, "2.4000000000000004e+01" );

        auto const rawText = read_file( rawFile );
        BOOST_CHECK_EQUAL( rawText.find( "epsilon_components_order" ), std::string::npos );
        BOOST_CHECK_EQUAL( rawText.find( "Feel++ tensor field diagnostic" ), std::string::npos );

        std::remove( tensorFile.c_str() );
        std::remove( rawFile.c_str() );
    }
}

BOOST_AUTO_TEST_CASE( symmetric_tensor3d_formats_and_scaling )
{
    auto mesh = unitCube();
    auto Xh = Pdhms<1>( mesh );
    auto u = Xh->element();
    static_assert( HasTensorComponentApi<decltype( u )> );
    static_assert( HasEvaluateTensorApi<decltype( u )> );
    static_assert( HasEvaluateSymmetricApi<decltype( u )> );
    static_assert( HasTensorMatlabDiagnosticApi<decltype( u )> );

    auto uxx = u.tensorComponent( Component::X, Component::X );
    auto uxy = u.tensorComponent( Component::X, Component::Y );
    auto uxz = u.tensorComponent( Component::X, Component::Z );
    auto uyy = u.tensorComponent( Component::Y, Component::Y );
    auto uyz = u.tensorComponent( Component::Y, Component::Z );
    auto uzz = u.tensorComponent( Component::Z, Component::Z );
    uxx.on( _range=elements( mesh ), _expr=cst( 11. ) );
    uxy.on( _range=elements( mesh ), _expr=cst( 12. ) );
    uxz.on( _range=elements( mesh ), _expr=cst( 13. ) );
    uyy.on( _range=elements( mesh ), _expr=cst( 22. ) );
    uyz.on( _range=elements( mesh ), _expr=cst( 23. ) );
    uzz.on( _range=elements( mesh ), _expr=cst( 33. ) );

    auto ctx = Xh->context();
    node_type p0( 3 ), p1( 3 );
    p0( 0 ) = 0.25; p0( 1 ) = 0.25; p0( 2 ) = 0.25;
    p1( 0 ) = 0.75; p1( 1 ) = 0.75; p1( 2 ) = 0.75;
    ctx.add( p0 );
    ctx.add( p1 );

    auto xyComponent = u.tensorComponent( Component::X, Component::Y );
    auto yxComponent = u.tensorComponent( Component::Y, Component::X );
    auto yzComponent = u.tensorComponent( Component::Y, Component::Z );
    auto zyComponent = u.tensorComponent( Component::Z, Component::Y );
    auto componentCtx = xyComponent.functionSpace()->context();
    componentCtx.add( p0 );
    componentCtx.add( p1 );
    auto xy = xyComponent.evaluate( componentCtx );
    auto yx = yxComponent.evaluate( componentCtx );
    auto yz = yzComponent.evaluate( componentCtx );
    auto zy = zyComponent.evaluate( componentCtx );
    for ( int p = 0; p < ctx.nPoints(); ++p )
    {
        BOOST_CHECK_SMALL( xy( p ) - yx( p ), 1e-12 );
        BOOST_CHECK_SMALL( yz( p ) - zy( p ), 1e-12 );
    }

    auto storage = u.evaluateSymmetricStorage( ctx );
    check_rows( storage, ctx.nPoints(), { 11., 12., 13., 22., 23., 33. } );

    auto diagonalFirst = u.evaluateSymmetric(
        ctx, SymmetricTensorFormat{ SymmetricTensorOrder::DiagonalFirst,
                                    SymmetricTensorScaling::Tensor } );
    check_rows( diagonalFirst, ctx.nPoints(), { 11., 22., 33., 12., 13., 23. } );
    BOOST_CHECK_SMALL( diagonalFirst( 0, 3 ) - xy( 0 ), 1e-12 );

    auto engineering = u.evaluateSymmetric(
        ctx, SymmetricTensorFormat{ SymmetricTensorOrder::DiagonalFirst,
                                    SymmetricTensorScaling::EngineeringShear } );
    check_rows( engineering, ctx.nPoints(), { 11., 22., 33., 24., 26., 46. } );

    auto mandel = u.evaluateSymmetric(
        ctx, SymmetricTensorFormat{ SymmetricTensorOrder::DiagonalFirst,
                                    SymmetricTensorScaling::Mandel } );
    check_rows( mandel, ctx.nPoints(), { 11., 22., 33.,
                                         std::sqrt( 2. )*12.,
                                         std::sqrt( 2. )*13.,
                                         std::sqrt( 2. )*23. } );

    auto vtk = u.evaluateSymmetric(
        ctx, SymmetricTensorFormat{ SymmetricTensorOrder::VtkTensor6,
                                    SymmetricTensorScaling::Tensor } );
    check_rows( vtk, ctx.nPoints(), { 11., 22., 33., 12., 23., 13. } );

    auto full = u.evaluateTensor( ctx );
    check_rows( full, ctx.nPoints(), { 11., 12., 13.,
                                       12., 22., 23.,
                                       13., 23., 33. } );
}

BOOST_AUTO_TEST_CASE( full_tensor_fields_evaluate_row_major )
{
    auto mesh = unitSquare();
    auto Xh = Pdhm<1>( mesh );
    auto u = Xh->element( mat<2,2>( cst( 1. ), cst( 2. ),
                                    cst( 3. ), cst( 4. ) ) );
    static_assert( HasTensorComponentApi<decltype( u )> );
    static_assert( HasEvaluateTensorApi<decltype( u )> );
    static_assert( !HasEvaluateSymmetricApi<decltype( u )> );
    static_assert( !HasTensorMatlabDiagnosticApi<decltype( u )> );

    auto ctx = Xh->context();
    node_type p( 2 );
    p( 0 ) = 0.5; p( 1 ) = 0.5;
    ctx.add( p );

    auto full = u.evaluateTensor( ctx );
    check_rows( full, ctx.nPoints(), { 1., 2., 3., 4. } );
}

BOOST_AUTO_TEST_CASE( scalar_and_vector_fields_do_not_expose_tensor_apis )
{
    auto mesh = unitSquare();

    auto Sh = Pdh<1>( mesh );
    auto scalar = Sh->element( cst( 1. ) );
    static_assert( !HasTensorComponentApi<decltype( scalar )> );
    static_assert( !HasEvaluateTensorApi<decltype( scalar )> );
    static_assert( !HasEvaluateSymmetricApi<decltype( scalar )> );
    static_assert( !HasTensorMatlabDiagnosticApi<decltype( scalar )> );
    auto scalarCtx = Sh->context();
    node_type p( 2 );
    p( 0 ) = 0.5; p( 1 ) = 0.5;
    scalarCtx.add( p );
    auto scalarValues = scalar.evaluate( scalarCtx );
    BOOST_REQUIRE_EQUAL( scalarValues.size(), scalarCtx.nPoints() );

    auto Vh = Pdhv<1>( mesh );
    auto vector = Vh->element( vec( cst( 1. ), cst( 2. ) ) );
    static_assert( !HasTensorComponentApi<decltype( vector )> );
    static_assert( !HasEvaluateTensorApi<decltype( vector )> );
    static_assert( !HasEvaluateSymmetricApi<decltype( vector )> );
    static_assert( !HasTensorMatlabDiagnosticApi<decltype( vector )> );
    auto vectorCtx = Vh->context();
    vectorCtx.add( p );
    auto vectorValues = vector.evaluate( vectorCtx );
    BOOST_REQUIRE_EQUAL( vectorValues.size(), vectorCtx.nPoints()*decltype( vector )::nComponents );
}

BOOST_AUTO_TEST_SUITE_END()
