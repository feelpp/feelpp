/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <cmath>

#include <feel/feelcore/environment.hpp>
#include <feel/feelopt/solveroptimizationfield.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

/**
 * Minimize a variational field-tracking functional with PETSc TAO.
 *
 * The example illustrates the field-native optimization interface: the
 * unknown is a Feel++ element, while gradient and Hessian callbacks assemble
 * ordinary `form1` and `form2` expressions through lightweight proxies.
 *
 * @param argc command-line argument count
 * @param argv command-line argument values
 * @return zero when TAO converges to the expected field, nonzero otherwise
 */
int main( int argc, char** argv )
{
    using namespace Feel;
    using Feel::cout;

    Environment environment(
        _argc = argc, _argv = argv,
        _about = about( _name = "qs_optimization_field",
                        _author = "Feel++ Consortium",
                        _email = "feelpp-devel@feelpp.org" ) );

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto control = Xh->element( "control" );
    auto target = Xh->element( expr( "1+x+y:x:y" ), "target" );
    auto trial = Xh->element( "u" );
    auto test = Xh->element( "v" );

    auto optimization = optimizationField(
        _space = Xh, _name = "field-tracking", _backend = "petsc" );
    optimization
        .algorithm( "nls" )
        .maxIterations( 20 )
        .gradientTolerance( 1e-10 )
        .objectiveGradient(
            [mesh, target, test]( auto const& field, auto& gradient )
            {
                auto const mismatch = idv( field ) - idv( target );
                gradient = integrate(
                    _range = elements( mesh ),
                    _expr = mismatch * id( test ) );
                return 0.5 * integrate(
                                 _range = elements( mesh ),
                                 _expr = mismatch * mismatch )
                                 .evaluate()( 0, 0 );
            } )
        .hessian(
            [mesh, trial, test]( auto const&, auto& hessian )
            {
                hessian = integrate(
                    _range = elements( mesh ),
                    _expr = idt( trial ) * id( test ) );
            } );

    auto const result = optimization.solve( control );
    auto const trackingError = normL2(
        _range = elements( mesh ), _expr = idv( control ) - idv( target ) );

    double constexpr objectiveTolerance = 1e-12;
    double constexpr gradientTolerance = 1e-8;
    double constexpr trackingTolerance = 1e-6;
    bool const resultIsCorrect =
        result.converged && std::isfinite( result.objective ) &&
        std::abs( result.objective ) < objectiveTolerance &&
        std::isfinite( result.gradientNorm ) &&
        result.gradientNorm < gradientTolerance &&
        std::isfinite( trackingError ) && trackingError < trackingTolerance;

    cout << "TAO " << result.solverType << ": " << result.reason
         << ", iterations=" << result.iterations
         << ", objective=" << result.objective
         << ", gradient residual=" << result.gradientNorm
         << ", ||control-target||_L2=" << trackingError << '\n';

    auto exporter = Feel::exporter( _mesh = mesh,
                                    _name = "qs_optimization_field" );
    exporter->add( "target", target );
    exporter->add( "control", control );
    exporter->save();

    return resultIsCorrect ? 0 : 1;
}
