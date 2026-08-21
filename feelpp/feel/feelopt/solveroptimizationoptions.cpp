/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <feel/feelcore/environment.hpp>
#include <feel/feelopt/solveroptimizationoptions.hpp>

namespace Feel
{

po::options_description
solveroptimization_options( std::string const& prefix )
{
    po::options_description options(
        "Optimization solver " + prefix + " options" );
    options.add_options()
        ( prefixvm( prefix, "tao-type" ).c_str(),
          po::value<std::string>()->default_value( "lmvm" ),
          "optimization algorithm (lmvm, blmvm, cg, nls, ntr, ...)" )
        ( prefixvm( prefix, "tao-maxit" ).c_str(),
          po::value<int>()->default_value( 1000 ),
          "maximum number of optimization iterations" )
        ( prefixvm( prefix, "tao-maxfcn" ).c_str(),
          po::value<int>()->default_value( 10000 ),
          "maximum number of objective evaluations" )
        ( prefixvm( prefix, "tao-gatol" ).c_str(),
          po::value<double>()->default_value( 1e-8 ),
          "absolute gradient tolerance" )
        ( prefixvm( prefix, "tao-grtol" ).c_str(),
          po::value<double>()->default_value( 1e-8 ),
          "objective-relative gradient tolerance" )
        ( prefixvm( prefix, "tao-gttol" ).c_str(),
          po::value<double>()->default_value( 0.0 ),
          "initial-gradient-relative tolerance" )
        ( prefixvm( prefix, "tao-steptol" ).c_str(),
          po::value<double>()->default_value( 0.0 ),
          "step or trust-region-radius tolerance" )
        ( prefixvm( prefix, "tao-monitor" ).c_str(),
          po::value<bool>()->default_value( false ),
          "print the optimization iteration monitor" )
        ( prefixvm( prefix, "tao-converged-reason" ).c_str(),
          po::value<bool>()->default_value( false ),
          "print the optimization convergence reason" )
        ( prefixvm( prefix, "tao-view" ).c_str(),
          po::value<bool>()->default_value( false ),
          "print the configured optimization solver" );
    return options;
}

} // namespace Feel
