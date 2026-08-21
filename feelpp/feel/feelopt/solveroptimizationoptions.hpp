/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#ifndef FEELPP_FEELOPT_SOLVEROPTIMIZATIONOPTIONS_HPP
#define FEELPP_FEELOPT_SOLVEROPTIMIZATIONOPTIONS_HPP 1

#include <string>

#include <boost/program_options/options_description.hpp>

namespace Feel
{

namespace po = boost::program_options;

/**
 * Define backend-independent optimization solver options.
 *
 * The @p prefix permits independent TAO-backed solver instances while the
 * option names remain part of the Feel++ command-line interface.
 *
 * @param prefix solver-specific option prefix
 * @return optimization option description
 */
po::options_description
solveroptimization_options( std::string const& prefix = "" );

} // namespace Feel

#endif // FEELPP_FEELOPT_SOLVEROPTIMIZATIONOPTIONS_HPP
