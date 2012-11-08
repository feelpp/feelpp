/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-05-23

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
/**
   \file oseen.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-05-23
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/oseendata.hpp>

namespace Feel
{

/**
 * \return the command lines options of the oseen solver
 */
po::options_description oseen_options( std::string const& prefix,
                                       OseenDefaults defaults )
{
    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += "-";

    po::options_description _oseen( "Oseen " + prefix + " options" );
    _oseen.add_options()
    ( ( _prefix+"oseen-bc-coeff-diff" ).c_str(),
      po::value<double>()->default_value( defaults.BC_COEFF_DIFF ),
      "coefficient for diffusive terms of weak Dirichlet conditions" )

    ( ( _prefix+"oseen-bc-coeff-conv" ).c_str(),
      po::value<double>()->default_value( defaults.BC_COEFF_CONV ),
      "coefficient for convective terms of weak Dirichlet conditions" )

    ( ( _prefix+"oseen-stab-coeff-div" ).c_str(),
      po::value<double>()->default_value( defaults.STAB_COEFF_DIV ),
      "coefficient for convective terms of weak Dirichlet conditions" )

    ( ( _prefix+"oseen-stab-coeff-p" ).c_str(),
      po::value<double>()->default_value( defaults.STAB_COEFF_P ),
      "coefficient for convective terms of weak Dirichlet conditions" )

    ( ( _prefix+"oseen-eps-compress" ).c_str(),
      po::value<double>()->default_value( defaults.EPS_COMPRESS ),
      "coefficient for convective terms of weak Dirichlet conditions" )

    ( ( _prefix+"oseen-divdiv-coeff" ).c_str(),
      po::value<double>()->default_value( defaults.DIVDIV_COEFF ),
      "coefficient for convective terms of weak Dirichlet conditions" )

    ( ( _prefix+"oseen-weak-dirichlet" ).c_str(),
      po::value<bool>()->default_value( defaults.WEAK_DIRICHLET ),
      "use weak Dirichlet BC imposition (0=strong, 1=weak)" )

    ( ( _prefix+"oseen-export-matlab" ).c_str(),
      po::value<bool>()->default_value( defaults.EXPORT_MATLAB ),
      "export left and right hand side in matlab format" )

    ;
    return _oseen;
}

} // Feel
