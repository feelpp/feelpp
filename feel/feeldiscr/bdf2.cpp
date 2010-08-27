/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-06-10

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file bdf2.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-06-10
 */
#include <boost/timer.hpp>

#include <feel/feeldiscr/bdf2.hpp>

namespace Feel
{

/**
 * \return the command lines options for BDF
 */
po::options_description bdf_options( std::string const& prefix )
{
    std::string _prefix = prefix;
    if ( !_prefix.empty() )
        _prefix += "-";

    po::options_description _options( "BDF (Backward Differences time discretization) options (" + prefix + ")");
    _options.add_options()
        // solver options
        ((_prefix+"bdf-time-initial").c_str(), Feel::po::value<double>()->default_value( 0.0 ), "initial time")
        ((_prefix+"bdf-time-final").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "final time")
        ((_prefix+"bdf-time-step").c_str(), Feel::po::value<double>()->default_value( 1.0 ), "time step")
        ((_prefix+"bdf-time-order").c_str(), Feel::po::value<int>()->default_value( 1 ), "order in time")
        ((_prefix+"bdf-time-strategy").c_str(), Feel::po::value<int>()->default_value( 0 ), "strategy, 0=constant time steps, 1=adaptive time steps")
        ;
    return _options;
}

}

