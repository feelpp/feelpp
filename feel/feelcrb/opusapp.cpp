/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-06-08

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
#include <feel/feelcore/feel.hpp>

#include <boost/program_options/options_description.hpp>

namespace Feel
{
/**
 * \return the command lines options of the petsc backend
 */
po::options_description opusapp_options( std::string const& prefix )
{
    std::string _prefix = prefix;

    if ( !_prefix.empty() )
        _prefix += ".";

    po::options_description _options( "Application " + prefix + " options" );
    _options.add_options()
        ( ( _prefix+"run.mode" ).c_str(), Feel::po::value<int>()->default_value( 2 ), "execution mode: pfem=0, scm=1, crb=2" )
        ( ( _prefix+"run.sampling.size" ).c_str(), Feel::po::value<int>()->default_value( 10 ), "size of sampling in parameter space" )
        ( ( _prefix+"run.sampling.mode" ).c_str(), Feel::po::value<int>()->default_value( 10 ), "type of sampling in parameter space (random=0, equidistributed=1)" )
        ( (_prefix+"how-compute-unkown-for-eim").c_str(), Feel::po::value<std::string>()->default_value("CRB-with-ad"), "CRB-with-ad or FEM-with-ad or FEM-without-ad (ad:affine decomposition)")
        ( (_prefix+"export-solution").c_str(), Feel::po::value<bool>()->default_value(false), "export solution and error for each parameters if true")
        ( (_prefix+"run.parameter").c_str(), Feel::po::value<std::vector<double>> ()->multitoken(), "a command line parameter for the online part")
    ;

    return _options;
}

}
