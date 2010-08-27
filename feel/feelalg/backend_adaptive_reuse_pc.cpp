/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-05-23

  Copyright (C) 2007, 2009 Université Joseph Fourier (Grenoble I)

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
   \file backend_adaptive_reuse_pc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-05-23
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/backend_adaptive_reuse_pc.hpp>

namespace Feel
{

/**
 * \return the command lines options of the adaptive reuse pc backend
 */
po::options_description backend_adaptive_reuse_pc_options( std::string const& prefix,
                                                           BackendAdaptiveReusePCdefaults defaults )
{
    std::string _prefix = prefix;
    if ( !_prefix.empty() )
        _prefix += "-";

    po::options_description _options("Backend adaptive reuse pc(arpc) " + prefix + " solver options");
    _options.add_options()
        ((_prefix+"arpc-maxiter").c_str(),
         po::value<int>()->default_value( defaults.maxiter ),
         "maximum number of iterations")
        ;
    return _options;
}

} // Feel


