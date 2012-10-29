/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-09-30

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file pod.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2011-09-30
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelcrb/pod.hpp>


namespace Feel
{
Feel::po::options_description
podOptions( std::string const& prefix )
{
    Feel::po::options_description podoptions( "POD Options" );
    podoptions.add_options()
    ( "pod.store-pod-matrix"   , Feel::po::value<bool>()->default_value( false ), "indicate if we store the pod matrix on a file" )
    ( "pod.store-pod-matrix-format-octave"   , Feel::po::value<bool>()->default_value( false ), "indicate if we store the pod matrix on a file with octave format" )
    ;

    return podoptions;
}
}

