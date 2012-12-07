/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-18

  Copyright (C) 2005,2006,2009 EPFL

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
   \file application.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-18
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>

#if defined(FEELPP_HAS_TAU)
#include <Profile/Profiler.h>
#endif /* FEELPP_HAS_TAU */


namespace Feel
{
FEELPP_NO_EXPORT
po::options_description
serialOptions()
{
    po::options_description serial( "Serial application options" );
    return serial;
}

Application::Application( int argc,
                          char** argv,
                          AboutData const& ad )
    :
    super( argc, argv, ad, serialOptions(), true )
{

#if defined(FEELPP_HAS_TAU)
    TAU_PROFILE_SET_NODE( _S_process_id );
#endif /* FEELPP_HAS_TAU */
}


Application::Application( int argc,
                          char** argv,
                          AboutData const& ad,
                          po::options_description const& od )
    :
    super( argc, argv, ad, serialOptions().add( od ), true )

{
#if defined(FEELPP_HAS_TAU)
    TAU_PROFILE_SET_NODE( _S_process_id );
#endif /* FEELPP_HAS_TAU */
}

Application::Application( Application const& a )
    :
    super( a )
{
}

Application::~Application()
{
}

}
