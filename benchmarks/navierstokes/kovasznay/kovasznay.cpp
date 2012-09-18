/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-06-18

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file stokes.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-06-18
 */
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description stokesoptions( "Kovasznay Benchmark options" );
    stokesoptions.add_options()
    ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
    ( "mu", Feel::po::value<double>()->default_value( 1.0/40. ), "viscosity coefficient (default from Sherwin/Karnyadakis book)" )
    ( "beta", Feel::po::value<double>()->default_value( 0.0 ), "convection coefficient" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
    ( "weak", "use weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "penalisation", Feel::po::value<double>()->default_value( 1 ), "penalisation parameter for equal order approximation" )
    ( "stab-p", Feel::po::value<bool>()->default_value( true ), "0 = no stabilisation for equal order approx., 1 = stabilisation for equal order approx." )
    ( "stab-div", Feel::po::value<bool>()->default_value( false ), "0 = no stabilisation for divergence, 1 = stabilisation of divergence." )

    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return stokesoptions.add( Feel::feel_options() ) ;
}


/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "kovasznay" ,
                           "kovasznay" ,
                           "0.2",
                           "Kovasznay benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2009-2010 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}


