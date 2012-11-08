/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-06-06

  Copyright (C) 2007-2011 Université Joseph Fourier (Grenoble I)

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
   \file options.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-06-06
 */
#include <options.hpp>
#include <feel/feelalg/backendgmm.hpp>
#include <feel/feelalg/backend_adaptive_reuse_pc.hpp>
#include <feel/feelalg/backendpetsc.hpp>
#include <feel/feeldiscr/oseen.hpp>

Feel::po::options_description
makeOptions()
{
    Feel::po::options_description kovasznayoptions( "Kovasznay options" );
    kovasznayoptions.add_options()
    ( "nu", Feel::po::value<double>()->default_value( 1.0 ),
      "viscosity" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ),
      "first h value to start convergence" )
    ( "fixpointtol", Feel::po::value<double>()->default_value( 0.025 ),
      "Convergence tolerance for fixed point sub-iterations" )
    ( "maxsubiter", Feel::po::value<int>()->default_value( 1 ),
      "maximal number of fixed point sub-iterations" )
    ( "stabtheta", Feel::po::value<double>()->default_value( 1.5 ),
      "stabilization decoupling coefficient" )
    ( "doexport", Feel::po::value<int>()->default_value( 0 ),
      "stride for result export (0=no export)" )
    ;

    // modification of defaults from oseen and backends
    // to defaults more sensible for Kovasznay
    Feel::OseenDefaults oseenDefs;
    oseenDefs.STAB_COEFF_P = 0.01; // P2P2
    Feel::BackendGmmDefaults beGmmOsDefs;
    //     beGmmOsDefs.solver_type = "gmres";
    beGmmOsDefs.fillin = 20;
    beGmmOsDefs.threshold = 1.e-4;
    Feel::BackendGmmDefaults beGmmSyDefs;
    //     beGmmOsDefs.solver_type = "cg";
    beGmmSyDefs.pc_type = "id";

    return kovasznayoptions
           .add( Feel::oseen_options( std::string(), oseenDefs ) )
           .add( Feel::backend_adaptive_reuse_pc_options( "oseen" ) )
           .add( Feel::backendgmm_options( "oseen", beGmmOsDefs ) )
           .add( Feel::backend_adaptive_reuse_pc_options( "symm" ) )
           .add( Feel::backendgmm_options( "symm", beGmmSyDefs ) );
}


Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "kovasznay" ,
                           "kovasznay" ,
                           "0.1",
                           "Kovasznay benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2011 Université Joseph Fourier" );

    about.addAuthor( "Christoph Winkelmann", "developer",
                     "christoph.winkelmann@epfl.ch", "" );
    about.addAuthor( "Christophe Prud'homme", "developer",
                     "christophe.prudhomme@feelpp.org", "" );
    return about;

}

