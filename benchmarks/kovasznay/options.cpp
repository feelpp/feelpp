/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-06

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file options.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-06-06
 */
#include <options.hpp>
#include <life/lifealg/backendgmm.hpp>
#include <life/lifealg/backend_adaptive_reuse_pc.hpp>
#include <life/lifealg/backendpetsc.hpp>
#include <life/lifediscr/oseen.hpp>

Life::po::options_description
makeOptions()
{
    Life::po::options_description kovasznayoptions("Kovasznay options");
    kovasznayoptions.add_options()
        ("nu", Life::po::value<double>()->default_value( 1.0 ),
         "viscosity")
        ("hsize", Life::po::value<double>()->default_value( 0.5 ),
         "first h value to start convergence")
        ("fixpointtol", Life::po::value<double>()->default_value( 0.025 ),
         "Convergence tolerance for fixed point sub-iterations")
        ("maxsubiter", Life::po::value<int>()->default_value( 1 ),
         "maximal number of fixed point sub-iterations")
        ("stabtheta", Life::po::value<double>()->default_value( 1.5 ),
         "stabilization decoupling coefficient")
        ("export", Life::po::value<int>()->default_value( 0 ),
         "stride for result export (0=no export)")
        ;

    // modification of defaults from oseen and backends
    // to defaults more sensible for Kovasznay
    Life::OseenDefaults oseenDefs;
    oseenDefs.STAB_COEFF_P = 0.01; // P2P2
    Life::BackendGmmDefaults beGmmOsDefs;
//     beGmmOsDefs.solver_type = "gmres";
    beGmmOsDefs.fillin = 20;
    beGmmOsDefs.threshold = 1.e-4;
    Life::BackendGmmDefaults beGmmSyDefs;
//     beGmmOsDefs.solver_type = "cg";
    beGmmSyDefs.pc_type = "id";

    return kovasznayoptions
        .add( Life::oseen_options( std::string(), oseenDefs ) )
        .add( Life::backend_adaptive_reuse_pc_options( "oseen" ) )
        .add( Life::backendgmm_options( "oseen", beGmmOsDefs ) )
        .add( Life::backend_adaptive_reuse_pc_options( "symm" ) )
        .add( Life::backendgmm_options( "symm", beGmmSyDefs ) );
}

inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "kovasznay" ,
                           "kovasznay" ,
                           "0.1",
                           "two fluid problem",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2006 EPFL");

    about.addAuthor("Christoph Winkelmann", "developer",
                    "christoph.winkelmann@epfl.ch", "");
    return about;

}

