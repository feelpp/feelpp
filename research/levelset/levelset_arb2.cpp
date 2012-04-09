/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2007-06-19

  Copyright (C) 2007 EPFL

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
   \file levelset_arb2.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-06-19
 */

#include "levelset.hpp"

namespace Feel
{

void
LevelSet::advReactUpdateBdf2( AdvReact<space_p_type, imOrder, ENTITY>& advReact,
                              double dt,
                              double sign,
                              const element_p_type& vx,
                              const element_p_type& vy,
                              const element_p_type& phi,
                              const element_p_type& phio,
                              bool updateStabilization )
{
    using namespace Feel::vf;

    AUTO( beta, idv( vx )*oneX()+idv( vy )*oneY() );

    // coefficients for BDF2
    const double bdf2_0 =  1.5/dt;
    const double bdf2_1 =  2.0/dt;
    const double bdf2_2 = -0.5/dt;

    advReact.update( /* sigma = */ bdf2_0,
                                   /* beta  = */ sign*beta,
                                   /* f     = */ ( idv( phi )*bdf2_1 +
                                           idv( phio )*bdf2_2 ),
                                   /* g     = */ ( idv( phi )
                                           - sign*dt * ( gradv( phi )*( beta ) ) ),
                                   /* updtJ = */ updateStabilization
                   );
}

} // Feel
