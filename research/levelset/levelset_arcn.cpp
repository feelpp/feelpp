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
   \file levelset_arcn.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-06-19
 */

#include "levelset.hpp"

namespace Feel
{

void
LevelSet::advReactUpdateCN( AdvReact<space_p_type, imOrder, ENTITY>& advReact,
                            double dt,
                            double theta,
                            double sign,
                            const element_p_type& vx,
                            const element_p_type& vy,
                            const element_p_type& phi,
                            bool updateStabilization )
{
    using namespace Feel::vf;

    AUTO( beta, idv( vx )*oneX()+idv( vy )*oneY() );

    advReact.update( /* sigma = */ 1.0/dt,
                                   /* beta  = */ theta*sign*beta,
                                   /* f     = */ ( idv( phi )/dt
                                           - ( 1.0-theta )*sign
                                           * ( gradv( phi )*( beta ) ) ),
                                   /* g     = */ ( idv( phi )
                                           - sign*dt * ( gradv( phi )*( beta ) ) ),
                                   /* updtJ = */ updateStabilization
                   );
}


} // Feel
