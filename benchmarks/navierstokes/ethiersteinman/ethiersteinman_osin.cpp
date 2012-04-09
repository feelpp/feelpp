/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2007-04-08

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
   \file ethiersteinman_osin.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-04-08
 */

#include "ethiersteinman.hpp"

namespace Feel
{

void
EthierSteinman::oseenUpdateInit( Oseen<space_u_type, space_p_type, imOrder, ENTITY>& oseen,
                                 mesh_ptr_type mesh,
                                 value_type dt )
{
    using namespace Feel::vf;

    oseen.update( /* itRan = */ elements( *mesh ),
                                /* sigma = */ 1.0/dt,
                                /* nuInc = */ M_mu,
                                /* nuAbs = */ 0.0,
                                /* beta  = */ oneX()-oneX(),
                                /* f     = */ oneX()-oneX(),
                                /* c     = */ 0.0,
                                /* g     = */ oneX()-oneX(),
                                /* noSlip= */ 0.0,
                                /* updtJ = */ false );
}

} // Feel
