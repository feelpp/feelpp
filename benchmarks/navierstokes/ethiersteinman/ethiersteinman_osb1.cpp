/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2007-04-08

  Copyright (C) 2007 EPFL

  This library is free software; you can redistribute it and/or
x  modify it under the terms of the GNU Lesser General Public
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
   \file ethiersteinman_osb1.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-04-08
 */

#include "ethiersteinman.hpp"

namespace Feel
{

void
EthierSteinman::oseenUpdateBdf1( Oseen<space_u_type, space_p_type, imOrder, ENTITY>& oseen,
                                 value_type dt,
                                 element_u_type& uxn,
                                 element_u_type& uyn,
                                 element_u_type& uzn,
                                 element_u_type& ux,
                                 element_u_type& uy,
                                 element_u_type& uz,
                                 element_u_type& ux0,
                                 element_u_type& uy0,
                                 element_u_type& uz0,
                                 bool updateStabilization,
                                 element_p_type& pn )
{
    using namespace Feel::vf;

    value_type epsCompress = this->vm()["epscompress"].as<double>();

    oseen.update( /* itRan = */ marked2elements( *uxn.functionSpace()->mesh(),
                                1 ),
                                /* sigma = */ 0.0,
                                /* nuInc = */ 0.0,
                                /* nuAbs = */ M_mu,
                                /* beta  = */ ( idv( uxn )*oneX()
                                        + idv( uyn )*oneY()
                                        //                                   + idv(uzn)*oneZ()
                                              ),
                                /* f     = */ ( ( idv( ux )/dt )*oneX()
                                        + ( idv( uy )/dt )*oneY()
                                        //                                   + (idv(uz)/dt)*oneZ()
                                              ),
                                /* c     = */ epsCompress*idv( pn ),
                                /* g     = */ ( idv( ux0 )*oneX()
                                        + idv( uy0 )*oneY()
                                        //                                   + idv(uz0)*oneZ()
                                              ),
                                /* noSlip= */ 1.0,
                                /* updtJ = */ updateStabilization );
}

} // Feel
