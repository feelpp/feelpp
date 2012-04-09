/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2007-05-11

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
   \file kovasznay_osit.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-05-11
 */

#include "kovasznay.hpp"

namespace Feel
{

void
Kovasznay::oseenUpdateIter( Oseen<space_type, imOrder, ENTITY>& oseen,
                            element_U_type& u,
                            bool updateStabilization )
{
    using namespace Feel::vf;

    value_type pi = 4.0 * math::atan( value_type( 1.0 ) );
    value_type lambda = 1./( 2.*M_nu ) - std::sqrt( 1./( 4.*M_nu*M_nu ) + 4.*pi*pi );
    AUTO( uxe, 1. - exp( lambda * Px() ) * cos( 2.*pi*Py() ) );
    AUTO( uye, lambda/( 2.*pi ) * exp( lambda * Px() ) * sin( 2.*pi*Py() ) );

    oseen.update( /* itRan = */ marked2elements( *u.functionSpace()->mesh(), 1 ),
                                /* sigma = */ 0.0,
                                /* nuInc = */ 0.0,
                                /* nuAbs = */ constant( M_nu ),
                                /* beta  = uxe*oneX() + uye*oneY(), */
                                /* beta  = */ idv( u ),
                                /* f     = */ 0.0*oneX(),
                                /* c     = */ 0.0,
                                /* g     = */ uxe*oneX() + uye*oneY(),
                                /* noSlip= */ 1.0,
                                /* updtJ = */ updateStabilization );

}

} // Feel
