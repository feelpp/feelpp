/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-11-13

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
   \file feelslepc.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-11-13
 */
#ifndef __FeelSLEPc_H
#define __FeelSLEPc_H 1

#include <feel/feelcore/feelpetsc.hpp>

#if defined( FEELPP_HAS_SLEPC )
# include <slepc.h>
#endif /* FEELPP_HAS_SLEPC */



namespace Feel
{
namespace SLEPc
{
FEELPP_STRONG_INLINE int EPSDestroy ( EPS& eps )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    return ::EPSDestroy( &eps );
#else
    return ::EPSDestroy( eps );
#endif
}

} // slepc
} // Feel

#endif /* __FeelSLEPc_H */

