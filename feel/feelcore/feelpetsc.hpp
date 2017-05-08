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
   \file feelpetsc.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-11-13
 */
#ifndef __FeelPETSc_H
#define __FeelPETSc_H 1

#include <feel/feelcore/feel.hpp>

#if defined( FEELPP_HAS_PETSC_H )
#include <petsc.h>
#include <petscerror.h>
#if defined( FEELPP_HAS_SLEPC )
# include <slepc.h>
#endif /* FEELPP_HAS_SLEPC */

#endif /* FEELPP_HAS_PETSC_H */

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
typedef  PetscBool PetscTruth;
#else
typedef PetscTruth  PetscBool;
#endif

#if !defined(PETSC_VERSION_LT) || !defined(PETSC_VERSION_GT) || !defined(PETSC_VERSION_GE)
#define PETSC_VERSION_LESS_THAN(major,minor,subminor)                   \
    ((PETSC_VERSION_MAJOR < (major) ||                                  \
      (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR < (minor) || \
                                          (PETSC_VERSION_MINOR == (minor) && \
                                           PETSC_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)

#define PETSC_VERSION_GREATER_THAN(major,minor,subminor)                \
    ((PETSC_VERSION_MAJOR > (major) ||                                  \
      (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR > (minor) || \
                                          (PETSC_VERSION_MINOR == (minor) && \
                                           PETSC_VERSION_SUBMINOR > (subminor))))) ? 1 : 0)

#define PETSC_VERSION_GREATER_OR_EQUAL_THAN(major,minor,subminor)       \
    ((PETSC_VERSION_MAJOR > (major) ||                                  \
      (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR > (minor) || \
                                          (PETSC_VERSION_MINOR == (minor) && \
                                           PETSC_VERSION_SUBMINOR >= (subminor))))) ? 1 : 0)
#else
#define PETSC_VERSION_LESS_THAN(major,minor,subminor)                   \
    PETSC_VERSION_LT(major,minor,subminor)

#define PETSC_VERSION_GREATER_THAN(major,minor,subminor)                \
    PETSC_VERSION_GT(major,minor,subminor)

#define PETSC_VERSION_GREATER_OR_EQUAL_THAN(major,minor,subminor)       \
    PETSC_VERSION_GE(major,minor,subminor)
#endif

namespace Feel
{
namespace PETSc
{
FEELPP_STRONG_INLINE int VecDestroy( Vec& vec )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    return ::VecDestroy( &vec );
#else
    return ::VecDestroy( vec );
#endif
}
FEELPP_STRONG_INLINE int VecScatterDestroy( VecScatter& scatter )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    return ::VecScatterDestroy( &scatter );
#else
    return ::VecScatterDestroy( scatter );
#endif


}
FEELPP_STRONG_INLINE int MatDestroy( Mat& mat )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    return ::MatDestroy( &mat );
#else
    return ::MatDestroy( mat );
#endif
}

FEELPP_STRONG_INLINE int ISDestroy( IS& is )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    return ::ISDestroy( &is );
#else
    return ::ISDestroy( is );
#endif
}
FEELPP_STRONG_INLINE int KSPDestroy ( KSP& ksp )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    return ::KSPDestroy( &ksp );
#else
    return ::KSPDestroy( ksp );
#endif
}
FEELPP_STRONG_INLINE int PCDestroy ( PC& pc )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    return ::PCDestroy( &pc );
#else
    return ::PCDestroy( pc );
#endif
}

FEELPP_STRONG_INLINE int MatNullSpaceDestroy( MatNullSpace& nullsp )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    return ::MatNullSpaceDestroy( &nullsp );
#else
    return ::MatNullSpaceDestroy( nullsp );
#endif
}
FEELPP_STRONG_INLINE int SNESDestroy ( SNES& snes )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    return ::SNESDestroy( &snes );
#else
    return ::SNESDestroy( snes );
#endif
}

FEELPP_STRONG_INLINE int PetscViewerDestroy ( PetscViewer& petsc_viewer )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    return ::PetscViewerDestroy( &petsc_viewer );
#else
    return ::PetscViewerDestroy( petsc_viewer );
#endif
}


} // PETSc
} // Feel

#endif /* __FeelPETSc_H */
