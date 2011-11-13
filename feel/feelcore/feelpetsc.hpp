/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2011-11-13
 */
#ifndef __FeelPETSc_H
#define __FeelPETSc_H 1

#include <feel/feelcore/feel.hpp>

#if defined( HAVE_PETSC_H )
extern "C"
{
#include <petsc.h>
#include <petscerror.h>
}
#if defined( FEEL_HAVE_SLEPC )
# include <slepc.h>
#endif /* HAVE_SLEPC */

#endif /* HAVE_PETSC_H */

#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
typedef  PetscBool PetscTruth;
#else
typedef PetscTruth  PetscBool;
#endif

namespace Feel
{
namespace PETSc
{
FEEL_STRONG_INLINE int VecDestroy( Vec& vec )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ::VecDestroy( &vec );
#else
    ::VecDestroy( vec );
#endif
}
FEEL_STRONG_INLINE int VecScatterDestroy( VecScatter& scatter )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ::VecScatterDestroy(&scatter);
#else
    ::VecScatterDestroy(scatter);
#endif


}
FEEL_STRONG_INLINE int MatDestroy( Mat& mat )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ::MatDestroy( &mat );
#else
    ::MatDestroy( mat );
#endif
}

FEEL_STRONG_INLINE int ISDestroy( IS& is )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ::ISDestroy( &is );
#else
    ::ISDestroy( is );
#endif
}
FEEL_STRONG_INLINE int KSPDestroy (KSP& ksp )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ::KSPDestroy( &ksp );
#else
    ::KSPDestroy( ksp );
#endif
}
FEEL_STRONG_INLINE int PCDestroy (PC& pc )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ::PCDestroy( &pc);
#else
    ::PCDestroy( pc );
#endif
}

FEEL_STRONG_INLINE int MatNullSpaceDestroy(MatNullSpace& nullsp )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ::MatNullSpaceDestroy( &nullsp );
#else
    ::MatNullSpaceDestroy( nullsp );
#endif
}
FEEL_STRONG_INLINE int SNESDestroy (SNES& snes )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ::SNESDestroy( &snes);
#else
    ::SNESDestroy( snes );
#endif
}

FEEL_STRONG_INLINE int PetscViewerDestroy (PetscViewer& petsc_viewer)
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    ::PetscViewerDestroy( &petsc_viewer );
#else
    ::PetscViewerDestroy( petsc_viewer );
#endif
}


} // PETSc
} // Feel

#endif /* __FeelPETSc_H */
