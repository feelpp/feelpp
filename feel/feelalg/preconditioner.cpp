/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-01-16

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file preconditioner.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-01-16
 */
#define FEELPP_INSTANTIATE_PRECONDITIONER 1

#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelalg/preconditionerpetsc.hpp>

namespace Feel
{
template <typename T>
typename Preconditioner<T>::preconditioner_ptrtype
Preconditioner<T>::build( std::string const& name,
                          BackendType backend,
                          WorldComm const& worldComm )
{
    switch ( backend )
    {
    default:
    case BACKEND_PETSC:
    {
        return preconditioner_ptrtype( new PreconditionerPetsc<T>( name, worldComm ) );
    }
    }
    return preconditioner_ptrtype();
}

template <>
typename Preconditioner<std::complex<double>>::preconditioner_ptrtype
Preconditioner<std::complex<double>>::build( std::string const& name,
                                             BackendType backend,
                                             WorldComm const& worldComm )
{
#if defined(PETSC_HAS_COMPLEX_SUPPORT ) 
    switch ( backend )
    {
    default:
    case BACKEND_PETSC:
        {
            return preconditioner_ptrtype( new PreconditionerPetsc<T>( name, worldComm ) );
        }
    }
#endif
    return preconditioner_ptrtype();
}
    
template <typename T>
FEELPP_STRONG_INLINE
void
Preconditioner<T>::setMatrix( sparse_matrix_ptrtype mat )
{
    if (M_is_initialized)
    {
        M_mat_has_changed = true;
        //this->clear();
    }

    //M_is_initialized = false;
    M_matrix = mat;
}

template <typename T>
void
Preconditioner<T>::setType ( const PreconditionerType pct )
{
    if (M_is_initialized && M_preconditioner_type!=pct)
    {
        this->clear();
    }

    M_preconditioner_type = pct;

}

template <typename T>
void
Preconditioner<T>::setMatSolverPackageType ( const MatSolverPackageType mspt )
{
    if (M_is_initialized && M_matSolverPackage_type!=mspt )
    {
        this->clear();
    }

    M_matSolverPackage_type  = mspt;
    //M_is_initialized = false;
}

template <typename T>
void
Preconditioner<T>::setPrecMatrixStructure( MatrixStructure mstruct  )
{
    M_prec_matrix_structure = mstruct;
}



template class Preconditioner<double>;
template class Preconditioner<std::complex<double>>;

} // Feel
