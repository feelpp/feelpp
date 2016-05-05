/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2015-06-09

  Copyright (C) 2015 Université Joseph Fourier (Grenoble I)

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
   \file preconditionerpetscfeelpp.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2015-06-09
 */

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3, 3, 0 ) && PETSC_VERSION_LESS_THAN( 3, 6, 0 )
#include <petsc-private/pcimpl.h>
#elif PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3, 3, 0 )
#include <petsc/private/pcimpl.h>
#else
#include <private/pcimpl.h>
#endif

typedef struct
{
    boost::shared_ptr<Feel::Preconditioner<double>> M_inHousePrec;
} PC_FEELPP;

#undef __FUNCT__
#define __FUNCT__ "PCSetUp_FEELPP"
static PetscErrorCode PCSetUp_FEELPP( PC pc )
{
    PC_FEELPP* pcfeelpp = (PC_FEELPP*)pc->data;
    PetscFunctionReturn( 0 );
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_FEELPP"
static PetscErrorCode PCApply_FEELPP( PC pc, Vec x, Vec y )
{
    PC_FEELPP* pcfeelpp = (PC_FEELPP*)pc->data;

    boost::shared_ptr<Feel::Preconditioner<double>> preconditioner = pcfeelpp->M_inHousePrec;
    // convert petsc Vec into feelpp Vector
    boost::shared_ptr<Feel::VectorPetsc<double>> x_vec;
    boost::shared_ptr<Feel::VectorPetsc<double>> y_vec;
    if ( preconditioner->worldComm().localSize() > 1 )
    {
        CHECK( preconditioner->matrix() ) << "matrix is not defined";
        x_vec.reset( new Feel::VectorPetscMPI<double>( x, preconditioner->matrix()->mapColPtr() ) );
        y_vec.reset( new Feel::VectorPetscMPI<double>( y, preconditioner->matrix()->mapRowPtr() ) );
    }
    else
    {
        x_vec.reset( new Feel::VectorPetsc<double>( x ) );
        y_vec.reset( new Feel::VectorPetsc<double>( y ) );
    }
    preconditioner->apply( *x_vec, *y_vec );

    PetscFunctionReturn( 0 );
}

#undef __FUNCT__
#define __FUNCT__ "PCReset_FEELPP"
static PetscErrorCode PCReset_FEELPP( PC pc )
{
    PetscFunctionReturn( 0 );
}

#undef __FUNCT__
#define __FUNCT__ "PCDestroy_FEELPP"
static PetscErrorCode PCDestroy_FEELPP( PC pc )
{
    PetscFunctionReturn( 0 );
}

#undef __FUNCT__
#define __FUNCT__ "PCSetFromOptions_FEELPP"
#if PETSC_VERSION_LESS_THAN( 3, 6, 0 )
static PetscErrorCode PCSetFromOptions_FEELPP( PC pc )
#else
static PetscErrorCode PCSetFromOptions_FEELPP( PetscOptions*, PC pc )
#endif
{
    PetscFunctionReturn( 0 );
}

#undef __FUNCT__
#define __FUNCT__ "PCView_FEELPP"
static PetscErrorCode PCView_FEELPP( PC pc, PetscViewer viewer )
{
    PetscFunctionReturn( 0 );
}

#undef __FUNCT__
#define __FUNCT__ "PCSetPrecond_FEELPP"
static PetscErrorCode PCSetPrecond_FEELPP( PC pc, boost::shared_ptr<Feel::Preconditioner<double>> const& pcInHouse )
{
    PC_FEELPP* pcfeelpp = (PC_FEELPP*)pc->data;
    pcfeelpp->M_inHousePrec = pcInHouse;
    PetscFunctionReturn( 0 );
}

#undef __FUNCT__
#define __FUNCT__ "PCCreate_FEELPP"
PETSC_EXTERN PetscErrorCode PCCreate_FEELPP( PC pc )
{
    PC_FEELPP* pcfeelpp;
    PetscErrorCode ierr;

    PetscFunctionBegin;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3, 5, 0 )
    ierr = PetscNewLog( pc, &pcfeelpp );
    CHKERRQ( ierr );
#else
    ierr = PetscNewLog( pc, PC_FEELPP, &pcfeelpp );
    CHKERRQ( ierr );
#endif
    pc->data = (void*)pcfeelpp;

    pc->ops->apply = PCApply_FEELPP;
    pc->ops->applytranspose = 0;
    pc->ops->setup = PCSetUp_FEELPP;
    pc->ops->reset = PCReset_FEELPP;
    pc->ops->destroy = PCDestroy_FEELPP;
    pc->ops->setfromoptions = PCSetFromOptions_FEELPP;
    pc->ops->view = PCView_FEELPP;
    pc->ops->applyrichardson = 0;
    PetscFunctionReturn( 0 );
}
