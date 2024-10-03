/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date:15 Dec 2017

 Copyright (C) 2017 Feel++ Consortium

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

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 ) && PETSC_VERSION_LESS_THAN( 3,6,0 )
#include <petsc-private/pcimpl.h>
#elif PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
#include <petsc/private/pcimpl.h>
#else
#include <private/pcimpl.h>
#endif

typedef struct {
    PetscBool allocated;
    KSP       kspPMM;
    Mat       matPMM;
} PC_PMM_Feelpp;


#undef __FUNCT__
#define __FUNCT__ "PCSetUp_PMM_Feelpp"
static PetscErrorCode PCSetUp_PMM_Feelpp(PC pc)
{
    PetscErrorCode ierr;
    PC_PMM_Feelpp *pcpmm = (PC_PMM_Feelpp*)pc->data;

    if ( !pcpmm->allocated )
    {
        ierr = KSPCreate(PetscObjectComm((PetscObject)pc),&pcpmm->kspPMM);CHKERRQ(ierr);
        ierr = PetscObjectIncrementTabLevel((PetscObject)pcpmm->kspPMM,(PetscObject)pc,1);CHKERRQ(ierr);
        ierr = KSPSetType(pcpmm->kspPMM,KSPPREONLY);CHKERRQ(ierr);
        ierr = KSPSetOptionsPrefix(pcpmm->kspPMM,((PetscObject)pc)->prefix);CHKERRQ(ierr);
        ierr = KSPAppendOptionsPrefix(pcpmm->kspPMM,"pmm_");CHKERRQ(ierr);
        ierr = KSPSetFromOptions(pcpmm->kspPMM);CHKERRQ(ierr);
        pcpmm->allocated = PETSC_TRUE;
    }

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
    ierr = KSPSetOperators(pcpmm->kspPMM,pcpmm->matPMM,pcpmm->matPMM);CHKERRQ(ierr);
#else
    ierr = KSPSetOperators(pcpmm->kspPMM,pcpmm->matPMM,pcpmm->matPMM,SAME_PRECONDITIONER);CHKERRQ(ierr);
#endif

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_PMM_Feelpp"
static PetscErrorCode PCApply_PMM_Feelpp(PC pc,Vec x,Vec y)
{
    PetscErrorCode ierr;
    PC_PMM_Feelpp         *pcpmm = (PC_PMM_Feelpp*)pc->data;
    ierr = KSPSolve(pcpmm->kspPMM,x,y);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCReset_PMM_Feelpp"
static PetscErrorCode PCReset_PMM_Feelpp(PC pc)
{
    PC_PMM_Feelpp *pcpmm = (PC_PMM_Feelpp*)pc->data;
    KSPDestroy(&pcpmm->kspPMM);
    MatDestroy(&pcpmm->matPMM);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCDestroy_PMM_Feelpp"
static PetscErrorCode PCDestroy_PMM_Feelpp(PC pc)
{
    PCReset_PMM_Feelpp(pc);
    PetscFree(pc->data);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCView_PMM_Feelpp"
static PetscErrorCode PCView_PMM_Feelpp(PC pc,PetscViewer viewer)
{
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCGetKSP_PMM_Feelpp"
static PetscErrorCode PCGetKSP_PMM_Feelpp( PC pc, KSP& ksp )
{
    PC_PMM_Feelpp *pcpmm = (PC_PMM_Feelpp*)pc->data;
    ksp = pcpmm->kspPMM;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetPressureMassMatrix_PMM_Feelpp"
static PetscErrorCode PCSetPressureMassMatrix_PMM_Feelpp(PC pc, Mat mat )
{
    PetscErrorCode ierr;
    PC_PMM_Feelpp *pcpmm = (PC_PMM_Feelpp*)pc->data;

    //PetscFunctionBegin;
    // increases the reference count for that object by one
    if ( mat ) { PetscObjectReference((PetscObject) mat); }
    MatDestroy(&pcpmm->matPMM);
    pcpmm->matPMM = mat;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCCreate_PMM_Feelpp"
PETSC_EXTERN PetscErrorCode PCCreate_PMM_Feelpp(PC pc)
{
  PC_PMM_Feelpp         *pcpmm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,18,0 )
  ierr     = PetscNew(&pcpmm);CHKERRQ(ierr);  
#elif PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
  ierr     = PetscNewLog(pc,&pcpmm);CHKERRQ(ierr);
#else
  ierr     = PetscNewLog(pc,PC_PMM_Feelpp,&pcpmm);CHKERRQ(ierr);
#endif
  pc->data = (void*)pcpmm;

  pc->ops->apply           = PCApply_PMM_Feelpp;
  pc->ops->applytranspose  = 0;
  pc->ops->setup           = PCSetUp_PMM_Feelpp;
  pc->ops->reset           = PCReset_PMM_Feelpp;
  pc->ops->destroy         = PCDestroy_PMM_Feelpp;
  pc->ops->view            = PCView_PMM_Feelpp;
  pc->ops->applyrichardson = 0;
  PetscFunctionReturn(0);
}
