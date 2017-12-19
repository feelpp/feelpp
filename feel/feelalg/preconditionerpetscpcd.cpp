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
    KSP       kspA,kspQ;
    Mat       matA, matF, matAp, matQ, matMv, matBTBt;
    Vec       x1, x2, MvDiag;
    int       pcdOrder;
    char      pcdLaplacianType[256];
} PC_PCD_Feelpp;


#undef __FUNCT__
#define __FUNCT__ "PCSetUp_PCD_Feelpp"
static PetscErrorCode PCSetUp_PCD_Feelpp(PC pc)
{
    PetscErrorCode ierr;
    PC_PCD_Feelpp *pcpcd = (PC_PCD_Feelpp*)pc->data;

    if ( !pcpcd->allocated )
    {
        ierr = KSPCreate(PetscObjectComm((PetscObject)pc),&pcpcd->kspA);CHKERRQ(ierr);
        ierr = PetscObjectIncrementTabLevel((PetscObject)pcpcd->kspA,(PetscObject)pc,1);CHKERRQ(ierr);
        ierr = KSPSetType(pcpcd->kspA,KSPPREONLY);CHKERRQ(ierr);
        ierr = KSPSetOptionsPrefix(pcpcd->kspA,((PetscObject)pc)->prefix);CHKERRQ(ierr);
        ierr = KSPAppendOptionsPrefix(pcpcd->kspA,"pcd_A");CHKERRQ(ierr);
        ierr = KSPSetFromOptions(pcpcd->kspA);CHKERRQ(ierr);

        ierr = KSPCreate(PetscObjectComm((PetscObject)pc),&pcpcd->kspQ);CHKERRQ(ierr);
        ierr = PetscObjectIncrementTabLevel((PetscObject)pcpcd->kspQ,(PetscObject)pc,1);CHKERRQ(ierr);
        ierr = KSPSetType(pcpcd->kspQ,KSPPREONLY);CHKERRQ(ierr);
        ierr = KSPSetOptionsPrefix(pcpcd->kspQ,((PetscObject)pc)->prefix);CHKERRQ(ierr);
        ierr = KSPAppendOptionsPrefix(pcpcd->kspQ,"pcd_Q");CHKERRQ(ierr);
        ierr = KSPSetFromOptions(pcpcd->kspQ);CHKERRQ(ierr);

#if PETSC_VERSION_LESS_THAN(3,6,0)
        ierr = MatGetVecs(pcpcd->matF,&pcpcd->x2,&pcpcd->x1);CHKERRQ(ierr);
#else
        ierr = MatCreateVecs(pcpcd->matF,&pcpcd->x2,&pcpcd->x1);CHKERRQ(ierr);
#endif

        pcpcd->allocated = PETSC_TRUE;
    }

    PetscBool isBTBt;
    PetscStrcmp(pcpcd->pcdLaplacianType,"BTBt",&isBTBt);
    if ( isBTBt )
    {
        if (!pcpcd->MvDiag )
        {
#if PETSC_VERSION_LESS_THAN(3,6,0)
            ierr = MatGetVecs(pcpcd->matMv,&pcpcd->MvDiag,NULL);CHKERRQ(ierr);
#else
            ierr = MatCreateVecs(pcpcd->matMv,&pcpcd->MvDiag,NULL);CHKERRQ(ierr);
#endif
        }

        Mat B, C;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
        ierr = MatSchurComplementGetSubMatrices(pc->mat,NULL,NULL,&B,&C,NULL);CHKERRQ(ierr);
#else
        ierr = MatSchurComplementGetSubmatrices(pc->mat,NULL,NULL,&B,&C,NULL);CHKERRQ(ierr);
#endif
        ierr = MatGetDiagonal(pcpcd->matMv,pcpcd->MvDiag);CHKERRQ(ierr);
        ierr = VecReciprocal(pcpcd->MvDiag);CHKERRQ(ierr);
        // diag(F)^-1 * B
        ierr =  MatDiagonalScale( B, pcpcd->MvDiag ,NULL);CHKERRQ(ierr);
        // C* diag(F)^-1 * B
        if ( !pcpcd->matBTBt )
            ierr = MatMatMult(C,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&pcpcd->matBTBt);
        else
            ierr = MatMatMult(C,B,MAT_REUSE_MATRIX,PETSC_DEFAULT,&pcpcd->matBTBt);
        CHKERRQ(ierr);
        ierr = VecReciprocal(pcpcd->MvDiag);CHKERRQ(ierr);
        ierr =  MatDiagonalScale( B, pcpcd->MvDiag ,NULL);CHKERRQ(ierr);
        pcpcd->matAp = pcpcd->matBTBt;
    }
    else
        pcpcd->matAp = pcpcd->matA;

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
    ierr = KSPSetOperators(pcpcd->kspA,pcpcd->matAp,pcpcd->matAp);CHKERRQ(ierr);
    ierr = KSPSetOperators(pcpcd->kspQ,pcpcd->matQ,pcpcd->matQ);CHKERRQ(ierr);
#else
    ierr = KSPSetOperators(pcpcd->kspA,pcpcd->matAp,pcpcd->matAp,SAME_PRECONDITIONER);CHKERRQ(ierr);
    ierr = KSPSetOperators(pcpcd->kspQ,pcpcd->matQ,pcpcd->matQ,SAME_PRECONDITIONER);CHKERRQ(ierr);
#endif

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_PCD_Feelpp"
static PetscErrorCode PCApply_PCD_Feelpp(PC pc,Vec x,Vec y)
{
    PetscErrorCode ierr;
    PC_PCD_Feelpp         *pcpcd = (PC_PCD_Feelpp*)pc->data;

    if ( pcpcd->pcdOrder == 1 )
    {
        ierr = KSPSolve(pcpcd->kspQ,x,pcpcd->x1);CHKERRQ(ierr);
        ierr = MatMult(pcpcd->matF,pcpcd->x1,pcpcd->x2);CHKERRQ(ierr);
        ierr = KSPSolve(pcpcd->kspA,pcpcd->x2,y);CHKERRQ(ierr);
    }
    else
    {
        ierr = KSPSolve(pcpcd->kspA,x,pcpcd->x1);CHKERRQ(ierr);
        ierr = MatMult(pcpcd->matF,pcpcd->x1,pcpcd->x2);CHKERRQ(ierr);
        ierr = KSPSolve(pcpcd->kspQ,pcpcd->x2,y);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCReset_PCD_Feelpp"
static PetscErrorCode PCReset_PCD_Feelpp(PC pc)
{
    PC_PCD_Feelpp *pcpcd = (PC_PCD_Feelpp*)pc->data;
    KSPDestroy(&pcpcd->kspA);
    KSPDestroy(&pcpcd->kspQ);
    MatDestroy(&pcpcd->matA);
    MatDestroy(&pcpcd->matQ);
    MatDestroy(&pcpcd->matF);
    MatDestroy(&pcpcd->matBTBt);
    MatDestroy(&pcpcd->matMv);
    VecDestroy(&pcpcd->x1);
    VecDestroy(&pcpcd->x2);
    VecDestroy(&pcpcd->MvDiag);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCDestroy_PCD_Feelpp"
static PetscErrorCode PCDestroy_PCD_Feelpp(PC pc)
{
    PCReset_PCD_Feelpp(pc);
    PetscFree(pc->data);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCView_PCD_Feelpp"
static PetscErrorCode PCView_PCD_Feelpp(PC pc,PetscViewer viewer)
{
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCGetKSP_A_PCD_Feelpp"
static PetscErrorCode PCGetKSP_A_PCD_Feelpp( PC pc, KSP& ksp )
{
    PC_PCD_Feelpp *pcpcd = (PC_PCD_Feelpp*)pc->data;
    ksp = pcpcd->kspA;
    PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "PCGetKSP_Q_PCD_Feelpp"
static PetscErrorCode PCGetKSP_Q_PCD_Feelpp( PC pc, KSP& ksp )
{
    PC_PCD_Feelpp *pcpcd = (PC_PCD_Feelpp*)pc->data;
    ksp = pcpcd->kspQ;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetMatA_PCD_Feelpp"
static PetscErrorCode PCSetMatA_PCD_Feelpp(PC pc, Mat mat )
{
    PetscErrorCode ierr;
    PC_PCD_Feelpp *pcpcd = (PC_PCD_Feelpp*)pc->data;
    // increases the reference count for that object by one
    if ( mat ) { PetscObjectReference((PetscObject) mat); }
    MatDestroy(&pcpcd->matA);
    pcpcd->matA = mat;
    PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "PCSetMatQ_PCD_Feelpp"
static PetscErrorCode PCSetMatQ_PCD_Feelpp(PC pc, Mat mat )
{
    PetscErrorCode ierr;
    PC_PCD_Feelpp *pcpcd = (PC_PCD_Feelpp*)pc->data;
    // increases the reference count for that object by one
    if ( mat ) { PetscObjectReference((PetscObject) mat); }
    MatDestroy(&pcpcd->matQ);
    pcpcd->matQ = mat;
    PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "PCSetMatF_PCD_Feelpp"
static PetscErrorCode PCSetMatF_PCD_Feelpp(PC pc, Mat mat )
{
    PetscErrorCode ierr;
    PC_PCD_Feelpp *pcpcd = (PC_PCD_Feelpp*)pc->data;
    // increases the reference count for that object by one
    if ( mat ) { PetscObjectReference((PetscObject) mat); }
    MatDestroy(&pcpcd->matF);
    pcpcd->matF = mat;
    PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "PCSetMatMv_PCD_Feelpp"
static PetscErrorCode PCSetMatMv_PCD_Feelpp(PC pc, Mat mat )
{
    PetscErrorCode ierr;
    PC_PCD_Feelpp *pcpcd = (PC_PCD_Feelpp*)pc->data;
    // increases the reference count for that object by one
    if ( mat ) { PetscObjectReference((PetscObject) mat); }
    MatDestroy(&pcpcd->matMv);
    pcpcd->matMv = mat;
    PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "PCSetOrder_PCD_Feelpp"
static PetscErrorCode PCSetOrder_PCD_Feelpp(PC pc, int order )
{
    PetscErrorCode ierr;
    PC_PCD_Feelpp *pcpcd = (PC_PCD_Feelpp*)pc->data;
    pcpcd->pcdOrder = order;
    PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "PCSetOrder_PCD_Feelpp"
static PetscErrorCode PCSetLaplacianType_PCD_Feelpp(PC pc, const char * type )
{
    PetscErrorCode ierr;
    PC_PCD_Feelpp *pcpcd = (PC_PCD_Feelpp*)pc->data;
    PetscStrcpy(pcpcd->pcdLaplacianType, type);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCCreate_PCD_Feelpp"
PETSC_EXTERN PetscErrorCode PCCreate_PCD_Feelpp(PC pc)
{
  PC_PCD_Feelpp         *pcpcd;
  PetscErrorCode ierr;

  PetscFunctionBegin;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
  ierr     = PetscNewLog(pc,&pcpcd);CHKERRQ(ierr);
#else
  ierr     = PetscNewLog(pc,PC_PCD_Feelpp,&pcpcd);CHKERRQ(ierr);
#endif
  pc->data = (void*)pcpcd;

  pc->ops->apply           = PCApply_PCD_Feelpp;
  pc->ops->applytranspose  = 0;
  pc->ops->setup           = PCSetUp_PCD_Feelpp;
  pc->ops->reset           = PCReset_PCD_Feelpp;
  pc->ops->destroy         = PCDestroy_PCD_Feelpp;
  pc->ops->view            = PCView_PCD_Feelpp;
  pc->ops->applyrichardson = 0;

  pcpcd->pcdOrder = 1;
  PetscStrcpy(pcpcd->pcdLaplacianType, "Laplacian");
  pcpcd->MvDiag = NULL;
  pcpcd->matBTBt = NULL;
  PetscFunctionReturn(0);
}
