
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 ) && PETSC_VERSION_LESS_THAN( 3,6,0 )
#include <petsc-private/pcimpl.h>
#elif PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
#include <petsc/private/pcimpl.h>
#else
#include <private/pcimpl.h>
#endif

typedef struct {
  PetscBool allocated;
  PetscBool scalediag;
  KSP       kspL;
  Vec       scale;
  Vec       x0,y0,x1;
  Mat       L;             /* keep a copy to reuse when obtained with L = A10*A01 */
  Mat       massMatrix;
} PC_LSC2;

#undef __FUNCT__
#define __FUNCT__ "PCLSC2Allocate_Private"
static PetscErrorCode PCLSC2Allocate_Private(PC pc)
{
  PC_LSC2         *lsc = (PC_LSC2*)pc->data;
  Mat            A;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (lsc->allocated) PetscFunctionReturn(0);
  ierr = KSPCreate(PetscObjectComm((PetscObject)pc),&lsc->kspL);CHKERRQ(ierr);
  ierr = PetscObjectIncrementTabLevel((PetscObject)lsc->kspL,(PetscObject)pc,1);CHKERRQ(ierr);
  ierr = KSPSetType(lsc->kspL,KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(lsc->kspL,((PetscObject)pc)->prefix);CHKERRQ(ierr);
  ierr = KSPAppendOptionsPrefix(lsc->kspL,"lsc_");CHKERRQ(ierr);
  ierr = KSPSetFromOptions(lsc->kspL);CHKERRQ(ierr);
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
  ierr = MatSchurComplementGetSubMatrices(pc->mat,&A,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
#else
  ierr = MatSchurComplementGetSubmatrices(pc->mat,&A,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
#endif
#if PETSC_VERSION_LESS_THAN(3,6,0)
  ierr = MatGetVecs(A,&lsc->x0,&lsc->y0);CHKERRQ(ierr);
  ierr = MatGetVecs(pc->pmat,&lsc->x1,NULL);CHKERRQ(ierr);
#else
  ierr = MatCreateVecs(A,&lsc->x0,&lsc->y0);CHKERRQ(ierr);
  ierr = MatCreateVecs(pc->pmat,&lsc->x1,NULL);CHKERRQ(ierr);
#endif
  if (lsc->scalediag) {
    ierr = VecDuplicate(lsc->x0,&lsc->scale);CHKERRQ(ierr);
  }
  lsc->allocated = PETSC_TRUE;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetUp_LSC2"
static PetscErrorCode PCSetUp_LSC2(PC pc)
{
  PC_LSC2         *lsc = (PC_LSC2*)pc->data;
  Mat            L,Lp,B,C;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PCLSC2Allocate_Private(pc);CHKERRQ(ierr);


  if (lsc->scale) {
      if ( lsc->massMatrix )
      {
          ierr = MatGetDiagonal(lsc->massMatrix,lsc->scale);CHKERRQ(ierr);
      }
      else
      {
          Mat Ap;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
          ierr = MatSchurComplementGetSubMatrices(pc->mat,NULL,&Ap,NULL,NULL,NULL);CHKERRQ(ierr);
#else
          ierr = MatSchurComplementGetSubmatrices(pc->mat,NULL,&Ap,NULL,NULL,NULL);CHKERRQ(ierr);
#endif
          ierr = MatGetDiagonal(Ap,lsc->scale);CHKERRQ(ierr); /* Should be the mass matrix, but we don't have plumbing for that yet */
      }
      ierr = VecReciprocal(lsc->scale);CHKERRQ(ierr);
  }


  ierr = PetscObjectQuery((PetscObject)pc->mat,"LSC_L",(PetscObject*)&L);CHKERRQ(ierr);
  if (!L) {ierr = PetscObjectQuery((PetscObject)pc->pmat,"LSC_L",(PetscObject*)&L);CHKERRQ(ierr);}
  ierr = PetscObjectQuery((PetscObject)pc->pmat,"LSC_Lp",(PetscObject*)&Lp);CHKERRQ(ierr);
  if (!Lp) {ierr = PetscObjectQuery((PetscObject)pc->mat,"LSC_Lp",(PetscObject*)&Lp);CHKERRQ(ierr);}
  if (!L) {
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
    ierr = MatSchurComplementGetSubMatrices(pc->mat,NULL,NULL,&B,&C,NULL);CHKERRQ(ierr);
#else
    ierr = MatSchurComplementGetSubmatrices(pc->mat,NULL,NULL,&B,&C,NULL);CHKERRQ(ierr);
#endif
    if (!lsc->L) {
        if (lsc->scale)
        {
            // diag(F)^-1 * B
            ierr =  MatDiagonalScale( B, lsc->scale ,NULL);CHKERRQ(ierr);
        }
        // C* diag(F)^-1 * B
        ierr = MatMatMult(C,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&lsc->L);CHKERRQ(ierr);
        // revert B and lsc->scale
        if (lsc->scale)
        {
            ierr = VecReciprocal(lsc->scale);CHKERRQ(ierr);
            ierr = MatDiagonalScale( B, lsc->scale ,NULL);CHKERRQ(ierr);
            ierr = VecReciprocal(lsc->scale);CHKERRQ(ierr);
        }
    } else {
        if (lsc->scale)
        {
            // diag(F)^-1 * B
            ierr =  MatDiagonalScale( B, lsc->scale ,NULL);CHKERRQ(ierr);
        }
        // C* diag(F)^-1 * B
        ierr = MatMatMult(C,B,MAT_REUSE_MATRIX,PETSC_DEFAULT,&lsc->L);CHKERRQ(ierr);
        // revert B and lsc->scale
        if (lsc->scale)
        {
            ierr = VecReciprocal(lsc->scale);CHKERRQ(ierr);
            ierr = MatDiagonalScale( B, lsc->scale ,NULL);CHKERRQ(ierr);
            ierr = VecReciprocal(lsc->scale);CHKERRQ(ierr);
        }
    }
    Lp = L = lsc->L;
  }
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
  ierr = KSPSetOperators(lsc->kspL,L,Lp);CHKERRQ(ierr);
#else
  ierr = KSPSetOperators(lsc->kspL,L,Lp,SAME_PRECONDITIONER);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_LSC2"
static PetscErrorCode PCApply_LSC2(PC pc,Vec x,Vec y)
{
  PC_LSC2         *lsc = (PC_LSC2*)pc->data;
  Mat            A,B,C;
  PetscErrorCode ierr;

  PetscFunctionBegin;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
  ierr = MatSchurComplementGetSubMatrices(pc->mat,&A,NULL,&B,&C,NULL);CHKERRQ(ierr);
#else
  ierr = MatSchurComplementGetSubmatrices(pc->mat,&A,NULL,&B,&C,NULL);CHKERRQ(ierr);
#endif
  ierr = KSPSolve(lsc->kspL,x,lsc->x1);CHKERRQ(ierr);
  ierr = MatMult(B,lsc->x1,lsc->x0);CHKERRQ(ierr);
  if (lsc->scale) {
    ierr = VecPointwiseMult(lsc->x0,lsc->x0,lsc->scale);CHKERRQ(ierr);
  }
  ierr = MatMult(A,lsc->x0,lsc->y0);CHKERRQ(ierr);
  if (lsc->scale) {
    ierr = VecPointwiseMult(lsc->y0,lsc->y0,lsc->scale);CHKERRQ(ierr);
  }
  ierr = MatMult(C,lsc->y0,lsc->x1);CHKERRQ(ierr);
  ierr = KSPSolve(lsc->kspL,lsc->x1,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCReset_LSC2"
static PetscErrorCode PCReset_LSC2(PC pc)
{
  PC_LSC2         *lsc = (PC_LSC2*)pc->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecDestroy(&lsc->x0);CHKERRQ(ierr);
  ierr = VecDestroy(&lsc->y0);CHKERRQ(ierr);
  ierr = VecDestroy(&lsc->x1);CHKERRQ(ierr);
  ierr = VecDestroy(&lsc->scale);CHKERRQ(ierr);
  ierr = KSPDestroy(&lsc->kspL);CHKERRQ(ierr);
  ierr = MatDestroy(&lsc->L);CHKERRQ(ierr);
  ierr = MatDestroy(&lsc->massMatrix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCDestroy_LSC2"
static PetscErrorCode PCDestroy_LSC2(PC pc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PCReset_LSC2(pc);CHKERRQ(ierr);
  ierr = PetscFree(pc->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetFromOptions_LSC2"
#if PETSC_VERSION_LESS_THAN(3,6,0)
static PetscErrorCode PCSetFromOptions_LSC2(PC pc)
#else
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,7,0)
static PetscErrorCode PCSetFromOptions_LSC2(PetscOptionItems*, PC pc)
#else
static PetscErrorCode PCSetFromOptions_LSC2(PetscOptions*, PC pc)
#endif
#endif
{
#if 0
    PC_LSC2         *lsc = (PC_LSC2*)pc->data;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = PetscOptionsHead("LSC2 options");CHKERRQ(ierr);
    {
        ierr = PetscOptionsBool("-pc_lsc_scale_diag","Use diagonal of velocity block (A) for scaling","None",lsc->scalediag,&lsc->scalediag,NULL);CHKERRQ(ierr);
    }
    ierr = PetscOptionsTail();CHKERRQ(ierr);
#endif
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCView_LSC2"
static PetscErrorCode PCView_LSC2(PC pc,PetscViewer viewer)
{
  PC_LSC2         *jac = (PC_LSC2*)pc->data;
  PetscErrorCode ierr;
  PetscBool      iascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = KSPView(jac->kspL,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    if (jac->scalediag)
        PetscViewerASCIIPrintf(viewer, "use scalediag\n" );
    else
        PetscViewerASCIIPrintf(viewer, "not use scalediag\n" );

    if ( jac->massMatrix )
        PetscViewerASCIIPrintf(viewer, "has a mass matrix register\n" );

  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetMassMatrix_LSC2"
static PetscErrorCode PCSetMassMatrix_LSC2(PC pc, Mat mat )
{
  PC_LSC2         *lsc = (PC_LSC2*)pc->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  // increases the reference count for that object by one
  if ( mat ) { PetscObjectReference((PetscObject) mat); }
  MatDestroy(&lsc->massMatrix);
  lsc->massMatrix = mat;
  PetscFunctionReturn(0);
}

static PetscErrorCode PCLSC2GetKSP( PC pc, KSP& ksp )
{
    PC_LSC2 *mylsc = (PC_LSC2*)pc->data;
    ksp = mylsc->kspL;
    PetscFunctionReturn(0);
 }

 static PetscErrorCode PCLSC2SetScaleDiag( PC pc, PetscBool scalediag )
 {
     PC_LSC2 *mylsc = (PC_LSC2*)pc->data;
     mylsc->scalediag = scalediag;
     PetscFunctionReturn(0);
 }


/*MC
     PCLSC2 - Preconditioning for Schur complements, based on Least Squares Commutators

   Options Database Key:
.    -pc_lsc_scale_diag - Use the diagonal of A for scaling

   Level: intermediate

   Notes:
   This preconditioner will normally be used with PCFieldSplit to precondition the Schur complement, but
   it can be used for any Schur complement system.  Consider the Schur complement

.vb
   S = A11 - A10 inv(A00) A01
.ve

   PCLSC2 currently doesn't do anything with A11, so let's assume it is 0.  The idea is that a good approximation to
   inv(S) is given by

.vb
   inv(A10 A01) A10 A00 A01 inv(A10 A01)
.ve

   The product A10 A01 can be computed for you, but you can provide it (this is
   usually more efficient anyway).  In the case of incompressible flow, A10 A10 is a Laplacian, call it L.  The current
   interface is to hang L and a preconditioning matrix Lp on the preconditioning matrix.

   If you had called KSPSetOperators(ksp,S,Sp), S should have type MATSCHURCOMPLEMENT and Sp can be any type you
   like (PCLSC2 doesn't use it directly) but should have matrices composed with it, under the names "LSC2_L" and "LSC2_Lp".
   For example, you might have setup code like this

.vb
   PetscObjectCompose((PetscObject)Sp,"LSC2_L",(PetscObject)L);
   PetscObjectCompose((PetscObject)Sp,"LSC2_Lp",(PetscObject)Lp);
.ve

   And then your Jacobian assembly would look like

.vb
   PetscObjectQuery((PetscObject)Sp,"LSC2_L",(PetscObject*)&L);
   PetscObjectQuery((PetscObject)Sp,"LSC2_Lp",(PetscObject*)&Lp);
   if (L) { assembly L }
   if (Lp) { assemble Lp }
.ve

   With this, you should be able to choose LSC2 preconditioning, using e.g. ML's algebraic multigrid to solve with L

.vb
   -fieldsplit_1_pc_type lsc -fieldsplit_1_lsc_pc_type ml
.ve

   Since we do not use the values in Sp, you can still put an assembled matrix there to use normal preconditioners.

   References:
+  Elman, Howle, Shadid, Shuttleworth, and Tuminaro, Block preconditioners based on approximate commutators, 2006.
-  Silvester, Elman, Kay, Wathen, Efficient preconditioning of the linearized Navier-Stokes equations for incompressible flow, 2001.

   Concepts: physics based preconditioners, block preconditioners

.seealso:  PCCreate(), PCSetType(), PCType (for list of available types), PC, Block_Preconditioners, PCFIELDSPLIT,
           PCFieldSplitGetSubKSP(), PCFieldSplitSetFields(), PCFieldSplitSetType(), PCFieldSplitSetIS(), PCFieldSplitSetSchurPre(),
           MatCreateSchurComplement()
M*/

#undef __FUNCT__
#define __FUNCT__ "PCCreate_LSC2"
PETSC_EXTERN PetscErrorCode PCCreate_LSC2(PC pc)
{
  PC_LSC2         *lsc;
  PetscErrorCode ierr;

  PetscFunctionBegin;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,18,0 )
  ierr     = PetscNew(&lsc);CHKERRQ(ierr);  
#elif PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
  ierr     = PetscNewLog(pc,&lsc);CHKERRQ(ierr);
#else
  ierr     = PetscNewLog(pc,PC_LSC2,&lsc);CHKERRQ(ierr);
#endif
  pc->data = (void*)lsc;

  pc->ops->apply           = PCApply_LSC2;
  pc->ops->applytranspose  = 0;
  pc->ops->setup           = PCSetUp_LSC2;
  pc->ops->reset           = PCReset_LSC2;
  pc->ops->destroy         = PCDestroy_LSC2;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,7,0 )
#else
  pc->ops->setfromoptions  = PCSetFromOptions_LSC2;
#endif
  pc->ops->view            = PCView_LSC2;
  pc->ops->applyrichardson = 0;
  PetscFunctionReturn(0);
}
