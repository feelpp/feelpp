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

typedef struct {
    PetscBool allocated;
    PetscBool scalediag;
    KSP       kspL;
    Vec       scale;
    Vec       x0,y0,x1;
    Mat L;             /* keep a copy to reuse when obtained with L = A10*A01 */
} PC_LSC;

namespace PetscImpl
{

static PetscErrorCode PCLSCGetKSP( PC pc, KSP& ksp )
{
    PC_LSC *mylsc = (PC_LSC*)pc->data;
    ksp = mylsc->kspL;
    PetscFunctionReturn(0);
}

static PetscErrorCode PCLSCSetScaleDiag( PC pc, PetscBool scalediag )
{
    PC_LSC *mylsc = (PC_LSC*)pc->data;
    mylsc->scalediag = scalediag;
    PetscFunctionReturn(0);
}

} // namespace PetscImpl

//------------------------------------------------------------------------------//

typedef struct {
  PetscInt cycles;                             /* Type of cycle to run: 1 V 2 W */
  PetscInt level;                              /* level = 0 coarsest level */
  PetscInt levels;                             /* number of active levels used */
  Vec      b;                                  /* Right hand side */
  Vec      x;                                  /* Solution */
  Vec      r;                                  /* Residual */

  PetscErrorCode (*residual)(Mat,Vec,Vec,Vec);

  Mat           A;                             /* matrix used in forming residual*/
  KSP           smoothd;                       /* pre smoother */
  KSP           smoothu;                       /* post smoother */
  Mat           interpolate;
  Mat           restrct;                       /* restrict is a reserved word in C99 and on Cray */
  Vec           rscale;                        /* scaling of restriction matrix */
  PetscLogEvent eventsmoothsetup;              /* if logging times for each level */
  PetscLogEvent eventsmoothsolve;
  PetscLogEvent eventresidual;
  PetscLogEvent eventinterprestrict;
} PC_MG_Levels;

typedef struct {
  PCMGType  am;                               /* Multiplicative, additive or full */
  PetscInt  cyclesperpcapply;                 /* Number of cycles to use in each PCApply(), multiplicative only*/
  PetscInt  maxlevels;                        /* total number of levels allocated */
  PetscInt  galerkin;                         /* use Galerkin process to compute coarser matrices, 0=no, 1=yes, 2=yes but computed externally */
  PetscBool usedmfornumberoflevels;           /* sets the number of levels by getting this information out of the DM */

  PetscInt     nlevels;
  PC_MG_Levels **levels;
  PetscInt     default_smoothu;               /* number of smooths per level if not over-ridden */
  PetscInt     default_smoothd;               /*  with calls to KSPSetTolerances() */
  PetscReal    rtol,abstol,dtol,ttol;         /* tolerances for when running with PCApplyRichardson_MG */

  void          *innerctx;                   /* optional data for preconditioner, like PCEXOTIC that inherits off of PCMG */
  PetscLogStage stageApply;
} PC_MG;


#if defined(PETSC_HAVE_ML)

EXTERN_C_BEGIN
/* HAVE_CONFIG_H flag is required by ML include files */
#if !defined(HAVE_CONFIG_H)
#define HAVE_CONFIG_H
#endif
#include <ml_include.h>
#include <ml_viz_stats.h>
EXTERN_C_END

typedef enum {PCML_NULLSPACE_AUTO,PCML_NULLSPACE_USER,PCML_NULLSPACE_BLOCK,PCML_NULLSPACE_SCALAR} PCMLNullSpaceType;
    //static const char *const PCMLNullSpaceTypes[] = {"AUTO","USER","BLOCK","SCALAR","PCMLNullSpaceType","PCML_NULLSPACE_",0};

/* The context (data structure) at each grid level */
typedef struct {
  Vec x,b,r;                  /* global vectors */
  Mat A,P,R;
  KSP ksp;
  Vec coords;                 /* projected by ML, if PCSetCoordinates is called; values packed by node */
} GridCtx;

/* The context used to input PETSc matrix into ML at fine grid */
typedef struct {
  Mat         A;       /* Petsc matrix in aij format */
  Mat         Aloc;    /* local portion of A to be used by ML */
  Vec         x,y;
  ML_Operator *mlmat;
  PetscScalar *pwork;  /* tmp array used by PetscML_comm() */
} FineGridCtx;
/* Private context for the ML preconditioner */
typedef struct {
  ML                *ml_object;
  ML_Aggregate      *agg_object;
  GridCtx           *gridctx;
  FineGridCtx       *PetscMLdata;
  PetscInt          Nlevels,MaxNlevels,MaxCoarseSize,CoarsenScheme,EnergyMinimization,MinPerProc,PutOnSingleProc,RepartitionType,ZoltanScheme;
  PetscReal         Threshold,DampingFactor,EnergyMinimizationDropTol,MaxMinRatio,AuxThreshold;
  PetscBool         SpectralNormScheme_Anorm,BlockScaling,EnergyMinimizationCheap,Symmetrize,OldHierarchy,KeepAggInfo,Reusable,Repartition,Aux;
  PetscBool         reuse_interpolation;
  PCMLNullSpaceType nulltype;
  PetscMPIInt       size; /* size of communicator for pc->pmat */
  PetscInt          dim;  /* data from PCSetCoordinates(_ML) */
  PetscInt          nloc;
  PetscReal         *coords; /* ML has a grid object for each level: the finest grid will point into coords */
} PC_ML;

namespace PetscImpl
{

static PetscErrorCode PCMLSetMaxNlevels( PC pc, PetscInt maxNLevel  )
{
    PC_MG *mg    = (PC_MG*)pc->data;
    PC_ML *pcml  = (PC_ML*)mg->innerctx;
    pcml->MaxNlevels = maxNLevel;
    PetscFunctionReturn(0);
}
static PetscErrorCode PCMLSetReuseInterpolation( PC pc, PetscBool reuse_interpolation )
{
    PC_MG *mg    = (PC_MG*)pc->data;
    PC_ML *pcml  = (PC_ML*)mg->innerctx;
    pcml->reuse_interpolation = reuse_interpolation;
    PetscFunctionReturn(0);
}
static PetscErrorCode PCMLSetKeepAggInfo( PC pc, PetscBool KeepAggInfo )
{
    PC_MG *mg    = (PC_MG*)pc->data;
    PC_ML *pcml  = (PC_ML*)mg->innerctx;
    pcml->KeepAggInfo = KeepAggInfo;
    PetscFunctionReturn(0);
}
static PetscErrorCode PCMLSetReusable( PC pc, PetscBool Reusable )
{
    PC_MG *mg    = (PC_MG*)pc->data;
    PC_ML *pcml  = (PC_ML*)mg->innerctx;
    pcml->Reusable = Reusable;
    PetscFunctionReturn(0);
}
static PetscErrorCode PCMLSetOldHierarchy( PC pc, PetscBool OldHierarchy )
{
    PC_MG *mg    = (PC_MG*)pc->data;
    PC_ML *pcml  = (PC_ML*)mg->innerctx;
    pcml->OldHierarchy = OldHierarchy;
    PetscFunctionReturn(0);
}

} // namespace PetscImpl


#endif // PETSC_HAVE_ML

#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )

#include <HYPRE_struct_mv.h>
#include <HYPRE_struct_ls.h>
#include <_hypre_struct_mv.h>
#include <HYPRE_sstruct_mv.h>
#include <HYPRE_sstruct_ls.h>
#include <_hypre_sstruct_mv.h>

/*
   Private context (data structure) for the  preconditioner.
*/
#if defined(PETSC_HAVE_HYPRE)
#if PETSC_VERSION_LESS_THAN( 3,6,0 )
typedef struct {
  HYPRE_Solver   hsolver;
  HYPRE_IJMatrix ij;
  HYPRE_IJVector b,x;

  HYPRE_Int (*destroy)(HYPRE_Solver);
  HYPRE_Int (*solve)(HYPRE_Solver,HYPRE_ParCSRMatrix,HYPRE_ParVector,HYPRE_ParVector);
  HYPRE_Int (*setup)(HYPRE_Solver,HYPRE_ParCSRMatrix,HYPRE_ParVector,HYPRE_ParVector);

  MPI_Comm comm_hypre;
  char     *hypre_type;

  /* options for Pilut and BoomerAMG*/
  PetscInt maxiter;
  double   tol;

  /* options for Pilut */
  PetscInt factorrowsize;

  /* options for ParaSails */
  PetscInt nlevels;
  double   threshhold;
  double   filter;
  PetscInt sym;
  double   loadbal;
  PetscInt logging;
  PetscInt ruse;
  PetscInt symt;

  /* options for Euclid */
  PetscBool bjilu;
  PetscInt  levels;

  /* options for Euclid and BoomerAMG */
  PetscBool printstatistics;

  /* options for BoomerAMG */
  PetscInt  cycletype;
  PetscInt  maxlevels;
  double    strongthreshold;
  double    maxrowsum;
  PetscInt  gridsweeps[3];
  PetscInt  coarsentype;
  PetscInt  measuretype;
  PetscInt  relaxtype[3];
  double    relaxweight;
  double    outerrelaxweight;
  PetscInt  relaxorder;
  double    truncfactor;
  PetscBool applyrichardson;
  PetscInt  pmax;
  PetscInt  interptype;
  PetscInt  agg_nl;
  PetscInt  agg_num_paths;
  PetscInt  nodal_coarsen;
  PetscBool nodal_relax;
  PetscInt  nodal_relax_levels;
} PC_HYPRE;
#else //PETSC_VERSION_LESS_THAN(3,6,0)
typedef struct {
  HYPRE_Solver   hsolver;
  HYPRE_IJMatrix ij;
  HYPRE_IJVector b,x;

  HYPRE_Int (*destroy)(HYPRE_Solver);
  HYPRE_Int (*solve)(HYPRE_Solver,HYPRE_ParCSRMatrix,HYPRE_ParVector,HYPRE_ParVector);
  HYPRE_Int (*setup)(HYPRE_Solver,HYPRE_ParCSRMatrix,HYPRE_ParVector,HYPRE_ParVector);
  HYPRE_Int (*setdgrad)(HYPRE_Solver,HYPRE_ParCSRMatrix);
  HYPRE_Int (*setdcurl)(HYPRE_Solver,HYPRE_ParCSRMatrix);
  HYPRE_Int (*setcoord)(HYPRE_Solver,HYPRE_ParVector,HYPRE_ParVector,HYPRE_ParVector);
  HYPRE_Int (*setdim)(HYPRE_Solver,HYPRE_Int);

  MPI_Comm comm_hypre;
  char     *hypre_type;

  /* options for Pilut and BoomerAMG*/
  PetscInt maxiter;
  double   tol;

  /* options for Pilut */
  PetscInt factorrowsize;

  /* options for ParaSails */
  PetscInt nlevels;
  double   threshhold;
  double   filter;
  PetscInt sym;
  double   loadbal;
  PetscInt logging;
  PetscInt ruse;
  PetscInt symt;

  /* options for BoomerAMG */
  PetscBool printstatistics;

  /* options for BoomerAMG */
  PetscInt  cycletype;
  PetscInt  maxlevels;
  double    strongthreshold;
  double    maxrowsum;
  PetscInt  gridsweeps[3];
  PetscInt  coarsentype;
  PetscInt  measuretype;
  PetscInt  relaxtype[3];
  double    relaxweight;
  double    outerrelaxweight;
  PetscInt  relaxorder;
  double    truncfactor;
  PetscBool applyrichardson;
  PetscInt  pmax;
  PetscInt  interptype;
  PetscInt  agg_nl;
  PetscInt  agg_num_paths;
  PetscInt  nodal_coarsen;
  PetscBool nodal_relax;
  PetscInt  nodal_relax_levels;

  /* options for AS (Auxiliary Space preconditioners) */
  PetscInt  as_print;
  PetscInt  as_max_iter;
  PetscReal as_tol;
  PetscInt  as_relax_type;
  PetscInt  as_relax_times;
  PetscReal as_relax_weight;
  PetscReal as_omega;
  PetscInt  as_amg_alpha_opts[5]; /* AMG coarsen type, agg_levels, relax_type, interp_type, Pmax for vector Poisson (AMS) or Curl problem (ADS) */
  PetscReal as_amg_alpha_theta;   /* AMG strength for vector Poisson (AMS) or Curl problem (ADS) */
  PetscInt  as_amg_beta_opts[5];  /* AMG coarsen type, agg_levels, relax_type, interp_type, Pmax for scalar Poisson (AMS) or vector Poisson (ADS) */
  PetscReal as_amg_beta_theta;    /* AMG strength for scalar Poisson (AMS) or vector Poisson (ADS)  */
  PetscInt  ams_cycle_type;
  PetscInt  ads_cycle_type;

  /* additional data */
  HYPRE_IJVector coords[3];
  HYPRE_IJVector constants[3];
  HYPRE_IJMatrix G;
  HYPRE_IJMatrix C;
  HYPRE_IJMatrix alpha_Poisson;
  HYPRE_IJMatrix beta_Poisson;
  PetscBool      ams_beta_is_zero;
} PC_HYPRE;
#endif //PETSC_VERSION_LESS_THAN( 3,6,0 )
#endif //defined(PETSC_HAVE_HYPRE)


namespace PetscImpl
{
#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_LESS_THAN( 3,6,0 )
static PetscErrorCode PCHYPRE_EUCLIDSetLevels( PC pc, PetscInt levels  )
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    PetscErrorCode ierr;

    jac->levels = levels;

    ierr = HYPRE_EuclidSetLevel( jac->hsolver,levels );
    //ierr = HYPRE_EuclidSetSparseA( jac->hsolver,2.22045e-14 );
    //ierr = HYPRE_EuclidSetILUT( jac->hsolver,2.22045e-14 );
    //ierr = HYPRE_EuclidSetStats( jac->hsolver,1);
    //ierr = HYPRE_EuclidSetMem( jac->hsolver,1);
    //ierr = HYPRE_EuclidDestroy( jac->hsolver );
    PetscFunctionReturn(0);
}
#endif

static PetscErrorCode PCHYPRE_BOOMERAMGSetMaxIter( PC pc, PetscInt maxIter  )
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->maxiter = maxIter;
    PetscErrorCode ierr;
    ierr = HYPRE_BoomerAMGSetMaxIter( jac->hsolver, jac->maxiter);
    PetscFunctionReturn(0);
}
static PetscErrorCode PCHYPRE_BOOMERAMGSetTol( PC pc, double tol )
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->tol=tol;
    PetscErrorCode ierr;
    ierr = HYPRE_BoomerAMGSetTol( jac->hsolver, jac->tol);
    PetscFunctionReturn(0);
}
static PetscErrorCode PCHYPRE_BOOMERAMGSetCycleType( PC pc, int cycleType ) 
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->cycletype = cycleType;
    PetscErrorCode ierr;
    ierr = HYPRE_BoomerAMGSetCycleType( jac->hsolver, jac->cycletype);
    PetscFunctionReturn(0);
}
static PetscErrorCode PCHYPRE_BOOMERAMGSetMaxLevels( PC pc, int maxLevels )
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->maxlevels=maxLevels;
    PetscErrorCode ierr;
    ierr = HYPRE_BoomerAMGSetMaxLevels( jac->hsolver, jac->maxlevels);
    PetscFunctionReturn(0);
}
static PetscErrorCode PCHYPRE_BOOMERAMGSetCoarsenType( PC pc, int coarsenType )
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->coarsentype=coarsenType;
    PetscErrorCode ierr;
    ierr = HYPRE_BoomerAMGSetCoarsenType( jac->hsolver, jac->coarsentype);
    PetscFunctionReturn(0);
}
static PetscErrorCode PCHYPRE_BOOMERAMGSetStrongThreshold( PC pc, double strongT )
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->strongthreshold = strongT;
    PetscErrorCode ierr;
    ierr = HYPRE_BoomerAMGSetStrongThreshold( jac->hsolver, jac->strongthreshold);
    PetscFunctionReturn(0);
}
static PetscErrorCode PCHYPRE_BOOMERAMGSetAggNumLevels( PC pc, int aggNl )
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->agg_nl = aggNl;
    PetscErrorCode ierr;
    ierr = HYPRE_BoomerAMGSetAggNumLevels( jac->hsolver, jac->agg_nl);
    PetscFunctionReturn(0);
}
static PetscErrorCode PCHYPRE_BOOMERAMGSetRelaxTypeAll( PC pc, int relaxTypeAll )
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->relaxtype[0] = jac->relaxtype[1]  = relaxTypeAll;
    PetscErrorCode ierr;
    ierr = HYPRE_BoomerAMGSetRelaxType( jac->hsolver, relaxTypeAll);
    jac->relaxtype[2] = 9;
    PetscFunctionReturn(0);
}
static PetscErrorCode PCHYPRE_BOOMERAMGSetInterpolationType( PC pc, int interpolationType )
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->interptype = interpolationType;
    PetscErrorCode ierr;
    ierr = HYPRE_BoomerAMGSetInterpType( jac->hsolver, interpolationType);
    PetscFunctionReturn(0);
}


#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,6,0 )
static PetscErrorCode PCHYPRE_AMSSetPrintLevel(PC pc, PetscInt printLevel)
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->as_print=printLevel;
    PetscErrorCode ierr;
    ierr = HYPRE_AMSSetPrintLevel(jac->hsolver,jac->as_print);
    CHKERRQ(ierr);
    PetscFunctionReturn(0.);
}
static PetscErrorCode PCHYPRE_AMSSetMaxIter(PC pc, PetscInt maxIter)
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->as_max_iter=maxIter;
    PetscErrorCode ierr;
    ierr = HYPRE_AMSSetMaxIter(jac->hsolver,jac->as_max_iter);
    CHKERRQ(ierr);
    PetscFunctionReturn(0.);
}
static PetscErrorCode PCHYPRE_AMSSetCycleType(PC pc, PetscInt cycleType)
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->ams_cycle_type=cycleType;
    PetscErrorCode ierr;
    ierr = HYPRE_AMSSetCycleType(jac->hsolver,jac->ams_cycle_type);
    CHKERRQ(ierr);
    PetscFunctionReturn(0.);
}
static PetscErrorCode PCHYPRE_AMSSetTol(PC pc, double tol)
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->as_tol=tol;
    PetscErrorCode ierr;
    ierr = HYPRE_AMSSetTol(jac->hsolver,jac->as_tol);
    CHKERRQ(ierr);
    PetscFunctionReturn(0.);
}
static PetscErrorCode PCHYPRE_AMSSetSmoothingOptions(PC pc, PetscInt relaxType, PetscInt relaxTimes, double relaxWeight, double omega)
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    PetscErrorCode ierr;
    jac->as_relax_type=relaxType;
    jac->as_relax_times=relaxTimes;
    jac->as_relax_weight=relaxWeight;
    jac->as_omega=omega;
    ierr = HYPRE_AMSSetSmoothingOptions(jac->hsolver,jac->as_relax_type,
            jac->as_relax_times,
            jac->as_relax_weight,
            jac->as_omega);
    CHKERRQ(ierr);
    PetscFunctionReturn(0.);
}
static PetscErrorCode PCHYPRE_AMSSetCoordinateVectors(PC pc,Vec x, Vec y, Vec z)
{
    PetscErrorCode ierr;
    PetscScalar *x_v;
    PetscScalar *y_v;
    PetscScalar *z_v;
    ierr = VecGetArray(x, &x_v); CHKERRQ(ierr);
    ierr = VecGetArray(y, &y_v); CHKERRQ(ierr);
    ierr = VecGetArray(z, &z_v); CHKERRQ(ierr);
    PetscReal *coord;
    PetscInt nloc;
    ierr = VecGetLocalSize(x, &nloc);
    CHKERRQ(ierr);
    coord = new PetscReal[3*nloc];
    for(int i = 0; i < 3*nloc; i++)
    {
      coord[i+0] = x_v[i];
      coord[i+1] = y_v[i];
      coord[i+2] = z_v[i];
    }
    ierr = PCSetCoordinates(pc,3, nloc, coord); CHKERRQ(ierr);
    ierr = VecRestoreArray(x, &x_v); CHKERRQ(ierr);
    ierr = VecRestoreArray(y, &y_v); CHKERRQ(ierr);
    ierr = VecRestoreArray(z, &z_v); CHKERRQ(ierr);
    delete [] coord;
    PetscFunctionReturn(0);
}
static PetscErrorCode PCHYPRE_AMSSetAlphaPoissonMatrix_HYPRE(PC pc, Mat G)
{
    PetscErrorCode ierr;
    ierr = PCHYPRESetAlphaPoissonMatrix(pc, G);
    CHKERRQ(ierr);
    PetscFunctionReturn(0.);
}
static PetscErrorCode PCHYPRE_AMSSetBetaPoissonMatrix_HYPRE(PC pc, Mat G)
{
    PetscErrorCode ierr;
    ierr =  PCHYPRESetBetaPoissonMatrix(pc, G);
    CHKERRQ(ierr);
    PetscFunctionReturn(0.);
}
#endif // PETSC_HAVE_HYPRE && PETSC >= 3.6
} // namespace PetscImpl
#endif // PETSC_HAVE_HYPRE


#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,4,0 )

typedef struct _PC_FieldSplitLink *PC_FieldSplitLink;
struct _PC_FieldSplitLink {
  KSP               ksp;
  Vec               x,y,z;
  char              *splitname;
  PetscInt          nfields;
  PetscInt          *fields,*fields_col;
  VecScatter        sctx;
  IS                is,is_col;
  PC_FieldSplitLink next,previous;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,7,0)
  PetscLogEvent     event;
#endif
};

typedef struct {
  PCCompositeType type;
  PetscBool       defaultsplit;                    /* Flag for a system with a set of 'k' scalar fields with the same layout (and bs = k) */
  PetscBool       splitdefined;                    /* Flag is set after the splits have been defined, to prevent more splits from being added */
  PetscInt        bs;                              /* Block size for IS and Mat structures */
  PetscInt        nsplits;                         /* Number of field divisions defined */
  Vec             *x,*y,w1,w2;
  Mat             *mat;                            /* The diagonal block for each split */
  Mat             *pmat;                           /* The preconditioning diagonal block for each split */
  Mat             *Afield;                         /* The rows of the matrix associated with each split */
  PetscBool       issetup;

  /* Only used when Schur complement preconditioning is used */
  Mat                       B;                     /* The (0,1) block */
  Mat                       C;                     /* The (1,0) block */
  Mat                       schur;                 /* The Schur complement S = A11 - A10 A00^{-1} A01, the KSP here, kspinner, is H_1 in [El08] */
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,5,0)
  Mat                       schurp;                /* Assembled approximation to S built by MatSchurComplement to be used as a preconditioning matrix when solving with S */
#endif
  Mat                       schur_user;            /* User-provided preconditioning matrix for the Schur complement */
  PCFieldSplitSchurPreType  schurpre;              /* Determines which preconditioning matrix is used for the Schur complement */
  PCFieldSplitSchurFactType schurfactorization;
  KSP                       kspschur;              /* The solver for S */
  KSP                       kspupper;              /* The solver for A in the upper diagonal part of the factorization (H_2 in [El08]) */
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,8,0)
  PetscScalar               schurscale;            /* Scaling factor for the Schur complement solution with DIAG factorization */
#endif

  PC_FieldSplitLink         head;
#if PETSC_VERSION_LESS_THAN(3,8,0)
  PetscBool                 reset;                  /* indicates PCReset() has been last called on this object, hack */
#endif
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,7,0)
  PetscBool                 isrestrict;             /* indicates PCFieldSplitRestrictIS() has been last called on this object, hack */
#endif
  PetscBool                 suboptionsset;          /* Indicates that the KSPSetFromOptions() has been called on the sub-KSPs */
  PetscBool                 dm_splits;              /* Whether to use DMCreateFieldDecomposition() whenever possible */
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,5,0)
  PetscBool                 diag_use_amat;          /* Whether to extract diagonal matrix blocks from Amat, rather than Pmat (weaker than -pc_use_amat) */
  PetscBool                 offdiag_use_amat;       /* Whether to extract off-diagonal matrix blocks from Amat, rather than Pmat (weaker than -pc_use_amat) */
#endif
} PC_FieldSplit;


namespace PetscImpl
{

/*
 * Warning : need to rewrite a petsc function for schur complement (in fieldsplit.c)
 * Allow to fix use of inner solver different of outer solver A^{-1}
 */
static PetscErrorCode  PCFieldSplitGetSubKSP_FieldSplit_Schur(PC pc,PetscInt *n,KSP **subksp)
{
  PC_FieldSplit  *jac = (PC_FieldSplit*)pc->data;
  PetscErrorCode ierr;

  //PetscFunctionBegin;
  ierr = PetscMalloc(jac->nsplits*sizeof(KSP),subksp);CHKERRQ(ierr);
  //ierr = MatSchurComplementGetKSP(jac->schur,*subksp);CHKERRQ(ierr);

  (*subksp)[0] = jac->head->ksp;
  (*subksp)[1] = jac->kspschur;
  if (n) *n = jac->nsplits;
  PetscFunctionReturn(0);
}

static PetscErrorCode PCFieldSplit_GetKSPInnerSchur( PC pc, KSP &ksp )
{
    PC_FieldSplit     *jac = (PC_FieldSplit*)pc->data;
    PetscErrorCode    ierr;

    CHECK( jac->schur ) << "mat schur complement is not initialized\n";

    KSP kspA = jac->head->ksp, kspInner = NULL;
    ierr  = MatSchurComplementGetKSP(jac->schur, &kspInner);CHKERRQ(ierr);
    if ( kspInner == NULL || kspInner==kspA )
    {
        ierr = KSPCreate(PetscObjectComm((PetscObject)pc), &ksp);CHKERRQ(ierr);
#if PETSC_VERSION_LESS_THAN(3,5,0)
        ierr = KSPSetOperators(ksp,jac->mat[0],jac->pmat[0],pc->flag);CHKERRQ(ierr);
#else
        ierr = KSPSetOperators(ksp,jac->mat[0],jac->pmat[0]/*,pc->flag*/);CHKERRQ(ierr);
#endif
        ierr = MatSchurComplementSetKSP(jac->schur,ksp);CHKERRQ(ierr);CHKERRQ(ierr);
    }
    else
    {
        ksp = kspInner;
        //ierr = KSPSetOperators(ksp,jac->mat[0],jac->pmat[0],pc->flag);CHKERRQ(ierr);
    }
    //KSPView( ksp, PETSC_VIEWER_STDOUT_WORLD );
    PetscFunctionReturn(0);
}

static PetscErrorCode PCFieldSplit_GetKSPUpperSchur( PC pc, KSP &ksp )
{
    PC_FieldSplit     *jac = (PC_FieldSplit*)pc->data;
    PetscErrorCode    ierr;

    KSP kspA = jac->head->ksp;
    if ( jac->kspupper == NULL || jac->kspupper == kspA )
    {
        ierr = KSPCreate(PetscObjectComm((PetscObject)pc), &ksp);CHKERRQ(ierr);
#if PETSC_VERSION_LESS_THAN(3,5,0)
        ierr = KSPSetOperators(ksp,jac->mat[0],jac->pmat[0],pc->flag);CHKERRQ(ierr);
#else
        ierr = KSPSetOperators(ksp,jac->mat[0],jac->pmat[0]/*,pc->flag*/);CHKERRQ(ierr);
#endif
        // need in PCApply_FieldSplit_Schur
        ierr = VecDuplicate(jac->head->x, &jac->head->z);CHKERRQ(ierr);

        jac->kspupper = ksp;
    }
    else
    {
        ksp = jac->kspupper;
    }

    PetscFunctionReturn(0);
}


static PetscErrorCode PCFieldSplit_GetMatSchurComplement( PC pc, Mat &schur )
{
    PC_FieldSplit     *jac = (PC_FieldSplit*)pc->data;
    PetscErrorCode    ierr;
    if ( !jac->schur )
        SETERRQ(((PetscObject) pc)->comm,PETSC_ERR_PLIB,"jac->schur is not initialized");

    schur = jac->schur;
    PetscFunctionReturn(0);
}

static PetscErrorCode PCFieldSplit_UpdateMatPrecondSchurComplement( PC pc, Mat schurMatPrecond )
{
    PC_FieldSplit     *jac = (PC_FieldSplit*)pc->data;
    PetscErrorCode    ierr;

    //jac->schur_user = schurMatPrecond;
#if PETSC_VERSION_LESS_THAN(3,5,0)
    ierr = KSPSetOperators(jac->kspschur,jac->schur,schurMatPrecond,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
#else
    ierr = KSPSetOperators(jac->kspschur,jac->schur,schurMatPrecond/*,DIFFERENT_NONZERO_PATTERN*/);CHKERRQ(ierr);
#endif

    PetscFunctionReturn(0);
}

} // namespace PetscImpl

#endif // PETSc Version

typedef struct {
  KSP          ksp;
  PC           pc;                   /* actual preconditioner used on each processor */
  Vec          xsub,ysub;            /* vectors of a subcommunicator to hold parallel vectors of PetscObjectComm((PetscObject)pc) */
  Vec          xdup,ydup;            /* parallel vector that congregates xsub or ysub facilitating vector scattering */
  Mat          pmats;                /* matrix and optional preconditioner matrix belong to a subcommunicator */
  VecScatter   scatterin,scatterout; /* scatter used to move all values to each processor group (subcommunicator) */
  PetscBool    useparallelmat;
  PetscSubcomm psubcomm;
  PetscInt     nsubcomm;           /* num of data structure PetscSubcomm */
} PC_Redundant;





#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )

namespace PetscImpl
{

#ifdef __FUNCT__
#undef __FUNCT__
#endif
#define __FUNCT__ "PCSetUp_Redundant"
static PetscErrorCode PCSetUp_Redundant(PC pc)
{
  PC_Redundant   *red = (PC_Redundant*)pc->data;
  PetscErrorCode ierr;
  PetscInt       mstart,mend,mlocal,M;
  PetscMPIInt    size;
  MPI_Comm       comm,subcomm;
  Vec            x;
  const char     *prefix;

  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)pc,&comm);CHKERRQ(ierr);
  
  /* if pmatrix set by user is sequential then we do not need to gather the parallel matrix */
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  if (size == 1) red->useparallelmat = PETSC_FALSE;

  if (!pc->setupcalled) {
    PetscInt mloc_sub;
    if (!red->psubcomm) {
      ierr = PetscSubcommCreate(comm,&red->psubcomm);CHKERRQ(ierr);
      ierr = PetscSubcommSetNumber(red->psubcomm,red->nsubcomm);CHKERRQ(ierr);
      ierr = PetscSubcommSetType(red->psubcomm,PETSC_SUBCOMM_CONTIGUOUS);CHKERRQ(ierr);
      /* enable runtime switch of psubcomm type, e.g., '-psubcomm_type interlaced */
      ierr = PetscSubcommSetFromOptions(red->psubcomm);CHKERRQ(ierr);
      ierr = PetscLogObjectMemory((PetscObject)pc,sizeof(PetscSubcomm));CHKERRQ(ierr);

      /* create a new PC that processors in each subcomm have copy of */
#if PETSC_VERSION_LESS_THAN(3,6,0)
      subcomm = red->psubcomm->comm;
#else
      subcomm = red->psubcomm->child;
#endif

      ierr = KSPCreate(subcomm,&red->ksp);CHKERRQ(ierr);
      ierr = PetscObjectIncrementTabLevel((PetscObject)red->ksp,(PetscObject)pc,1);CHKERRQ(ierr);
      ierr = PetscLogObjectParent((PetscObject)pc,(PetscObject)red->ksp);CHKERRQ(ierr);
      ierr = KSPSetType(red->ksp,KSPPREONLY);CHKERRQ(ierr);
      ierr = KSPGetPC(red->ksp,&red->pc);CHKERRQ(ierr);
      ierr = PCSetType(red->pc,PCLU);CHKERRQ(ierr);

      ierr = PCGetOptionsPrefix(pc,&prefix);CHKERRQ(ierr);
      ierr = KSPSetOptionsPrefix(red->ksp,prefix);CHKERRQ(ierr);
      ierr = KSPAppendOptionsPrefix(red->ksp,"redundant_");CHKERRQ(ierr);
    } else {
#if PETSC_VERSION_LESS_THAN(3,6,0)
      subcomm = red->psubcomm->comm;
#else
      subcomm = red->psubcomm->child;
#endif
    }

    if (red->useparallelmat) {
      /* grab the parallel matrix and put it into processors of a subcomminicator */
#if PETSC_VERSION_LESS_THAN(3,6,0)
      ierr = MatGetRedundantMatrix(pc->pmat,red->psubcomm->n,subcomm,MAT_INITIAL_MATRIX,&red->pmats);CHKERRQ(ierr);
#else
      ierr = MatCreateRedundantMatrix(pc->pmat,red->psubcomm->n,subcomm,MAT_INITIAL_MATRIX,&red->pmats);CHKERRQ(ierr);
#endif
      ierr = KSPSetOperators(red->ksp,red->pmats,red->pmats);CHKERRQ(ierr);
       
      /* get working vectors xsub and ysub */
#if PETSC_VERSION_LESS_THAN(3,6,0)
      ierr = MatGetVecs(red->pmats,&red->xsub,&red->ysub);CHKERRQ(ierr);
#else
      ierr = MatCreateVecs(red->pmats,&red->xsub,&red->ysub);CHKERRQ(ierr);
#endif

      /* create working vectors xdup and ydup.
       xdup concatenates all xsub's contigously to form a mpi vector over dupcomm  (see PetscSubcommCreate_interlaced())
       ydup concatenates all ysub and has empty local arrays because ysub's arrays will be place into it.
       Note: we use communicator dupcomm, not PetscObjectComm((PetscObject)pc)! */
      ierr = MatGetLocalSize(red->pmats,&mloc_sub,NULL);CHKERRQ(ierr);
      ierr = VecCreateMPI(red->psubcomm->dupparent,mloc_sub,PETSC_DECIDE,&red->xdup);CHKERRQ(ierr);
      ierr = VecCreateMPIWithArray(red->psubcomm->dupparent,1,mloc_sub,PETSC_DECIDE,NULL,&red->ydup);CHKERRQ(ierr);

      /* create vecscatters */
      if (!red->scatterin) { /* efficiency of scatterin is independent from psubcomm_type! */
        IS       is1,is2;
        PetscInt *idx1,*idx2,i,j,k;

#if PETSC_VERSION_LESS_THAN(3,6,0)
        ierr = MatGetVecs(pc->pmat,&x,0);CHKERRQ(ierr);
#else
        ierr = MatCreateVecs(pc->pmat,&x,0);CHKERRQ(ierr);
#endif
        ierr = VecGetSize(x,&M);CHKERRQ(ierr);
        ierr = VecGetOwnershipRange(x,&mstart,&mend);CHKERRQ(ierr);
        mlocal = mend - mstart;
        ierr = PetscMalloc2(red->psubcomm->n*mlocal,&idx1,red->psubcomm->n*mlocal,&idx2);CHKERRQ(ierr);
        j    = 0;
        for (k=0; k<red->psubcomm->n; k++) {
          for (i=mstart; i<mend; i++) {
            idx1[j]   = i;
            idx2[j++] = i + M*k;
          }
        }
        ierr = ISCreateGeneral(comm,red->psubcomm->n*mlocal,idx1,PETSC_COPY_VALUES,&is1);CHKERRQ(ierr);
        ierr = ISCreateGeneral(comm,red->psubcomm->n*mlocal,idx2,PETSC_COPY_VALUES,&is2);CHKERRQ(ierr);
        ierr = VecScatterCreate(x,is1,red->xdup,is2,&red->scatterin);CHKERRQ(ierr);
        ierr = ISDestroy(&is1);CHKERRQ(ierr);
        ierr = ISDestroy(&is2);CHKERRQ(ierr);

        /* Impl below is good for PETSC_SUBCOMM_INTERLACED (no inter-process communication) and PETSC_SUBCOMM_CONTIGUOUS (communication within subcomm) */
        ierr = ISCreateStride(comm,mlocal,mstart+ red->psubcomm->color*M,1,&is1);CHKERRQ(ierr);
        ierr = ISCreateStride(comm,mlocal,mstart,1,&is2);CHKERRQ(ierr);
        ierr = VecScatterCreate(red->xdup,is1,x,is2,&red->scatterout);CHKERRQ(ierr);
        ierr = ISDestroy(&is1);CHKERRQ(ierr);
        ierr = ISDestroy(&is2);CHKERRQ(ierr);
        ierr = PetscFree2(idx1,idx2);CHKERRQ(ierr);
        ierr = VecDestroy(&x);CHKERRQ(ierr);
      }
    } else { /* !red->useparallelmat */
      ierr = KSPSetOperators(red->ksp,pc->mat,pc->pmat);CHKERRQ(ierr);
    }
  } else { /* pc->setupcalled */
    if (red->useparallelmat) {
      MatReuse       reuse;
      /* grab the parallel matrix and put it into processors of a subcomminicator */
      /*--------------------------------------------------------------------------*/
      if (pc->flag == DIFFERENT_NONZERO_PATTERN) {
        /* destroy old matrices */
        ierr  = MatDestroy(&red->pmats);CHKERRQ(ierr);
        reuse = MAT_INITIAL_MATRIX;
      } else {
        reuse = MAT_REUSE_MATRIX;
      }
#if PETSC_VERSION_LESS_THAN(3,6,0)
      ierr = MatGetRedundantMatrix(pc->pmat,red->psubcomm->n,red->psubcomm->comm,reuse,&red->pmats);CHKERRQ(ierr);
#else
      ierr = MatCreateRedundantMatrix(pc->pmat,red->psubcomm->n,red->psubcomm->child,reuse,&red->pmats);CHKERRQ(ierr);
#endif
      ierr = KSPSetOperators(red->ksp,red->pmats,red->pmats);CHKERRQ(ierr);
    } else { /* !red->useparallelmat */
      ierr = KSPSetOperators(red->ksp,pc->mat,pc->pmat);CHKERRQ(ierr);
    }
  }

  if (pc->setfromoptionscalled) {
    ierr = KSPSetFromOptions(red->ksp);CHKERRQ(ierr);
  }
#if 0 // feelpp modif
  ierr = KSPSetUp(red->ksp);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}
#undef __FUNCT__
static PetscErrorCode PCRedundantChangeSetup(PC pc)
{
    pc->ops->setup          = PetscImpl::PCSetUp_Redundant;
    return(0);
}


} // namespace PetscImpl

#endif // PETSc version >= 3.5
