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
   \file preconditionerpetsc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-01-16
 */
#include <feel/feelalg/preconditionerpetsc.hpp>
#include <feel/feelalg/functionspetsc.hpp>
#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>


extern "C" {

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
#include <petsc-private/pcimpl.h>
#include <petsc-private/kspimpl.h>
#else
#include <private/pcimpl.h>
#include <private/kspimpl.h>
#endif

//------------------------------------------------------------------------------//
typedef struct {
    PetscBool allocated;
    PetscBool scalediag;
    KSP       kspL;
    Vec       scale;
    Vec       x0,y0,x1;
    Mat L;             /* keep a copy to reuse when obtained with L = A10*A01 */
} PC_LSC;

static void PCLSCGetKSP( PC pc, KSP& ksp )
{
    PC_LSC *mylsc = (PC_LSC*)pc->data;
    ksp = mylsc->kspL;
}
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

static void PCMLSetMaxNlevels( PC pc, PetscInt maxNLevel  )
{
    PC_MG *mg    = (PC_MG*)pc->data;
    PC_ML *pcml  = (PC_ML*)mg->innerctx;
    pcml->MaxNlevels = maxNLevel;
}
static void PCMLSetReuseInterpolation( PC pc, PetscBool reuse_interpolation )
{
    PC_MG *mg    = (PC_MG*)pc->data;
    PC_ML *pcml  = (PC_ML*)mg->innerctx;
    pcml->reuse_interpolation = reuse_interpolation;
}
static void PCMLSetKeepAggInfo( PC pc, PetscBool KeepAggInfo )
{
    PC_MG *mg    = (PC_MG*)pc->data;
    PC_ML *pcml  = (PC_ML*)mg->innerctx;
    pcml->KeepAggInfo = KeepAggInfo;
}
static void PCMLSetReusable( PC pc, PetscBool Reusable )
{
    PC_MG *mg    = (PC_MG*)pc->data;
    PC_ML *pcml  = (PC_ML*)mg->innerctx;
    pcml->Reusable = Reusable;
}
static void PCMLSetOldHierarchy( PC pc, PetscBool OldHierarchy )
{
    PC_MG *mg    = (PC_MG*)pc->data;
    PC_ML *pcml  = (PC_ML*)mg->innerctx;
    pcml->OldHierarchy = OldHierarchy;
}

#endif // PETSC_HAVE_ML

#if defined(PETSC_HAVE_HYPRE)

#include <HYPRE_struct_mv.h>
#include <HYPRE_struct_ls.h>
#include <_hypre_struct_mv.h>
#include <HYPRE_sstruct_mv.h>
#include <HYPRE_sstruct_ls.h>
#include <_hypre_sstruct_mv.h>

/*
   Private context (data structure) for the  preconditioner.
*/
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

static void PCHYPRE_EUCLIDSetLevels( PC pc, PetscInt levels  )
{
    PC_HYPRE       *jac = (PC_HYPRE*)pc->data;
    jac->levels = levels;
}

#endif // PETSC_HAVE_HYPRE



} // extern C


namespace Feel
{
template <typename T>
void PreconditionerPetsc<T>::apply( const Vector<T> & x, Vector<T> & y )
{
    if ( !this->M_is_initialized ) this->init();


    VectorPetsc<T> & x_pvec = dynamic_cast<VectorPetsc<T>&>( const_cast<Vector<T>&>( x ) );
    VectorPetsc<T> & y_pvec = dynamic_cast<VectorPetsc<T>&>( const_cast<Vector<T>&>( y ) );

    Vec x_vec = x_pvec.vec();
    Vec y_vec = y_pvec.vec();

    int ierr = PCApply( M_pc,x_vec,y_vec );
    CHKERRABORT( this->worldComm().globalComm(),ierr );
}


/*----------------------- inline functions ----------------------------------*/
template <typename T>
PreconditionerPetsc<T>::PreconditionerPetsc ( std::string const& name, WorldComm const& worldComm )
    :
    Preconditioner<T>( name, worldComm )
{
}
template <typename T>
PreconditionerPetsc<T>::PreconditionerPetsc ( PreconditionerPetsc const& p )
    :
    Preconditioner<T>( p )
{
}



template <typename T>
PreconditionerPetsc<T>::~PreconditionerPetsc ()
{
    this->clear ();
}


template <typename T>
void PreconditionerPetsc<T>::init ()
{
    CHECK( this->M_matrix ) << "ERROR: No matrix set for PreconditionerPetsc, but init() called" << "\n";
    this->M_matrix->close();
    indexsplit_ptrtype is;
    // Clear the preconditioner in case it has been created in the past
    if ( !this->M_is_initialized )
    {
        // Create the preconditioning object
        check( PCCreate( this->worldComm().globalComm(),&M_pc ) );
        check( PCSetFromOptions ( M_pc ) );
#if PETSC_VERSION_LESS_THAN(3,4,0)
        const PCType pc_type;
#else
        PCType pc_type;
#endif
        check( PCGetType ( M_pc, &pc_type ) );

        MatrixPetsc<T> * pmatrix = dynamic_cast<MatrixPetsc<T>*>( this->M_matrix.get() );

        M_mat = pmatrix->mat();

        if (this->M_preconditioner_type==FIELDSPLIT_PRECOND )
        {
            check( PCSetType( M_pc,( char* ) PCFIELDSPLIT ) );

            is = pmatrix->indexSplit();
            //is.showMe();

            std::string fieldsDefStr = option( _prefix=this->name(), _name="fieldsplit-fields" ).template as<std::string>();
            auto fieldsDef = IndexSplit::parseFieldsDef( fieldsDefStr );
            if ( fieldsDef.size() == 0 )
            {
                pmatrix->updatePCFieldSplit( M_pc );
            }
            else
            {
                //fieldsDef.showMe();
                auto isUsed = pmatrix->indexSplit()->applyFieldsDef( fieldsDef  );
                //isUsed.showMe();
                pmatrix->updatePCFieldSplit( M_pc,isUsed );
            }
        }

    }
    else if (this->M_mat_has_changed)
    {
        MatrixPetsc<T> * pmatrix = dynamic_cast<MatrixPetsc<T>*>( this->M_matrix.get() );
        M_mat = pmatrix->mat();
        if (this->M_preconditioner_type==FIELDSPLIT_PRECOND )
        {
            check( PCSetType( M_pc,( char* ) PCFIELDSPLIT ) );
            pmatrix->updatePCFieldSplit( M_pc );
            is=pmatrix->indexSplit();
        }
        this->M_mat_has_changed = false;
    }

    //check( PCSetOperators( M_pc,M_mat,M_mat, PetscGetMatStructureEnum(MatrixStructure::SAME_NONZERO_PATTERN) ) );
    //check( PCSetOperators( M_pc,M_mat,M_mat, PetscGetMatStructureEnum(MatrixStructure::DIFFERENT_NONZERO_PATTERN) ) );
    check( PCSetOperators( M_pc,M_mat,M_mat, PetscGetMatStructureEnum(this->M_prec_matrix_structure) ) );

    // Set the PCType.  Note: this used to be done *before* the call to
    // PCSetOperators(), and only when !M_is_initialized, but
    // 1.) Some preconditioners (those employing sub-preconditioners,
    // for example) have to call PCSetUp(), and can only do this after
    // the operators have been set.
    // 2.) It should be safe to call set_petsc_preconditioner_type()
    // multiple times.
    VLOG(2) << "prec : "  << this->M_preconditioner_type << "\n";
    setPetscPreconditionerType( this->M_preconditioner_type,this->M_matSolverPackage_type,M_pc,is,this->worldComm(),this->name() );

#if 0
    VLOG(2) << "mat solver package : "  << this->M_matSolverPackage_type << "("  << Environment::vm()["pc-factor-mat-solver-package-type"].template as<std::string>() << ")\n";
#if 0
    //?????
    std::string type =  Environment::vm()["pc-factor-mat-solver-package-type"].template as<std::string>();
    this->setMatSolverPackageType( matSolverPackageEnumType( type ) );
#endif
    if ( option( _prefix=this->name(), _name="pc-view" ).template as<bool>() )
        check( PCView( M_pc, PETSC_VIEWER_STDOUT_WORLD ) );
#endif
    this->M_is_initialized = true;
}


template <typename T>
void PreconditionerPetsc<T>::clear ()
{
    LOG(INFO) << "PreconditionerPetsc<T>::clear\n";

    if ( this-> M_is_initialized )
    {
        this->M_is_initialized = false;

        PetscTruth is_petsc_initialized;
        PetscInitialized( &is_petsc_initialized );
        if ( is_petsc_initialized )
        {
            LOG(INFO) << "calling PCDestroy\n";
            PETSc::PCDestroy( M_pc );
        }
    }

}

template <typename T>
void PreconditionerPetsc<T>::view() const
{
    this->check( PCView( M_pc, PETSC_VIEWER_STDOUT_WORLD ) );
}



template <typename T>
void PreconditionerPetsc<T>::setPetscPreconditionerType ( const PreconditionerType & preconditioner_type,
                                                          const MatSolverPackageType & matSolverPackage_type,
                                                          PC & pc,
                                                          WorldComm const& worldComm,
                                                          std::string const& name )
{
    indexsplit_ptrtype is;
    setPetscPreconditionerType(preconditioner_type,matSolverPackage_type,pc,is,worldComm,name);
}


void
SetPCType( PC& pc, const PreconditionerType & preconditioner_type, const MatSolverPackageType & matSolverPackage_type,
           WorldComm const& worldComm )
{
    int ierr = 0;

    switch ( preconditioner_type )
    {
    case IDENTITY_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCNONE );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case CHOLESKY_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCCHOLESKY );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case ICC_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCICC );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case ILU_PRECOND:
    {
        if ( worldComm.globalSize() == 1 )
        {
            ierr = PCSetType ( pc, ( char* ) PCILU );
            CHKERRABORT( worldComm.globalComm(),ierr );
        }
#if defined(PETSC_HAVE_HYPRE) //#ifdef FEELPP_HAS_PETSC_HYPRE
        else if ( matSolverPackage_type == MATSOLVER_EUCLID )
        {
            ierr = PCSetType( pc,( char* ) PCHYPRE );
            CHKERRABORT( worldComm.globalComm(),ierr );
            ierr = PCHYPRESetType( pc, "euclid" );
            CHKERRABORT( worldComm.globalComm(),ierr );
        }
        else if ( matSolverPackage_type == MATSOLVER_PILUT )
        {
            ierr = PCSetType( pc,( char* ) PCHYPRE );
            CHKERRABORT( worldComm.globalComm(),ierr );
            ierr = PCHYPRESetType( pc, "pilut" );
            CHKERRABORT( worldComm.globalComm(),ierr );
        }
#endif
        else
        {
#if defined(PETSC_HAVE_HYPRE) //#ifdef FEELPP_HAS_PETSC_HYPRE
            ierr = PCSetType( pc,( char* ) PCHYPRE );
            CHKERRABORT( worldComm.globalComm(),ierr );
            ierr = PCHYPRESetType( pc, "euclid" );
            CHKERRABORT( worldComm.globalComm(),ierr );
#else
            // But PETSc has no truly parallel ILU, instead you have to set
            // an actual parallel preconditioner (e.g. block Jacobi) and then
            // assign ILU sub-preconditioners.
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
            ierr = PCSetType ( pc, ( char* ) PCGASM );
#else
            ierr = PCSetType ( pc, ( char* ) PCASM );
#endif
            CHKERRABORT( worldComm.globalComm(),ierr );
#endif
        }

        break;
    }

    case LU_PRECOND:
    {

        // In serial, just set the LU preconditioner type
        //if (Feel::n_processors() == 1)
        // do be changed in parallel
        if ( worldComm.globalSize() == 1 || matSolverPackage_type == MATSOLVER_MUMPS || matSolverPackage_type == MATSOLVER_PASTIX )
        {
            ierr = PCSetType ( pc, ( char* ) PCLU );
            CHKERRABORT( worldComm.globalComm(),ierr );

        }

        else
        {
            // But PETSc has no truly parallel LU, instead you have to set
            // an actual parallel preconditioner (e.g. gasm) and then
            // assign LU sub-preconditioners.
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
            ierr = PCSetType ( pc, ( char* ) PCGASM );
#else
            ierr = PCSetType ( pc, ( char* ) PCASM );
#endif
            CHKERRABORT( worldComm.globalComm(),ierr );
        }

        break;
    }

    case ASM_PRECOND:
    {
        ierr = PCSetType ( pc, ( char* ) PCASM );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;
    }
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
    case GASM_PRECOND:
    {
        ierr = PCSetType ( pc, ( char* ) PCGASM );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;
    }
#endif
    case JACOBI_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCJACOBI );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case BLOCK_JACOBI_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCBJACOBI );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case SOR_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCSOR );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case EISENSTAT_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCEISENSTAT );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

#if !(PETSC_VERSION_LESS_THAN(2,1,2))
        // Only available for PETSC >= 2.1.2
    case USER_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCMAT );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;
#endif

    case SHELL_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCSHELL );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case FIELDSPLIT_PRECOND:
        ierr = PCSetType( pc,( char* ) PCFIELDSPLIT );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case LSC_PRECOND:
        ierr = PCSetType( pc,( char* ) PCLSC );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case ML_PRECOND:
        ierr = PCSetType( pc,( char* ) PCML );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case GAMG_PRECOND:
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2)
        ierr = PCSetType( pc,( char* ) PCGAMG );
        CHKERRABORT( worldComm.globalComm(),ierr );
#else
        LOG(ERROR) << "preconditioner GAMG is available from PETSc version >= 3.2 : ";
#endif
        break;

    case BOOMERAMG_PRECOND:
#if defined(PETSC_HAVE_HYPRE) // #ifdef FEELPP_HAS_PETSC_HYPRE
        ierr = PCSetType( pc,( char* ) PCHYPRE );
        CHKERRABORT( worldComm.globalComm(),ierr );
        ierr = PCHYPRESetType( pc, "boomeramg" );
        CHKERRABORT( worldComm.globalComm(),ierr );
#else
        LOG(ERROR) << "preconditioner boomeramg is available with HYPRE package";
#endif
        break;

    case NONE_PRECOND:
        ierr = PCSetType( pc,( char* ) PCNONE );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    default:
        std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
                  << preconditioner_type       << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }

}

template <typename T>
void PreconditionerPetsc<T>::setPetscPreconditionerType ( const PreconditionerType & preconditioner_type,
                                                          const MatSolverPackageType & matSolverPackage_type,
                                                          PC & pc,
                                                          indexsplit_ptrtype const& is,
                                                          WorldComm const& worldComm,
                                                          std::string const& name )
{
    int ierr = 0;
    SetPCType( pc, preconditioner_type, matSolverPackage_type, worldComm );
    // init with petsc option if given and not interfaced
    ierr = PCSetFromOptions( pc );
    CHKERRABORT( worldComm.globalComm(),ierr );

    // configure main preconditioner
    ConfigurePC( pc, is, worldComm, "", name );

    // prepare PC to use
    ierr = PCSetUp( pc );
    CHKERRABORT( worldComm.globalComm(),ierr );

    if ( option( _prefix=name, _name="pc-view" ).template as<bool>() )
    {
        ierr = PCView( pc, PETSC_VIEWER_STDOUT_WORLD );
        CHKERRABORT( worldComm.globalComm(),ierr );
    }

}



/**
 * ConfigurePC
 */
ConfigurePC::ConfigurePC( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                          WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( worldComm, sub, prefix ),
    M_useConfigDefaultPetsc( option(_name="pc-use-config-default-petsc",_prefix=prefix,_sub=sub).as<bool>() )
{
    run( pc,is );
}
ConfigurePC::ConfigurePC( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                          WorldComm const& worldComm, std::string const& sub, std::string const& prefix,
                          std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( worldComm, sub, prefix, prefixOverwrite ),
    M_useConfigDefaultPetsc( option(_name="pc-use-config-default-petsc",_prefix=prefix,_sub=sub).as<bool>() )
{
    run( pc,is );
}

void
ConfigurePC::run( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
    VLOG(2) << "configuring PC... (sub: " << this->sub() << ")";
    google::FlushLogFiles(google::INFO);
    const char* pctype;
    this->check( PCGetType ( pc, &pctype ) );
    VLOG(2) << "configuring PC (" << this->prefix() << "." << this->sub() << ")" << pctype <<  "\n";
    google::FlushLogFiles(google::INFO);

    if ( M_useConfigDefaultPetsc )
        return;

    if ( std::string(pctype) == "gasm" )
    {
        ConfigurePCGASM( pc, is, this->worldComm(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( std::string(pctype) == "asm" )
    {
        ConfigurePCASM( pc, is, this->worldComm(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( std::string(pctype) == "bjacobi" || std::string(pctype) == "block_jacobi" )
    {
        ConfigureSubPC( pc, is, this->worldComm(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( std::string(pctype) == "lu" )
    {
        ConfigurePCLU( pc, is, this->worldComm(), this->sub(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( std::string(pctype) == "ilu" )
    {
        ConfigurePCILU( pc, this->worldComm(), this->sub(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( std::string(pctype) == "sor" )
    {
        ConfigurePCSOR( pc, this->worldComm(), this->sub(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( std::string(pctype) == "gamg" )
    {
        ConfigurePCGAMG( pc, is, this->worldComm(), this->sub(), this->prefix() );
    }
    else if ( std::string(pctype) == "ml" )
    {
        ConfigurePCML( pc, is, this->worldComm(), this->sub(), this->prefix() );
    }
    else if ( std::string(pctype) == "fieldsplit" )
    {
        ConfigurePCFieldSplit( pc, is, this->worldComm(), this->sub(), this->prefix() );
    }
    else if ( std::string(pctype) == "lsc" )
    {
        ConfigurePCLSC( pc, is, this->worldComm(), this->sub(), this->prefix() );
    }
    else if ( std::string(pctype) == "hypre" )
    {
#if defined(PETSC_HAVE_HYPRE)
        const char* hypretype;
        this->check( PCHYPREGetType( pc, &hypretype ) );
        if ( std::string( hypretype ) == "euclid" )
            ConfigurePCHYPRE_EUCLID( pc, this->worldComm(), this->sub(), this->prefix() );
#if 0
        else if ( std::string( hypretype ) == "pilut" )
            ConfigurePCHYPRE_PILUT( pc, this->worldComm(), this->sub(), this->prefix() );
        else if ( std::string( hypretype ) == "boomeramg" )
            ConfigurePCHYPRE_BOOMERAMG( pc, this->worldComm(), this->sub(), this->prefix() );
#endif
#endif
    }

    VLOG(2) << "configuring PC " << pctype << " done\n";
    google::FlushLogFiles(google::INFO);
}





/**
 * ConfigureKSP
 */
ConfigureKSP::ConfigureKSP( KSP& ksp,WorldComm const& worldComm, std::string const& sub,std::string const& prefix )
    :
    ConfigurePCBase( worldComm,sub,prefix ),
    M_type( option(_name="ksp-type",_sub=sub,_prefix=prefix).as<std::string>() ),
    M_useConfigDefaultPetsc( option(_name="ksp-use-config-default-petsc",_prefix=prefix,_sub=sub).as<bool>() ),
    M_rtol( option(_name="ksp-rtol",_sub=sub,_prefix=prefix).as<double>() ),
    M_maxit( option(_name="ksp-maxit",_sub=sub,_prefix=prefix).as<size_type>() ),
    M_showMonitor( option(_name="ksp-monitor",_sub=sub,_prefix=prefix).as<bool>() ),
    M_kspView( option(_name="ksp-view",_sub=sub,_prefix=prefix).as<bool>() ),
    M_constantNullSpace( option(_name="constant-null-space",_sub=sub,_prefix=prefix).as<bool>() )
{
    runConfigureKSP( ksp );
}
ConfigureKSP::ConfigureKSP( KSP& ksp,WorldComm const& worldComm, std::string const& sub,std::string const& prefix,
                            std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( worldComm,sub,prefix,prefixOverwrite ),
    M_type( getOption<std::string>("ksp-type",prefix,sub,prefixOverwrite) ),
    M_useConfigDefaultPetsc( getOption<bool>("ksp-use-config-default-petsc",prefix,sub,prefixOverwrite) ),
    M_rtol( getOption<double>("ksp-rtol",prefix,sub,prefixOverwrite) ),
    M_maxit( getOption<size_type>("ksp-maxit",prefix,sub,prefixOverwrite) ),
    M_showMonitor( getOption<bool>("ksp-monitor",prefix,sub,prefixOverwrite) ),
    M_kspView( getOption<bool>("ksp-view",prefix,sub,prefixOverwrite) ),
    M_constantNullSpace( getOption<bool>("constant-null-space",prefix,sub,prefixOverwrite) )
{
    runConfigureKSP( ksp );
}
void
ConfigureKSP::runConfigureKSP( KSP& ksp )
{
    // set ksp type : gmres,cg,preonly,...
    this->check( KSPSetType( ksp, M_type.c_str() ) );

    if ( M_useConfigDefaultPetsc )
        return;

    // Norm that is passed in the Krylov convergence test routines
    //this->check( KSPSetNormType( ksp, KSP_NORM_DEFAULT /*KSP_NORM_PRECONDITIONED*/ ) );

    // set ksp tolerance
    this->check( KSPSetTolerances( ksp,M_rtol,PETSC_DEFAULT,PETSC_DEFAULT,M_maxit ) );
    // monitor
    if ( M_showMonitor )
        this->check( KSPMonitorSet( ksp,KSPMonitorDefault,PETSC_NULL,PETSC_NULL ) );
    // constant null space
    if ( M_constantNullSpace )
    {
        MatNullSpace nullsp;
        this->check( MatNullSpaceCreate( PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp ) );
        this->check( KSPSetNullSpace( ksp, nullsp ) );
        PETSc::MatNullSpaceDestroy( nullsp );
    }

}



/**
 * ConfigurePCLU
 */
ConfigurePCLU::ConfigurePCLU( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                              WorldComm const& worldComm, std::string const& sub, std::string const& prefix,
                              std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( worldComm,sub,prefix,prefixOverwrite ),
    //M_matSolverPackage( option(_name="pc-factor-mat-solver-package-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>() )
    M_matSolverPackage( getOption<std::string>("pc-factor-mat-solver-package-type",prefix,sub,prefixOverwrite) )
{
    VLOG(2) << "ConfigurePC : LU\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->matSolverPackage : " << M_matSolverPackage << "\n";
    google::FlushLogFiles(google::INFO);
    runConfigurePCLU( pc );
}
void
ConfigurePCLU::runConfigurePCLU( PC& pc )
{
    // set factor package
    this->check( PCFactorSetMatSolverPackage( pc, M_matSolverPackage.c_str() ) );
}
/**
 * ConfigurePCILU
 */
ConfigurePCILU::ConfigurePCILU( PC& pc, WorldComm const& worldComm,
                                std::string const& sub, std::string const& prefix,
                                std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( worldComm,sub,prefix,prefixOverwrite ),
    //M_levels( option(_name="pc-factor-levels",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>() ),
    //M_fill( option(_name="pc-factor-fill",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<double>() )
    M_levels( getOption<int>("pc-factor-levels",prefix,sub,prefixOverwrite) ),
    M_fill( getOption<double>("pc-factor-fill",prefix,sub,prefixOverwrite) )
{
    VLOG(2) << "ConfigurePC : ILU\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->levels : " << M_levels << "\n"
            << "  |->fill : " << M_fill << "\n";
    google::FlushLogFiles(google::INFO);
    runConfigurePCILU( pc );
}
void
ConfigurePCILU::runConfigurePCILU( PC& pc )
{
    // do we need to set the mat solver package for ilu ?
    //PetscPCFactorSetMatSolverPackage( pc, "petsc" );
    this->check( PCFactorSetLevels( pc, M_levels ) );
    this->check( PCFactorSetFill( pc, M_fill ) );
}

/**
 * ConfigurePCHYPRE_EUCLID
 */
ConfigurePCHYPRE_EUCLID::ConfigurePCHYPRE_EUCLID( PC& pc,
                                                  WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( worldComm,sub,prefix ),
    M_levels( option(_name="pc-factor-levels",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>() )
{
    VLOG(2) << "ConfigurePC : HYPRE_EUCLID(ILU)\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->levels : " << M_levels << "\n";
    google::FlushLogFiles(google::INFO);
    runConfigurePCHYPRE_EUCLID( pc );
}
void
ConfigurePCHYPRE_EUCLID::runConfigurePCHYPRE_EUCLID( PC& pc )
{
#if defined(PETSC_HAVE_HYPRE)
    PCHYPRE_EUCLIDSetLevels( pc, M_levels );
#endif
}


/**
 * ConfigurePCSOR
 */
ConfigurePCSOR::ConfigurePCSOR( PC& pc, WorldComm const& worldComm,
                                std::string const& sub, std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( worldComm,sub,prefix,prefixOverwrite ),
    //M_type( option(_name="pc-sor-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>() ),
    //M_omega( option(_name="pc-sor-omega",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<double>() ),
    //M_nIteration( option(_name="pc-sor-its",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>() ),
    //M_nLocalIteration( option(_name="pc-sor-lits",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>() )
    M_type( getOption<std::string>("pc-sor-type",prefix,sub,prefixOverwrite) ),
    M_omega( getOption<double>("pc-sor-omega",prefix,sub,prefixOverwrite) ),
    M_nIteration( getOption<int>("pc-sor-its",prefix,sub,prefixOverwrite) ),
    M_nLocalIteration( getOption<int>("pc-sor-lits",prefix,sub,prefixOverwrite) )
{
    VLOG(2) << "ConfigurePC : ILU\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->omega : " << M_omega << "\n";
    google::FlushLogFiles(google::INFO);
    runConfigurePCSOR( pc );
}
void
ConfigurePCSOR::runConfigurePCSOR( PC& pc )
{
    if ( M_type == "symmetric")
        this->check( PCSORSetSymmetric( pc, SOR_SYMMETRIC_SWEEP ) );
    else if ( M_type == "forward")
        this->check( PCSORSetSymmetric( pc, SOR_FORWARD_SWEEP ) );
    else if ( M_type == "backward")
        this->check( PCSORSetSymmetric( pc, SOR_BACKWARD_SWEEP ) );
    else if ( M_type == "local_symmetric")
        this->check( PCSORSetSymmetric( pc, SOR_LOCAL_SYMMETRIC_SWEEP ) );
    else if ( M_type == "local_forward")
        this->check( PCSORSetSymmetric( pc, SOR_LOCAL_FORWARD_SWEEP ) );
    else if ( M_type == "local_backward")
        this->check( PCSORSetSymmetric( pc, SOR_LOCAL_BACKWARD_SWEEP ) );

    this->check( PCSORSetOmega( pc, M_omega ) );
    this->check( PCSORSetIterations( pc, M_nIteration, M_nLocalIteration ) );
}


/**
 * ConfigurePCGASM
 */
ConfigurePCGASM::ConfigurePCGASM( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                                  WorldComm const& worldComm, std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( worldComm,"",prefix,prefixOverwrite ),
    //M_type( option(_name="pc-gasm-type",_prefix=prefix,_worldcomm=worldComm).as<std::string>() ),
    //M_overlap( option(_name="pc-gasm-overlap",_prefix=prefix,_worldcomm=worldComm).as<int>() )
    M_type( getOption<std::string>("pc-gasm-type",prefix,"",prefixOverwrite) ),
    M_overlap( getOption<int>("pc-gasm-overlap",prefix,"",prefixOverwrite) )
{
    VLOG(2) << "ConfigurePC : GASM\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->type : " << M_type  << "\n"
            << "  |->overlap : " << M_overlap << "\n";
    google::FlushLogFiles(google::INFO);
    runConfigurePCGASM( pc,is );
}
void
ConfigurePCGASM::runConfigurePCGASM( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    if ( M_type == "restrict" ) this->check( PCGASMSetType( pc, PC_GASM_RESTRICT ) );
    else if ( M_type == "basic" ) this->check( PCGASMSetType( pc, PC_GASM_BASIC ) );
    else if ( M_type == "interpolate" ) this->check( PCGASMSetType( pc, PC_GASM_INTERPOLATE ) );
    else if ( M_type == "none" ) this->check( PCGASMSetType( pc, PC_GASM_NONE ) );
    else CHECK( false ) << "invalid gasm type : " << M_type << "\n";

    this->check( PCGASMSetOverlap( pc, M_overlap ) );
#endif
    ConfigureSubPC( pc,is,this->worldComm(),this->prefix(),this->prefixOverwrite() );
}

/**
 * ConfigurePCGASM
 */
ConfigurePCASM::ConfigurePCASM( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                                WorldComm const& worldComm, std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( worldComm,"",prefix,prefixOverwrite ),
    //M_type( option(_name="pc-asm-type",_prefix=prefix,_worldcomm=worldComm).as<std::string>() ),
    //M_overlap( option(_name="pc-asm-overlap",_prefix=prefix,_worldcomm=worldComm).as<int>() )
    M_type( getOption<std::string>("pc-asm-type",prefix,"",prefixOverwrite) ),
    M_overlap( getOption<int>("pc-asm-overlap",prefix,"",prefixOverwrite) )
{
    VLOG(2) << "ConfigurePC : ASM\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->type : " << M_type  << "\n"
            << "  |->overlap : " << M_overlap << "\n";
    google::FlushLogFiles(google::INFO);
    runConfigurePCASM( pc,is );
}
void
ConfigurePCASM::runConfigurePCASM( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
    if ( M_type == "restrict" ) this->check( PCASMSetType( pc, PC_ASM_RESTRICT ) );
    else if ( M_type == "basic" ) this->check( PCASMSetType( pc, PC_ASM_BASIC ) );
    else if ( M_type == "interpolate" ) this->check( PCASMSetType( pc, PC_ASM_INTERPOLATE ) );
    else if ( M_type == "none" ) this->check( PCASMSetType( pc, PC_ASM_NONE ) );
    else CHECK( false ) << "invalid asm type : " << M_type << "\n";

    this->check( PCASMSetOverlap( pc, M_overlap ) );

    ConfigureSubPC( pc,is,this->worldComm(),this->prefix(),this->prefixOverwrite() );
}

/**
 * ConfigureSubPC
 */
ConfigureSubPC::ConfigureSubPC( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                                WorldComm const& worldComm, std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( worldComm,"",prefix, prefixOverwrite ),
    //M_subPCtype( option(_name="pc-type",_sub="sub",_prefix=prefix).as<std::string>() ),
    //M_subPCview( option(_name="pc-view",_sub="sub",_prefix=prefix).as<bool>() ),
    M_subPCtype( getOption<std::string>("pc-type",prefix,"sub",prefixOverwrite) ),
    M_subPCview( getOption<bool>("pc-view",prefix,"sub",prefixOverwrite) ),
    M_nBlock(0)
{
    this->check( PCSetUp( pc ) );
#if PETSC_VERSION_LESS_THAN(3,4,0)
    const PCType thepctype;
#else
    PCType thepctype;
#endif
    this->check( PCGetType( pc, &thepctype ) );
    M_subPCfromPCtype = std::string( thepctype );

    // To store array of local KSP contexts on this processor
    KSP* subksps;
    if ( M_subPCfromPCtype == "block_jacobi" || M_subPCfromPCtype == "bjacobi" )
        this->check( PCBJacobiGetSubKSP( pc, &M_nBlock, PETSC_NULL, &subksps ) );
    else if ( M_subPCfromPCtype == "asm" )
        this->check( PCASMGetSubKSP( pc, &M_nBlock, PETSC_NULL, &subksps ) );
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
    else if ( M_subPCfromPCtype == "gasm" )
        this->check( PCGASMGetSubKSP( pc, &M_nBlock, PETSC_NULL, &subksps ) );
#endif
    else CHECK( false ) << "invalid pctype " << M_subPCfromPCtype << "\n";

    VLOG(2) << "ConfigureSubPC : from "<< M_subPCfromPCtype <<"\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->nBlock    : " << M_nBlock << "\n"
            << "  |->subPCtype : " << M_subPCtype  << "\n";
    google::FlushLogFiles(google::INFO);

    for ( int i=0; i<M_nBlock; ++i )
    {
        runConfigureSubPC( subksps[i],is );
    }
}
void
ConfigureSubPC::runConfigureSubPC( KSP& ksp, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
    // configure coarse ksp
    //this->check( KSPSetFromOptions( ksp ) );
    ConfigureKSP kspConf( ksp, this->worldComm(), "sub", this->prefix(), this->prefixOverwrite() );
    this->check( KSPSetUp( ksp ) );

    PC subpc;
    // Get pointer to sub KSP object's PC
    this->check( KSPGetPC( ksp, &subpc ) );

    // configure sub-pc
    this->check( PCSetType( subpc, M_subPCtype.c_str() ) );
    ConfigurePC( subpc, is, this->worldComm(), "sub", this->prefix(), this->prefixOverwrite() );

    this->check( PCSetUp( subpc ) );

    if ( kspConf.kspView() )
        this->check( KSPView( ksp, PETSC_VIEWER_STDOUT_SELF ) );
    else if ( M_subPCview )
        this->check( PCView( subpc, PETSC_VIEWER_STDOUT_SELF ) );
}

/**
 * ConfigurePCML
 */
ConfigurePCML::ConfigurePCML( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                              WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( worldComm,sub,prefix ),
    M_mgType( option(_name="pc-mg-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>() ),
    M_nLevels( option(_name="pc-mg-levels",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>() ),
    M_mlReuseInterp( option(_name="pc-ml-reuse-interpolation",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<bool>() ),
    M_mlKeepAggInfo( option(_name="pc-ml-keep-agg-info",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<bool>() ),
    M_mlReusable( option(_name="pc-ml-reusable",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<bool>() ),
    M_mlOldHierarchy( option(_name="pc-ml-old-hierarchy",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<bool>() ),
    M_prefixMGCoarse( (boost::format( "%1%%2%mg-coarse" ) %prefixvm( prefix,"" ) %std::string((sub.empty())?"":sub+"-")  ).str() ),
    M_coarsePCtype( option(_name="pc-type",_prefix=M_prefixMGCoarse).as<std::string>() ),
    M_coarsePCview( option(_name="pc-view",_prefix=M_prefixMGCoarse).as<bool>() )
{
    VLOG(2) << "ConfigurePC : ML\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->mgType : " << M_mgType << "\n"
            << "  |->maxLevels : " << M_nLevels << "\n";
    google::FlushLogFiles(google::INFO);
    runConfigurePCML( pc,is );
}
void
ConfigurePCML::runConfigurePCML( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
#if 0
    // Sets the number of levels to use with MG.
    // Must be called before any other MG routine
    this->check( PCMGSetLevels( pc, M_nLevels, PETSC_NULL) );
    std::cout << " M_mgType " << M_mgType << "\n";
    if ( M_mgType=="multiplicative" ) this->check( PCMGSetType( pc, PC_MG_MULTIPLICATIVE ) );
    if ( M_mgType=="additive" ) this->check( PCMGSetType( pc, PC_MG_ADDITIVE ) );
    if ( M_mgType=="full" ) this->check( PCMGSetType( pc, PC_MG_FULL ) );
    if ( M_mgType=="kaskade" ) this->check( PCMGSetType( pc, PC_MG_KASKADE ) );
#endif
#if 0
    // warning, this function (seems) create 2 smoother up and down
    int smoothdown= option(_name="pc-mg-smoothdown",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>();
    ierr = PCMGSetNumberSmoothDown( pc, smoothdown );
    CHKERRABORT( worldComm.globalComm(),ierr );
#endif

#if defined(PETSC_HAVE_ML)
    PCMLSetMaxNlevels( pc, M_nLevels );
    PCMLSetReuseInterpolation( pc, static_cast<PetscBool>( M_mlReuseInterp ) );
    PCMLSetKeepAggInfo( pc, static_cast<PetscBool>( M_mlKeepAggInfo ) );
    PCMLSetReusable( pc, static_cast<PetscBool>( M_mlReusable ) );
    PCMLSetOldHierarchy( pc, static_cast<PetscBool>( M_mlOldHierarchy ) );
#endif

    // configure coarse solver
    configurePCMLCoarse( pc,is );
    // configure fine solvers
    ConfigurePCMGLevels( pc, is, this->worldComm(), this->sub(), this->prefix() );

    // must be called after setup ml pc because this one PCMGSetLevels and reset mg prec associated
    if ( M_mgType=="multiplicative" ) this->check( PCMGSetType( pc, PC_MG_MULTIPLICATIVE ) );
    else if ( M_mgType=="additive" ) this->check( PCMGSetType( pc, PC_MG_ADDITIVE ) );
    else if ( M_mgType=="full" ) this->check( PCMGSetType( pc, PC_MG_FULL ) );
    else if ( M_mgType=="kaskade" ) this->check( PCMGSetType( pc, PC_MG_KASKADE ) );
    else CHECK( false ) << "invalid mgType :" << M_mgType << "\n";

}
void
ConfigurePCML::configurePCMLCoarse( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
    this->check( PCSetUp( pc ) );
    // get coarse-ksp
    KSP coarseksp;
    this->check( PCMGGetCoarseSolve( pc, &coarseksp) );

    // configure coarse ksp
    this->check( KSPSetFromOptions( coarseksp ) );
    ConfigureKSP kspConf( coarseksp, this->worldComm(), "", M_prefixMGCoarse );
    this->check( KSPSetUp( coarseksp ) );
    // get coarse pc
    PC coarsepc;
    this->check( KSPGetPC( coarseksp, &coarsepc ) );
    // configure coarse pc
    this->check( PCSetFromOptions( coarsepc ) );
    this->check( PCSetType( coarsepc, M_coarsePCtype.c_str() ) );
    ConfigurePC( coarsepc, is, this->worldComm(), "", M_prefixMGCoarse );
    // setup pc (all do here because the setup of ml has not the same effect that classic setup)
    //this->check( PCSetUp( pc ) );

    //this->check( KSPSetUp( coarseksp ) );

#if 0
    // configure coarse ksp
    this->check( KSPSetFromOptions( coarseksp ) );
    ConfigureKSP kspConf( coarseksp, this->worldComm(), "", M_prefixMGCoarse );
    // setup coarse ksp
    this->check( KSPSetUp( coarseksp ) );
#endif
    // setup coarse pc
    this->check( PCSetUp( coarsepc ) );

    PetscViewer viewer = (this->sub().empty())? PETSC_VIEWER_STDOUT_WORLD : PETSC_VIEWER_STDOUT_SELF;
    if ( kspConf.kspView() )
        this->check( KSPView( coarseksp, viewer ) );
    else if ( M_coarsePCview )
        this->check( PCView( coarsepc, viewer ) );

}

/**
 * ConfigurePCGAMG
 */
ConfigurePCGAMG::ConfigurePCGAMG( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                                  WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( worldComm,sub,prefix ),
    M_mgType( option(_name="pc-mg-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>() ),
    M_gamgType( option(_name="pc-gamg-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>() ),
    M_nLevels( option(_name="pc-mg-levels",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>() ),
    M_procEqLim( option(_name="pc-gamg-proc-eq-lim",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>() ),
    M_coarseEqLim(option(_name="pc-gamg-coarse-eq-lim",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>() ),
    M_threshold( option(_name="pc-gamg-threshold",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<double>() ),
    M_prefixMGCoarse( (boost::format( "%1%%2%mg-coarse" ) %prefixvm( prefix,"" ) %std::string((sub.empty())?"":sub+"-")  ).str() ),
    M_coarsePCtype( option(_name="pc-type",_prefix=M_prefixMGCoarse).as<std::string>() ),
    M_coarsePCview( option(_name="pc-view",_prefix=M_prefixMGCoarse).as<bool>() )
{
    VLOG(2) << "ConfigurePC : GAMG\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->mgType : " << M_mgType << "\n"
            << "  |->maxLevels : " << M_nLevels << "\n";
    google::FlushLogFiles(google::INFO);
    runConfigurePCGAMG( pc, is );
}

void
ConfigurePCGAMG::runConfigurePCGAMG( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
    this->check( PCGAMGSetType( pc, M_gamgType.c_str() ) );
#endif
    // PCSetFromOptions is called here because PCGAMGSetType destroy all unless the type_name
    this->check( PCSetFromOptions( pc ) );
#if 0
    // Sets the number of levels to use with MG.
    // Must be called before any other MG routine
    this->check( PCMGSetLevels( pc, M_nLevels, PETSC_NULL) );

    if ( M_mgType=="multiplicative" ) this->check( PCMGSetType( pc, PC_MG_MULTIPLICATIVE ) );
    if ( M_mgType=="additive" ) this->check( PCMGSetType( pc, PC_MG_ADDITIVE ) );
    if ( M_mgType=="full" ) this->check( PCMGSetType( pc, PC_MG_FULL ) );
    if ( M_mgType=="kaskade" ) this->check( PCMGSetType( pc, PC_MG_KASKADE ) );
#endif
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
    //
    this->check( PCGAMGSetNlevels( pc,M_nLevels ) );
    // Set number of equations to aim for on coarse grids via processor reduction
    this->check( PCGAMGSetProcEqLim( pc, M_procEqLim ) );
    // Set max number of equations on coarse grids
    this->check( PCGAMGSetProcEqLim( pc, M_coarseEqLim ) );
    // Relative threshold to use for dropping edges in aggregation graph
    this->check( PCGAMGSetThreshold( pc, M_threshold ) );
#endif
    // setup sub-pc
    this->check( PCSetUp( pc ) );
    //this->check( PCView( pc, PETSC_VIEWER_STDOUT_WORLD ) );

    // configure coarse pc
    configurePCGAMGCoarse( pc,is );
    // configure level pc
    ConfigurePCMGLevels( pc, is, this->worldComm(), this->sub(), this->prefix() );

    // must be called after setup gamg pc because this one call PCMGSetLevels and reset mg prec associated
    if ( M_mgType=="multiplicative" ) this->check( PCMGSetType( pc, PC_MG_MULTIPLICATIVE ) );
    else if ( M_mgType=="additive" ) this->check( PCMGSetType( pc, PC_MG_ADDITIVE ) );
    else if ( M_mgType=="full" ) this->check( PCMGSetType( pc, PC_MG_FULL ) );
    else if ( M_mgType=="kaskade" ) this->check( PCMGSetType( pc, PC_MG_KASKADE ) );
    else CHECK( false ) << "invalid mgType :" << M_mgType << "\n";

}

void
ConfigurePCGAMG::configurePCGAMGCoarse( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
    // get coarse-ksp
    KSP coarseksp;
    this->check( PCMGGetCoarseSolve( pc, &coarseksp) );
    // configure coarse ksp
    this->check( KSPSetFromOptions( coarseksp ) );
    ConfigureKSP kspConf( coarseksp, this->worldComm(), "", M_prefixMGCoarse );
    // setup coarse ksp
    this->check( KSPSetUp( coarseksp ) );

    // get coarse pc
    PC coarsepc;
    this->check( KSPGetPC( coarseksp, &coarsepc ) );
    // configure coarse pc
    this->check( PCSetFromOptions( coarsepc ) );
    this->check( PCSetType( coarsepc, M_coarsePCtype.c_str() ) );
    ConfigurePC( coarsepc, is, this->worldComm(), "", M_prefixMGCoarse );
    // setup coarse pc
    this->check( PCSetUp( coarsepc ) );

    PetscViewer viewer = (this->sub().empty())? PETSC_VIEWER_STDOUT_WORLD : PETSC_VIEWER_STDOUT_SELF;
    if ( kspConf.kspView() )
        this->check( KSPView( coarseksp, viewer ) );
    else if ( M_coarsePCview )
        this->check( PCView( coarsepc, viewer ) );
}

/**
 * ConfigurePCMGLevels
 */
ConfigurePCMGLevels::ConfigurePCMGLevels( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                                          WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( worldComm,sub,prefix )
{
    //this->check( PCSetUp( pc ) );
    this->check( PCMGGetLevels( pc, &M_nLevels) );
    //std::cout << "M_nLevels " << M_nLevels << std::endl;

    M_prefixMGLevels.resize( M_nLevels-1 );
    M_mgLevelsPCtype.resize( M_nLevels-1 );
    M_mgLevelsPCview.resize( M_nLevels-1 );
    M_mgLevelsKSPview.resize( M_nLevels-1 );
    M_mgLevelsMatSolverPackage.resize( M_nLevels-1 );
    std::string mgctx = (sub.empty())? "mg-" : sub+"-mg-";

    // get generic option for all levels
    std::string prefixAllLevel = ( boost::format( "%1%%2%levels" ) %prefixvm( this->prefix(),"" ) %mgctx ).str();
    //M_mgLevelsKSPview[0] = option(_name="ksp-view",_prefix=prefixAllLevel).as<bool>();
    M_mgLevelsPCtype[0] = option(_name="pc-type",_prefix=prefixAllLevel).as<std::string>();
    M_mgLevelsPCview[0] = option(_name="pc-view",_prefix=prefixAllLevel).as<bool>();
    M_mgLevelsMatSolverPackage[0] = option(_name="pc-factor-mat-solver-package-type",_prefix=prefixAllLevel).as<std::string>();
    for ( int level=2; level<M_nLevels; ++level )
    {
        //std::string prefixAllLevel = ( boost::format( "%1%%2%levels" ) %prefixvm( this->prefix(),"" ) %mgctx ).str();
        //M_mgLevelsKSPview[level-1] = M_mgLevelsKSPview[0];
        M_prefixMGLevels[level-1] = prefixAllLevel;
        M_mgLevelsPCtype[level-1] = M_mgLevelsPCtype[0];
        M_mgLevelsPCview[level-1] = M_mgLevelsPCview[0];
        M_mgLevelsMatSolverPackage[level-1] = M_mgLevelsMatSolverPackage[0];
    }
    // overwrite specific options for each level < 5 (if given of course)
    for ( int level=1; level<std::min(M_nLevels,6); ++level )
    {
        std::string prefixCurLevel = ( boost::format( "%1%%2%levels%3%" ) %prefixvm( this->prefix(),"" ) %mgctx %level ).str();
        M_prefixMGLevels[level-1] = prefixCurLevel;
        if ( Environment::vm().count( prefixvm(prefixCurLevel,"ksp-view") ) )
            M_mgLevelsKSPview[level-1] = option(_name="ksp-view",_prefix=prefixCurLevel).as<bool>();
        if ( Environment::vm().count( prefixvm(prefixCurLevel,"pc-type") ) )
            M_mgLevelsPCtype[level-1] = option(_name="pc-type",_prefix=prefixCurLevel).as<std::string>();
        if ( Environment::vm().count( prefixvm(prefixCurLevel,"pc-view") ) )
            M_mgLevelsPCview[level-1] = option(_name="pc-view",_prefix=prefixCurLevel).as<bool>();
        if ( Environment::vm().count( prefixvm(prefixCurLevel,"pc-factor-mat-solver-package-type") ) )
            M_mgLevelsMatSolverPackage[level-1] = option(_name="pc-factor-mat-solver-package-type",_prefix=prefixCurLevel).as<std::string>();
    }
    // overwrite options for fine level
    std::string prefixFineLevel = ( boost::format( "%1%%2%fine-level" ) %prefixvm( this->prefix(),"" ) %mgctx ).str();
    //M_prefixMGLevels[M_nLevel-2] = prefixFineLevel;
    if ( Environment::vm().count( prefixvm(prefixFineLevel,"ksp-view") ) )
        M_mgLevelsKSPview[M_nLevels-2] = option(_name="ksp-view",_prefix=prefixFineLevel).as<bool>();
    if ( Environment::vm().count( prefixvm(prefixFineLevel,"pc-type") ) )
        M_mgLevelsPCtype[M_nLevels-2] = option(_name="pc-type",_prefix=prefixFineLevel).as<std::string>();
    if ( Environment::vm().count( prefixvm(prefixFineLevel,"pc-view") ) )
        M_mgLevelsPCview[M_nLevels-2] = option(_name="pc-view",_prefix=prefixFineLevel).as<bool>();
    if ( Environment::vm().count( prefixvm(prefixFineLevel,"pc-factor-mat-solver-package-type") ) )
        M_mgLevelsMatSolverPackage[M_nLevels-2] = option(_name="pc-factor-mat-solver-package-type",_prefix=prefixFineLevel).as<std::string>();

    // configure each levels
    for ( int level=1; level<M_nLevels; ++level )
        runConfigurePCMGLevels( pc,is,level );
}

void
ConfigurePCMGLevels::runConfigurePCMGLevels( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is, int level )
{
    std::string prefixCurLevel = M_prefixMGLevels[level-1];
    std::string mgctx = (this->sub().empty())? "mg-" : this->sub()+"-mg-";
    std::string prefixAllLevel = ( boost::format( "%1%%2%levels" ) %prefixvm( this->prefix(),"" ) %mgctx ).str();
    std::string prefixFineLevel = ( boost::format( "%1%%2%fine-level" ) %prefixvm( this->prefix(),"" ) %mgctx ).str();

    std::vector<std::string> prefixLevelOverwrite;
    if ( level<6 )
        prefixLevelOverwrite.push_back( prefixCurLevel );
    if ( level == (M_nLevels-1) )
        prefixLevelOverwrite.push_back( prefixFineLevel );
    //-------------------------------------------------------------------//
    KSP levelksp;
    // get ksp
    this->check( PCMGGetSmoother( pc, level, &levelksp ) );
    // init ksp from option
    this->check( KSPSetFromOptions( levelksp ) );

    ConfigureKSP kspConf( levelksp,this->worldComm(), "", prefixAllLevel, prefixLevelOverwrite );

#if 0
    // warning : use KSP_NORM_PRECONDITIONED and force convergence
    this->check( KSPSetNormType( levelksp, KSP_NORM_PRECONDITIONED ) );
    void *cctx;
    this->check( KSPDefaultConvergedCreate(&cctx) );
    this->check( KSPSetConvergenceTest( levelksp, KSPDefaultConverged, cctx, PETSC_NULL ) );
#endif
    // setup coarse ksp
    this->check( KSPSetUp( levelksp ) );

    //-------------------------------------------------------------------//
    // get level pc
    PC levelpc;
    this->check( KSPGetPC( levelksp, &levelpc ) );
    // configure level pc
    this->check( PCSetFromOptions( levelpc ) );
    this->check( PCSetType( levelpc, M_mgLevelsPCtype[level-1].c_str() ) );

    SetPCType( levelpc, pcTypeConvertStrToEnum( M_mgLevelsPCtype[level-1] ),
               matSolverPackageConvertStrToEnum( M_mgLevelsMatSolverPackage[level-1] ),
               this->worldComm() );
    ConfigurePC( levelpc, is, this->worldComm(), "", prefixAllLevel , prefixLevelOverwrite );
    // setup level pc
    this->check( PCSetUp( levelpc ) );

    PetscViewer viewer = (this->sub().empty())? PETSC_VIEWER_STDOUT_WORLD : PETSC_VIEWER_STDOUT_SELF;
    if ( kspConf.kspView() )
        this->check( KSPView( levelksp, viewer ) );
    else if ( M_mgLevelsPCview[level-1] )
        this->check( PCView( levelpc, viewer ) );
}

/**
 * ConfigurePCFieldSplit
 */
ConfigurePCFieldSplit::ConfigurePCFieldSplit( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                                              WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( worldComm,sub,prefix ),
    M_type( option(_name="fieldsplit-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>() ),
    M_schurFactType( option(_name="fieldsplit-schur-fact-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>() ),
    M_schurPrecond( option(_name="fieldsplit-schur-precondition",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>() )
{
    VLOG(2) << "ConfigurePC : FieldSplit\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->type : " << M_type << "\n";
    runConfigurePCFieldSplit( pc, is );
}
void
ConfigurePCFieldSplit::runConfigurePCFieldSplit( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
    /*const PetscInt ufields[] = {0,2},pfields[] = {1};
      this->check( PCFieldSplitSetFields( pc , NULL, 2, ufields,ufields) );
      this->check( PCFieldSplitSetFields( pc , NULL, 1, pfields,pfields) );*/

    PCCompositeType theFieldSplitType = PC_COMPOSITE_SCHUR;
    if ( M_type == "schur" ) theFieldSplitType = PC_COMPOSITE_SCHUR;
    else if ( M_type == "additive" ) theFieldSplitType = PC_COMPOSITE_ADDITIVE;
    else if ( M_type == "multiplicative" ) theFieldSplitType = PC_COMPOSITE_MULTIPLICATIVE;
    else if ( M_type == "symmetric-multiplicative" ) theFieldSplitType = PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE;
    else if ( M_type == "special" ) theFieldSplitType = PC_COMPOSITE_SPECIAL;
    this->check( PCFieldSplitSetType( pc, theFieldSplitType ) );

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
    if ( M_type == "schur" )
    {
        PCFieldSplitSchurFactType theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_FULL;
        if ( M_schurFactType == "diag")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_DIAG;
        else if ( M_schurFactType == "lower")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_LOWER;
        else if ( M_schurFactType == "upper")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_UPPER;
        else if ( M_schurFactType == "full")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_FULL;
        this->check( PCFieldSplitSetSchurFactType( pc,theSchurFactType ) );

        PCFieldSplitSchurPreType theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_SELF;
        if ( M_schurPrecond == "self")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_SELF;
        else if ( M_schurPrecond == "user")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_USER;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,4,0 )
        else if ( M_schurPrecond == "a11")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_A11;
#else
        else if ( M_schurPrecond == "a11")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_DIAG;
#endif
        this->check( PCFieldSplitSchurPrecondition( pc, theSchurPrecond, NULL ) );
    }
#endif

    // config sub ksp/pc
    ConfigurePCFieldSplit::ConfigureSubKSP( pc,is,this->worldComm(),this->sub(),this->prefix() );
}

/**
 * ConfigurePCFieldSplitSubKSP
 */
ConfigurePCFieldSplit::ConfigureSubKSP::ConfigureSubKSP( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                                                         WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( worldComm,sub,prefix ),
    M_nSplit(0)
{
    // call necessary before PCFieldSplitGetSubKSP
    this->check( PCSetUp( pc ) );

    // To store array of local KSP contexts on this processor
    KSP* subksps;
    this->check( PCFieldSplitGetSubKSP(pc,&M_nSplit,&subksps ) );

    M_prefixSplit.resize(M_nSplit);
    M_subPCview.resize(M_nSplit);
    M_subPCtype.resize(M_nSplit);
    M_subMatSolverPackage.resize(M_nSplit);
    for ( int i=0; i<M_nSplit; ++i )
    {
        std::string prefixSplit = prefixvm(this->prefix() , (boost::format( "fieldsplit-%1%" )  %i ).str() );
        M_prefixSplit[i] = prefixSplit;
        M_subPCview[i] = option(_name="pc-view",_prefix=prefixSplit).as<bool>();
        M_subPCtype[i] = option(_name="pc-type",_prefix=prefixSplit).as<std::string>();
        M_subMatSolverPackage[i] = option(_name="pc-factor-mat-solver-package-type",_prefix=prefixSplit).as<std::string>();
    }

    // Loop over sub-ksp objects
    for ( int splitId=0; splitId<M_nSplit; ++splitId )
    {
        VLOG(2) << "configure split " << splitId << " with prefix "<< M_prefixSplit[splitId] << "\n";
        google::FlushLogFiles(google::INFO);

        runConfigureSubKSP( subksps[splitId], is, splitId );
    }
}

void
ConfigurePCFieldSplit::ConfigureSubKSP::runConfigureSubKSP(KSP& ksp, PreconditionerPetsc<double>::indexsplit_ptrtype const& is, int splitId )
{
    std::string prefixSplit = M_prefixSplit[splitId];

#if 0
    Mat A00;Mat A01;Mat A10; Mat A11;
    this->check( PCFieldSplitGetSchurBlocks(pc,&A00,&A01,&A10, &A11) );
    if (i==0)
        this->check( KSPSetOperators( ksp, A00, A00,
                                      PetscGetMatStructureEnum(MatrixStructure::SAME_PRECONDITIONER)) );
#endif

    ConfigureKSP kspConf( ksp, this->worldComm(), "", prefixSplit );
    // setup ksp
    this->check( KSPSetUp( ksp ) );

    PC subpc;
    // get sub-pc
    this->check( KSPGetPC( ksp, &subpc ) );

    // configure sub PC
    //this->check( PCSetType( subpc, M_subPCtype[splitId].c_str() ) );
    SetPCType( subpc, pcTypeConvertStrToEnum( M_subPCtype[splitId] ),
               matSolverPackageConvertStrToEnum( M_subMatSolverPackage[splitId] ),
               this->worldComm() );

    // init pc from specific default option
    this->check( PCSetFromOptions( subpc ) );

    // in case of fieldsplit in fieldsplit, need to pass the is corresponding
    if ( M_subPCtype[splitId] == "fieldsplit" )
    {
        CHECK( is ) << "index split is not initialized\n";

        std::string fieldsDefStr = option(_name="fieldsplit-fields",_prefix=prefixSplit).as<std::string>();
        auto fieldsDef = IndexSplit::parseFieldsDef( fieldsDefStr );
        //std::cout << "fieldsplit fieldsDefStr " << fieldsDefStr << "\n";
        //auto isUsed = is.applyFieldsDef( IndexSplit::FieldsDef(  {  { 0 , { 0 } }, { 1 , { 2 } } } ) );
        //is->showMe();
        auto isUsed = is->applyFieldsDef( IndexSplit::FieldsDef( fieldsDef ) );
        //isUsed->showMe();

        std::vector<IS> isPetsc;
        PetscConvertIndexSplit( isPetsc ,*isUsed,this->worldComm());
        for ( int i = 0 ; i < isPetsc.size(); ++i )
        {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
            std::ostringstream os; os << i;
            this->check( PCFieldSplitSetIS( subpc,os.str().c_str(),isPetsc[i] ) );
#else
            this->check( PCFieldSplitSetIS( subpc,isPetsc[i] ) );
#endif
        }
        //TODO : maybe delete or store  isPetsc
    }

    if ( M_subPCtype[splitId] == "lsc" )
        CHECK( splitId==1 ) << "lsc must be use with only field 1, not " << splitId << "\n";


    // configure sub-pc
    ConfigurePC( subpc, is, this->worldComm(), "", prefixSplit );

    // setup sub-pc
    this->check( PCSetUp( subpc ) );

    //PetscViewer viewer = (this->sub().empty())? PETSC_VIEWER_STDOUT_WORLD : PETSC_VIEWER_STDOUT_SELF;
    PetscViewer viewer = PETSC_VIEWER_STDOUT_WORLD;
    if ( kspConf.kspView() )
        this->check( KSPView( ksp, viewer ) );
    else if ( M_subPCview[splitId] )
        this->check( PCView( subpc, viewer ) );
}


/**
 * ConfigurePCLSC
 */
ConfigurePCLSC::ConfigurePCLSC( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                                WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( worldComm,sub,prefix ),
    M_prefixLSC( prefixvm(this->prefix(),"lsc") ),
    M_subPCtype( option(_name="pc-type",_prefix=M_prefixLSC).as<std::string>() ),
    M_subMatSolverPackage( option(_name="pc-factor-mat-solver-package-type",_prefix=M_prefixLSC).as<std::string>() ),
    M_subPCview( option(_name="pc-view",_prefix=M_prefixLSC).as<bool>() )
{
    VLOG(2) << "ConfigurePC : LSC\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->prefixLSC : " << M_prefixLSC  << "\n"
            << "  |->subPCtype : " << M_subPCtype << "\n";
    google::FlushLogFiles(google::INFO);
    runConfigurePCLSC( pc, is );
}
void
ConfigurePCLSC::runConfigurePCLSC( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
    // setup sub-pc
    this->check( PCSetUp( pc ) );
    //-----------------------------------------------------------//
    // get sub-ksp
    KSP subksp;
    PCLSCGetKSP( pc, subksp );
    // configure sub-ksp
    this->check( KSPSetFromOptions( subksp ) );
    ConfigureKSP kspConf( subksp, this->worldComm(), "", M_prefixLSC );
    // setup sub-ksp
    this->check( KSPSetUp( subksp ) );
    //-----------------------------------------------------------//
    // get sub-pc
    PC subpc;
    this->check( KSPGetPC( subksp, &subpc ) );
    // configure sub-pc
    //this->check( PCSetType( subpc, M_subPCtype.c_str() ) );
    SetPCType( subpc, pcTypeConvertStrToEnum( M_subPCtype ),
               matSolverPackageConvertStrToEnum( M_subMatSolverPackage ),
               this->worldComm() );

    this->check( PCSetFromOptions( subpc ) );
    ConfigurePC( subpc, is, this->worldComm(), "", M_prefixLSC );
    // setup sub-pc
    this->check( PCSetUp( subpc ) );
    //-----------------------------------------------------------//
    // ksp and pc view
    if ( kspConf.kspView() )
        this->check( KSPView( subksp, PETSC_VIEWER_STDOUT_WORLD ) );
    else if ( M_subPCview )
        this->check( PCView( subpc, PETSC_VIEWER_STDOUT_WORLD ) );
}




#if 0
void
configurePCWithPetscCommandLineOption( std::string prefixFeelBase, std::string prefixPetscBase )
{
    // For catching PETSc error return codes
    int ierr = 0;

    std::string pctype = option(_name="pc-type",_prefix=prefixFeelBase ).as<std::string>();

    std::string option_pc_type = "-"+prefixPetscBase+"_pc_type";
    ierr = PetscOptionsClearValue( option_pc_type.c_str() );
    ierr = PetscOptionsInsertString( (option_pc_type+" "+pctype).c_str() );
    VLOG(2) << " configurePCWithPetscCommandLineOption with "<< option_pc_type << " "<< pctype << "\n";
    if ( pctype=="lu" )
    {
        std::string PCFMSPtype =  option(_name="pc-factor-mat-solver-package-type",_prefix=prefixFeelBase).as<std::string>();
        std::string option_pc_factor_mat_solver_package = "-"+prefixPetscBase+"_pc_factor_mat_solver_package";
        ierr = PetscOptionsClearValue( option_pc_factor_mat_solver_package.c_str() );
        ierr = PetscOptionsInsertString( (option_pc_factor_mat_solver_package+" "+PCFMSPtype).c_str() );
    }
    else if ( pctype=="gasm" )
    {
        int gasmoverlap = option(_name="pc-gasm-overlap",_prefix=prefixFeelBase).as<int>();

        std::string option_pc_gasm_overlap = "-"+prefixPetscBase+"_pc_gasm_overlap";
        std::string optionval_pc_gasm_overlap = ( boost::format( "%1% %2%") %option_pc_gasm_overlap %gasmoverlap ).str();
        ierr = PetscOptionsClearValue( option_pc_gasm_overlap.c_str() );
        ierr = PetscOptionsInsertString( optionval_pc_gasm_overlap.c_str() );

        std::string option_sub_pc_type = "-"+prefixPetscBase+"_sub_pc_type";
        std::string subpctype =  option(_name="pc-type",_sub="sub",_prefix=prefixFeelBase).as<std::string>();
        ierr = PetscOptionsClearValue( option_sub_pc_type.c_str() );
        ierr = PetscOptionsInsertString( (option_sub_pc_type+" "+subpctype).c_str() );

        if (subpctype=="lu")
        {
            std::string option_sub_pc_factor_mat_solver_package = "-"+prefixPetscBase+"_sub_pc_factor_mat_solver_package";
            std::string t = option(_name="pc-factor-mat-solver-package-type",_sub="sub",_prefix=prefixFeelBase).as<std::string>();
            ierr = PetscOptionsClearValue( option_sub_pc_factor_mat_solver_package.c_str() );
            ierr = PetscOptionsInsertString( (option_sub_pc_factor_mat_solver_package+" "+t).c_str() );
        }
    }

}
#endif



template class PreconditionerPetsc<double>;

}
