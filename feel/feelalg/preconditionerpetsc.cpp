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
        int ierr = PCCreate( this->worldComm().globalComm(),&M_pc );

        CHKERRABORT( this->worldComm().globalComm(),ierr );
        ierr = PCSetFromOptions ( M_pc );
        CHKERRABORT( this->worldComm().globalComm(),ierr );
#if PETSC_VERSION_LESS_THAN(3,4,0)
        const PCType pc_type;
#else
        PCType pc_type;
#endif
        ierr = PCGetType ( M_pc, &pc_type );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        LOG(INFO) << "preconditionerpetsc set as "  << pc_type << "\n";
        MatrixPetsc<T> * pmatrix = dynamic_cast<MatrixPetsc<T>*>( this->M_matrix.get() );

        M_mat = pmatrix->mat();

        if (this->M_preconditioner_type==FIELDSPLIT_PRECOND )
        {
            ierr = PCSetType( M_pc,( char* ) PCFIELDSPLIT );
            CHKERRABORT( this->worldComm(),ierr );

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
            int ierr = PCSetType( M_pc,( char* ) PCFIELDSPLIT );
            CHKERRABORT( this->worldComm(),ierr );
            pmatrix->updatePCFieldSplit( M_pc );
            is=pmatrix->indexSplit();
        }
        this->M_mat_has_changed = false;
    }

    //int ierr = PCSetOperators( M_pc,M_mat,M_mat, PetscGetMatStructureEnum(MatrixStructure::SAME_NONZERO_PATTERN) );
    //int ierr = PCSetOperators( M_pc,M_mat,M_mat, PetscGetMatStructureEnum(MatrixStructure::DIFFERENT_NONZERO_PATTERN) );
    int ierr = PCSetOperators( M_pc,M_mat,M_mat, PetscGetMatStructureEnum(this->M_prec_matrix_structure) );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Set the PCType.  Note: this used to be done *before* the call to
    // PCSetOperators(), and only when !M_is_initialized, but
    // 1.) Some preconditioners (those employing sub-preconditioners,
    // for example) have to call PCSetUp(), and can only do this after
    // the operators have been set.
    // 2.) It should be safe to call set_petsc_preconditioner_type()
    // multiple times.
    VLOG(2) << "prec : "  << this->M_preconditioner_type << "\n";
    setPetscPreconditionerType( this->M_preconditioner_type,this->M_matSolverPackage_type,M_pc,is,this->worldComm(),this->name() );
    VLOG(2) << "mat solver package : "  << this->M_matSolverPackage_type << "("  << Environment::vm()["pc-factor-mat-solver-package-type"].template as<std::string>() << ")\n";
    std::string type =  Environment::vm()["pc-factor-mat-solver-package-type"].template as<std::string>();
    this->setMatSolverPackageType( matSolverPackageEnumType( type ) );

    if ( option( _prefix=this->name(), _name="pc-view" ).template as<bool>() )
        check( PCView( M_pc, PETSC_VIEWER_STDOUT_WORLD ) );

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









/**
 * ConfigureKSP
 */
ConfigureKSP::ConfigureKSP( KSP& ksp,WorldComm const& worldComm, std::string const& sub,std::string const& prefix )
    :
    ConfigurePCBase( worldComm,sub,prefix ),
    M_type( option(_name="ksp-type",_sub=sub,_prefix=prefix).as<std::string>() ),
    M_rtol( option(_name="ksp-rtol",_sub=sub,_prefix=prefix).as<double>() ),
    M_maxit( option(_name="ksp-maxit",_sub=sub,_prefix=prefix).as<size_type>() ),
    M_showMonitor( option(_name="ksp-monitor",_sub=sub,_prefix=prefix).as<bool>() ),
    M_kspView( option(_name="ksp-view",_sub=sub,_prefix=prefix).as<bool>() )
{
    runConfigureKSP( ksp );
}
void
ConfigureKSP::runConfigureKSP( KSP& ksp )
{
    // set ksp type : gmres,cg,preonly,...
    this->check( KSPSetType( ksp, M_type.c_str() ) );
    // set ksp tolerance
    this->check( KSPSetTolerances( ksp,M_rtol,PETSC_DEFAULT,PETSC_DEFAULT,M_maxit ) );
    // monitor
    if ( M_showMonitor )
        this->check( KSPMonitorSet( ksp,KSPMonitorDefault,PETSC_NULL,PETSC_NULL ) );
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
    runConfigurePCML( pc,is );
}
void
ConfigurePCML::runConfigurePCML( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
    // Sets the number of levels to use with MG.
    // Must be called before any other MG routine
    this->check( PCMGSetLevels( pc, M_nLevels, PETSC_NULL) );

    if ( M_mgType=="multiplicative" ) this->check( PCMGSetType( pc, PC_MG_MULTIPLICATIVE ) );
    if ( M_mgType=="additive" ) this->check( PCMGSetType( pc, PC_MG_ADDITIVE ) );
    if ( M_mgType=="full" ) this->check( PCMGSetType( pc, PC_MG_FULL ) );
    if ( M_mgType=="kaskade" ) this->check( PCMGSetType( pc, PC_MG_KASKADE ) );

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

    // setup sub-pc
    //this->check( PCSetUp( pc ) );
    // this->check( PCView( pc, PETSC_VIEWER_STDOUT_WORLD ) );
    configurePCMLCoarse( pc,is );

    ConfigurePCMGLevels( pc, is, this->worldComm(), this->sub(), this->prefix() );
}
void
ConfigurePCML::configurePCMLCoarse( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is )
{
    //-----------------------------------------------------------//
    // get coarse-ksp
    KSP coarseksp;
    this->check( PCMGGetCoarseSolve( pc, &coarseksp) );

    // configure coarse ksp
    this->check( KSPSetFromOptions( coarseksp ) );

    ConfigureKSP kspConf( coarseksp, this->worldComm(), "", M_prefixMGCoarse );

    //this->check( KSPSetType( coarseksp, M_coarseKSPtype.c_str() ) );
    //KSPMonitorSet( coarseksp,KSPMonitorDefault,PETSC_NULL,PETSC_NULL );
    //KSPView(coarseksp,	PETSC_VIEWER_STDOUT_SELF);
    // setup coarse ksp
    //this->check( KSPSetUp( coarseksp ) );
    //-----------------------------------------------------------//
    // get coarse pc
    PC coarsepc;
    this->check( KSPGetPC( coarseksp, &coarsepc ) );
    // configure coarse pc
    this->check( PCSetFromOptions( coarsepc ) );
    this->check( PCSetType( coarsepc, M_coarsePCtype.c_str() ) );
    configurePC( coarsepc, is, this->worldComm(), "", M_prefixMGCoarse );
    //-----------------------------------------------------------//
    //this->check( KSPView(coarseksp, PETSC_VIEWER_STDOUT_SELF ) );
    this->check( PCSetUp( pc ) );
    // setup coarse ksp
    this->check( KSPSetUp( coarseksp ) );
    // setup coarse pc
    this->check( PCSetUp( coarsepc ) );
    //-----------------------------------------------------------//
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

    // Sets the number of levels to use with MG.
    // Must be called before any other MG routine
    this->check( PCMGSetLevels( pc, M_nLevels, PETSC_NULL) );

    if ( M_mgType=="multiplicative" ) this->check( PCMGSetType( pc, PC_MG_MULTIPLICATIVE ) );
    if ( M_mgType=="additive" ) this->check( PCMGSetType( pc, PC_MG_ADDITIVE ) );
    if ( M_mgType=="full" ) this->check( PCMGSetType( pc, PC_MG_FULL ) );
    if ( M_mgType=="kaskade" ) this->check( PCMGSetType( pc, PC_MG_KASKADE ) );

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
    //PreconditionerPetsc<double>::setPetscMGLevelsPreconditionerType( pc, is, this->worldComm(), this->prefix() );
    ConfigurePCMGLevels( pc, is, this->worldComm(), this->sub(), this->prefix() );
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

    //KSPView(coarseksp,	PETSC_VIEWER_STDOUT_SELF);
    // setup coarse ksp
    this->check( KSPSetUp( coarseksp ) );

    // get coarse pc
    PC coarsepc;
    this->check( KSPGetPC( coarseksp, &coarsepc ) );
    // configure coarse pc
    this->check( PCSetFromOptions( coarsepc ) );
    this->check( PCSetType( coarsepc, M_coarsePCtype.c_str() ) );
    configurePC( coarsepc, is, this->worldComm(), "", M_prefixMGCoarse );

    // setup coarse ksp
    //this->check( KSPSetUp( coarseksp ) );
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
    this->check( PCMGGetLevels( pc, &M_nLevels) );

    M_prefixMGLevels.resize( M_nLevels-1 );
    M_mgLevelsPCtype.resize( M_nLevels-1 );
    M_mgLevelsPCview.resize( M_nLevels-1 );
    M_mgLevelsKSPview.resize( M_nLevels-1 );
    std::string mgctx = (sub.empty())? "mg-" : sub+"-mg-";
    for ( int level=1; level<M_nLevels; ++level )
    {
        std::string prefixCurLevel = ( boost::format( "%1%%2%levels%3%" ) %prefixvm( this->prefix(),"" ) %mgctx %level ).str();
        M_prefixMGLevels[level-1] = prefixCurLevel;
        M_mgLevelsKSPview[level-1] = option(_name="ksp-view",_prefix=prefixCurLevel).as<bool>();
        M_mgLevelsPCtype[level-1] = option(_name="pc-type",_prefix=prefixCurLevel).as<std::string>();
        M_mgLevelsPCview[level-1] = option(_name="pc-view",_prefix=prefixCurLevel).as<bool>();
    }

    for ( int level=1; level<M_nLevels; ++level )
        runConfigurePCMGLevels( pc,is,level );
}

void
ConfigurePCMGLevels::runConfigurePCMGLevels( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is, int level )
{
    std::string prefixCurLevel = M_prefixMGLevels[level-1];

    //-------------------------------------------------------------------//
    KSP levelksp;
    // get ksp
    this->check( PCMGGetSmoother( pc, level, &levelksp ) );
    // init ksp from option
    this->check( KSPSetFromOptions( levelksp ) );

    ConfigureKSP kspConf( levelksp,this->worldComm(), "", prefixCurLevel );

    // warning : use KSP_NORM_PRECONDITIONED and force convergence
    this->check( KSPSetNormType( levelksp, KSP_NORM_PRECONDITIONED ) );
    void *cctx;
    this->check( KSPDefaultConvergedCreate(&cctx) );
    this->check( KSPSetConvergenceTest( levelksp, KSPDefaultConverged, cctx, PETSC_NULL ) );

    // setup coarse ksp
    this->check( KSPSetUp( levelksp ) );

    //-------------------------------------------------------------------//
    // get level pc
    PC levelpc;
    this->check( KSPGetPC( levelksp, &levelpc ) );
    // configure level pc
    this->check( PCSetFromOptions( levelpc ) );
    this->check( PCSetType( levelpc, M_mgLevelsPCtype[level-1].c_str() ) );
    configurePC( levelpc, is, this->worldComm(), "", prefixCurLevel );
    // setup level pc
    this->check( PCSetUp( levelpc ) );

    PetscViewer viewer = (this->sub().empty())? PETSC_VIEWER_STDOUT_WORLD : PETSC_VIEWER_STDOUT_SELF;
    if ( kspConf.kspView() )
        this->check( KSPView( levelksp, viewer ) );
    else if ( M_mgLevelsPCview[level-1] )
        this->check( PCView( levelpc, viewer ) );

}




void
configurePC( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
             WorldComm const& worldComm, std::string const& sub/* = ""*/, std::string const& prefix/* = ""*/ )
{
    VLOG(2) << "configuring PC... (sub: " << sub << ")";
    google::FlushLogFiles(google::INFO);
    const char* pctype;
    int ierr = PCGetType ( pc, &pctype );
    CHKERRABORT( worldComm.globalComm(),ierr );
    VLOG(2) << "configuring PC (" << prefix << "." << sub << ")" << pctype <<  "\n";
    google::FlushLogFiles(google::INFO);

    if ( std::string(pctype) == "gasm" || std::string(pctype) == "asm" || std::string(pctype) == "pbjacobi" )
    {
        if ( std::string(pctype) == "gasm" )
        {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
            std::string t = Environment::vm(_name="pc-gasm-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>();
            if ( t == "restrict" ) PCGASMSetType( pc, PC_GASM_RESTRICT );
            if ( t == "basic" ) PCGASMSetType( pc, PC_GASM_BASIC );
            if ( t == "interpolate" ) PCGASMSetType( pc, PC_GASM_INTERPOLATE );
            if ( t == "none" ) PCGASMSetType( pc, PC_GASM_NONE );

            int levels = Environment::vm(_name="pc-gasm-overlap",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>();
            ierr = PCGASMSetOverlap( pc, levels );
            CHKERRABORT( worldComm.globalComm(),ierr );
#endif
        }
        else if ( std::string(pctype) == "asm" )
        {
            std::string t = Environment::vm(_name="pc-asm-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>();
            if ( t == "restrict" ) PCASMSetType( pc, PC_ASM_RESTRICT );
            if ( t == "basic" ) PCASMSetType( pc, PC_ASM_BASIC );
            if ( t == "interpolate" ) PCASMSetType( pc, PC_ASM_INTERPOLATE );
            if ( t == "none" ) PCASMSetType( pc, PC_ASM_NONE );

            int levels = Environment::vm(_name="pc-asm-overlap",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>();
            ierr = PCASMSetOverlap( pc, levels );
            CHKERRABORT( worldComm.globalComm(),ierr );
        }

        PreconditionerPetsc<double>::setPetscSubpreconditionerType( pc, is, worldComm, prefix );
    }
    if ( std::string(pctype) == "lu" )
    {
        std::string t = Environment::vm(_name="pc-factor-mat-solver-package-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>();
        VLOG(2) << "mat solver package: " << t << "\n";
        google::FlushLogFiles(google::INFO);
        // set factor package
        ierr = PCFactorSetMatSolverPackage( pc, t.c_str() );
        CHKERRABORT( worldComm.globalComm(),ierr );

    }
    if ( std::string(pctype) == "ilu" )
    {
        // do we need to set the mat solver package for ilu ?
        //PetscPCFactorSetMatSolverPackage( pc, "petsc" );
        ierr = PCFactorSetLevels( pc, Environment::vm(_name="pc-factor-levels",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>() );
        CHKERRABORT( worldComm.globalComm(),ierr );
        ierr = PCFactorSetFill( pc, Environment::vm(_name="pc-factor-fill",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<double>() );
        CHKERRABORT( worldComm.globalComm(),ierr );
    }
    if ( std::string(pctype) == "gamg" )
    {
        ConfigurePCGAMG( pc, is, worldComm, sub, prefix );
    }
    if ( std::string(pctype) == "ml" )
    {
        ConfigurePCML( pc, is, worldComm, sub, prefix );
    }
    if ( std::string(pctype) == "fieldsplit" )
    {
        //std::cout << " precondi petsc\n";
        /*const PetscInt ufields[] = {0,2},pfields[] = {1};
        ierr = PCFieldSplitSetFields( pc , NULL, 2, ufields,ufields);
        CHKERRABORT( worldComm.globalComm(),ierr );
        ierr = PCFieldSplitSetFields( pc , NULL, 1, pfields,pfields);
        CHKERRABORT( worldComm.globalComm(),ierr );*/

        PCCompositeType theFieldSplitType = PC_COMPOSITE_SCHUR;
        std::string t = Environment::vm(_name="fieldsplit-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>();
        if ( t == "schur" ) theFieldSplitType = PC_COMPOSITE_SCHUR;
        if ( t == "additive" ) theFieldSplitType = PC_COMPOSITE_ADDITIVE;
        if ( t == "multiplicative" ) theFieldSplitType = PC_COMPOSITE_MULTIPLICATIVE;
        if ( t == "symmetric-multiplicative" ) theFieldSplitType = PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE;
        if ( t == "special" ) theFieldSplitType = PC_COMPOSITE_SPECIAL;
        ierr = PCFieldSplitSetType( pc, theFieldSplitType );
        CHKERRABORT( worldComm.globalComm(),ierr );

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
        if ( t == "schur" )
        {
            PCFieldSplitSchurFactType theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_FULL;
            std::string t2 = Environment::vm(_name="fieldsplit-schur-fact-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>();
            if (t2 == "diag")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_DIAG;
            if (t2 == "lower")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_LOWER;
            if (t2 == "upper")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_UPPER;
            if (t2 == "full")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_FULL;
            ierr = PCFieldSplitSetSchurFactType( pc,theSchurFactType );
            CHKERRABORT( worldComm.globalComm(),ierr );

            PCFieldSplitSchurPreType theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_SELF;
            std::string t3 = Environment::vm(_name="fieldsplit-schur-precondition",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>();
            if (t3 == "self")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_SELF;
            if (t3 == "user")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_USER;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,4,0 )
            if (t3 == "a11")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_A11;
#else
            if (t3 == "a11")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_DIAG;
#endif
            ierr = PCFieldSplitSchurPrecondition( pc, theSchurPrecond, NULL );
            CHKERRABORT( worldComm.globalComm(),ierr );
        }
#endif

        PreconditionerPetsc<double>::setPetscFieldSplitPreconditionerType( pc, is, worldComm, prefix );
    }

    VLOG(2) << "configuring PC " << pctype << " done\n";
    google::FlushLogFiles(google::INFO);
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

template <typename T>
void PreconditionerPetsc<T>::setPetscPreconditionerType ( const PreconditionerType & preconditioner_type,
                                                          const MatSolverPackageType & matSolverPackage_type,
                                                          PC & pc,
                                                          indexsplit_ptrtype const& is,
                                                          WorldComm const& worldComm,
                                                          std::string const& name )

{
    //mpi::communicator world;
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
        // In serial, just set the ILU preconditioner type
        //if (Feel::n_processors() == 1)
        // change in parallel version
        if ( worldComm.globalSize() == 1 )
        {
            ierr = PCSetType ( pc, ( char* ) PCILU );
            CHKERRABORT( worldComm.globalComm(),ierr );
        }

        else
        {
            // But PETSc has no truly parallel ILU, instead you have to set
            // an actual parallel preconditioner (e.g. block Jacobi) and then
            // assign ILU sub-preconditioners.
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
            ierr = PCSetType ( pc, ( char* ) PCGASM );
#else
            ierr = PCSetType ( pc, ( char* ) PCASM );
#endif
            CHKERRABORT( worldComm.globalComm(),ierr );

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

    case AMG_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCHYPRE );
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
        ierr = PCSetFromOptions( pc );
        CHKERRABORT( worldComm.globalComm(),ierr );

        break;

    case ML_PRECOND:
        ierr = PCSetType( pc,( char* ) PCML );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 2)
    case GAMG_PRECOND:
        ierr = PCSetType( pc,( char* ) PCGAMG );
        CHKERRABORT( worldComm.globalComm(),ierr );
#if 1
        // crash without this
        //ierr = PCSetFromOptions( pc );
        //CHKERRABORT( worldComm.globalComm(),ierr );
#endif
        break;
#else
        LOG(WARNING) << "PETSc GAMG is available from PETSc version >= 3.2";
#endif

    case NONE_PRECOND:
        ierr = PCSetType( pc,( char* ) PCNONE );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    default:
        std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
                  << preconditioner_type       << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }

    // Set additional options if we are doing AMG and
    // HYPRE is available
#ifdef FEELPP_HAS_PETSC_HYPRE

    if ( preconditioner_type == AMG_PRECOND )
    {
        ierr = PCHYPRESetType( pc, "boomeramg" );
        CHKERRABORT( worldComm.globalComm(),ierr );
    }

#endif

    // configure main preconditioner
    configurePC( pc, is, worldComm, "", name );

#if 0

    // configure sub pc/ksp
    if ( preconditioner_type == ASM_PRECOND ||
         preconditioner_type == GASM_PRECOND ||
         preconditioner_type == BLOCK_JACOBI_PRECOND )
        setPetscSubpreconditionerType( pc, worldComm, name );
    else if ( preconditioner_type == FIELDSPLIT_PRECOND )
        setPetscFieldSplitPreconditionerType( pc, is, worldComm, name );
#endif
    // prepare PC to use
    ierr = PCSetUp( pc );
    CHKERRABORT( worldComm.globalComm(),ierr );
}



template <typename T>
#if PETSC_VERSION_LESS_THAN(3,0,0)
void PreconditionerPetsc<T>::setPetscSubpreconditionerType( PC& pc, std::string const& prefix  )
#else
    void PreconditionerPetsc<T>::setPetscSubpreconditionerType( PC& pc, indexsplit_ptrtype const& is, WorldComm const& worldComm, std::string const& prefix )
#endif
{
    // For catching PETSc error return codes
    int ierr = 0;

    // All docs say must call KSPSetUp or PCSetUp before calling PCBJacobiGetSubKSP.
    // You must call PCSetUp after the preconditioner operators have been set, otherwise you get the:
    //
    // "Object is in wrong state!"
    // "Matrix must be set first."
    //
    // error messages...
    ierr = PCSetUp( pc );
    CHKERRABORT( worldComm.globalComm(),ierr );
#if PETSC_VERSION_LESS_THAN(3,4,0)
    const PCType thepctype;
#else
    PCType thepctype;
#endif
    ierr = PCGetType( pc, &thepctype );
    CHKERRABORT( worldComm.globalComm(),ierr );
    // To store array of local KSP contexts on this processor
    KSP* subksps;

    // the number of blocks on this processor
    int n_local = 1;

    // The global number of the first block on this processor.
    // This is not used, so we just pass PETSC_NULL instead.
    // int first_local;
    // Fill array of local KSP contexts
    VLOG(2) << "[setPetscSubpreconditionerType] set local preconditioner for preconditioner : " << thepctype << "\n";
    google::FlushLogFiles(google::INFO);
    if ( std::string( thepctype ) == "block_jacobi" || std::string( thepctype ) == "bjacobi" )
        ierr = PCBJacobiGetSubKSP( pc, &n_local, PETSC_NULL, &subksps );
    else if ( std::string( thepctype ) == "asm" )
        ierr = PCASMGetSubKSP( pc, &n_local, PETSC_NULL, &subksps );
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
    else if ( std::string( thepctype ) == "gasm" )
        ierr = PCGASMGetSubKSP( pc, &n_local, PETSC_NULL, &subksps );
#endif
    CHKERRABORT( worldComm.globalComm(),ierr );

    VLOG(2) << "number of sub ksp : " << n_local;
    if ( worldComm.globalSize() > 1 )
    {
        std::string subpctype = option(_name="pc-type",_sub="sub",_prefix=prefix).template as<std::string>();
        VLOG(2) << "subpctype: " << subpctype ;
        google::FlushLogFiles(google::INFO);
        // Loop over sub-ksp objects, set ILU preconditioner
        for ( int i=0; i<n_local; ++i )
        {
            // Get pointer to sub KSP object's PC
            PC subpc;

            ierr = KSPGetPC( subksps[i], &subpc );
            CHKERRABORT( worldComm.globalComm(),ierr );

            // Set requested type on the sub PC
            ierr = PCSetType( subpc, subpctype.c_str() );
            CHKERRABORT( worldComm.globalComm(),ierr );
            VLOG(2) << "sub pc " << i << "\n";
            google::FlushLogFiles(google::INFO);
            configurePC( subpc, is, worldComm, "sub", prefix );

            ierr = PCSetUp( subpc );
            CHKERRABORT( worldComm.globalComm(),ierr );

            if ( option(_name="pc-view",_sub="sub",_prefix=prefix).template as<bool>() )
            {
                ierr = PCView( subpc, PETSC_VIEWER_STDOUT_SELF );
                CHKERRABORT( worldComm.globalComm(),ierr );
            }

        }
    }
}


template <typename T>
void
PreconditionerPetsc<T>::setPetscFieldSplitPreconditionerType( PC& pc, indexsplit_ptrtype const& is,
                                                              WorldComm const& worldComm,
                                                              std::string const& prefix )
{
    // For catching PETSc error return codes
    int ierr = 0;

    // call necessary before PCFieldSplitGetSubKSP
    ierr = PCSetUp( pc );
    CHKERRABORT( worldComm.globalComm(),ierr );

    // the number of blocks on this processor
    int n_local;
    // To store array of local KSP contexts on this processor
    KSP* subksps;
    ierr = PCFieldSplitGetSubKSP(pc,&n_local,&subksps );
    CHKERRABORT( worldComm.globalComm(),ierr );

    // Loop over sub-ksp objects, set ILU preconditioner
    for ( int i=0; i<n_local; ++i )
    {

        std::string prefixSplit = prefixvm(prefix , (boost::format( "fieldsplit-%1%" )  %i ).str() );
        VLOG(2) << "configure split " << i << " with prefix "<< prefixSplit << "\n";
        google::FlushLogFiles(google::INFO);

        std::string subksptype = option(_name="ksp-type",_prefix=prefixSplit).template as<std::string>();
        VLOG(2) << " subksptype " << subksptype << "\n";

        ierr = KSPSetType ( subksps[i], subksptype.c_str() );
        CHKERRABORT( worldComm.globalComm(),ierr );
#if 0
        Mat A00;Mat A01;Mat A10; Mat A11;
        ierr = PCFieldSplitGetSchurBlocks(pc,&A00,&A01,&A10, &A11);
        CHKERRABORT( worldComm.globalComm(),ierr );
        if (i==0) {
        ierr = KSPSetOperators( subksps[i], A00, A00,
                                PetscGetMatStructureEnum(MatrixStructure::SAME_PRECONDITIONER));
        CHKERRABORT( worldComm.globalComm(),ierr ); }
#endif


        // Get pointer to sub KSP object's PC
        PC subpc;
        ierr = KSPGetPC( subksps[i], &subpc );
        CHKERRABORT( worldComm.globalComm(),ierr );

        // Set requested type on the sub PC
        std::string subpctype = option(_name="pc-type",_prefix=prefixSplit).template as<std::string>();
        ierr = PCSetType( subpc, subpctype.c_str() );
        CHKERRABORT( worldComm.globalComm(),ierr );

        // init pc from specific default option
        ierr = PCSetFromOptions( subpc );
        CHKERRABORT( worldComm.globalComm(),ierr );

        //LOG(INFO) << "configure split " << i << " (" << prefixSplit << ")" << subpctype <<  "\n";
        //google::FlushLogFiles(google::INFO);

        if ( subpctype == "fieldsplit" )
        {
            CHECK( is ) << "index split is not initialized\n";

            std::string fieldsDefStr = option(_name="fieldsplit-fields",_prefix=prefixSplit).template as<std::string>();
            auto fieldsDef = IndexSplit::parseFieldsDef( fieldsDefStr );
            //std::cout << "fieldsplit fieldsDefStr " << fieldsDefStr << "\n";
            //auto isUsed = is.applyFieldsDef( IndexSplit::FieldsDef(  {  { 0 , { 0 } }, { 1 , { 2 } } } ) );
            //is->showMe();
            auto isUsed = is->applyFieldsDef( IndexSplit::FieldsDef( fieldsDef ) );
            //isUsed->showMe();

            std::vector<IS> isPetsc;
            PetscConvertIndexSplit( isPetsc ,*isUsed,worldComm);
            for ( int i = 0 ; i < isPetsc.size(); ++i )
            {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
                std::ostringstream os; os << i;
                ierr=PCFieldSplitSetIS( subpc,os.str().c_str(),isPetsc[i] );
#else
                ierr=PCFieldSplitSetIS( subpc,isPetsc[i] );
#endif
                CHKERRABORT( worldComm,ierr );
            }

            //TODO : maybe delete or store  isPetsc

        }

        // configure sub PC
        //std::string prefixPetsc=(boost::format("fieldsplit_%1%_")%i ).str();
        configurePC( subpc, is, worldComm, "", prefixSplit/*, prefixPetsc*/ );

        // configure maybe sub sub PC
        const char* thesubpctype;
        ierr = PCGetType( subpc, &thesubpctype );
        CHKERRABORT( worldComm.globalComm(),ierr );
        if ( std::string( thesubpctype ) == "block_jacobi" || std::string( thesubpctype ) == "bjacobi" ||
             std::string( thesubpctype ) == "asm" || std::string( thesubpctype ) == "gasm" )
        {
            //setPetscSubpreconditionerType( subpc, worldComm, prefixSplit );
        }
        else if ( std::string(thesubpctype) == "lsc" )
        {
            CHECK( i==1 ) << "lsc must be use with only field 1, not " << i << "\n";

#if 0
            std::string prefixFeelBase = prefixvm(prefixSplit,"lsc");
            std::string prefixPetscBase = "fieldsplit_1_lsc";
            //std::cout << "USE LSC with " << prefixFeelBase << " et " << prefixPetscBase << "\n";
            configurePCWithPetscCommandLineOption( prefixFeelBase,prefixPetscBase );
#endif

            setPetscLSCPreconditionerType( subpc, is, worldComm, prefixSplit );

        }
        else if ( std::string(thesubpctype) == "gamg" )
        {
            //ierr = PCSetFromOptions( subpc );
            //CHKERRABORT( worldComm.globalComm(),ierr );
        }
        else if ( std::string(thesubpctype) == "fieldsplit" )
        {
            //setPetscFieldSplitPreconditionerType( subpc, is, worldComm, prefixvm( prefix, prefixSplit ) );
        }

        // setup sub-pc
        ierr = PCSetUp( subpc );
        CHKERRABORT( worldComm.globalComm(),ierr );

        if ( option(_name="pc-view",_prefix=prefixSplit).template as<bool>() )
        {
            ierr = PCView( subpc, PETSC_VIEWER_STDOUT_WORLD );
            CHKERRABORT( worldComm.globalComm(),ierr );
        }

    }

}


template <typename T>
void
PreconditionerPetsc<T>::setPetscLSCPreconditionerType( PC& pc,indexsplit_ptrtype const& is,
                                                       WorldComm const& worldComm,
                                                       std::string const& prefix )
{
    std::string prefixLSC = prefixvm(prefix,"lsc");
    int ierr = 0;
    //-----------------------------------------------------------//
    // setup sub-pc
    ierr = PCSetUp( pc );
    CHKERRABORT( worldComm.globalComm(),ierr );
    //-----------------------------------------------------------//
    // get sub-ksp
    KSP subksp;
    PCLSCGetKSP( pc, subksp );
    // configure sub-ksp
    ierr = KSPSetFromOptions( subksp );
    CHKERRABORT( worldComm.globalComm(),ierr );
    std::string subksptype = option(_name="ksp-type",_prefix=prefixLSC).template as<std::string>();
    ierr = KSPSetType( subksp, subksptype.c_str() );
    CHKERRABORT( worldComm.globalComm(),ierr );
    // setup sub-ksp
    ierr = KSPSetUp( subksp );
    CHKERRABORT( worldComm.globalComm(),ierr );
    //-----------------------------------------------------------//
    // get sub-pc
    PC subpc;
    ierr = KSPGetPC( subksp, &subpc );
    CHKERRABORT( worldComm.globalComm(),ierr );
    // configure sub-pc
    ierr = PCSetFromOptions( subpc );
    CHKERRABORT( worldComm.globalComm(),ierr );
    std::string subPCtype =  option(_name="pc-type",_prefix=prefixLSC).template as<std::string>();
    ierr = PCSetType( subpc, subPCtype.c_str() );
    CHKERRABORT( worldComm.globalComm(),ierr );
    configurePC( subpc, is, worldComm, "", prefixLSC );
    // setup sub-pc
    ierr = PCSetUp( subpc );
    CHKERRABORT( worldComm.globalComm(),ierr );
    //-----------------------------------------------------------//

    if ( option(_name="pc-view",_prefix=prefixLSC).template as<bool>() )
    {
        ierr = PCView( subpc, PETSC_VIEWER_STDOUT_WORLD );
        CHKERRABORT( worldComm.globalComm(),ierr );
    }

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
