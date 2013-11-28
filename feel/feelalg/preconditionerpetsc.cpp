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
            pmatrix->updatePCFieldSplit( M_pc );
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
    setPetscPreconditionerType( this->M_preconditioner_type,this->M_matSolverPackage_type,M_pc,this->worldComm(),this->name() );
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


void
configurePC( PC& pc, WorldComm const& worldComm, std::string sub = "", std::string prefix = "" )
{
    VLOG(2) << "configuring PC... (sub: " << sub << ")";
    google::FlushLogFiles(google::INFO);
    const char* pctype;
    int ierr = PCGetType ( pc, &pctype );
    CHKERRABORT( worldComm.globalComm(),ierr );
    VLOG(2) << "configuring PC (" << prefix << "." << sub << ")" << pctype <<  "\n";
    google::FlushLogFiles(google::INFO);
    if ( std::string(pctype) == "gasm" )
    {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
        std::string t = Environment::vm(_name="pc-gasm-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>();
        if ( t == "restrict" ) PCGASMSetType( pc, PC_GASM_RESTRICT );
        if ( t == "basic" ) PCGASMSetType( pc, PC_GASM_BASIC );
        if ( t == "interpolate" ) PCGASMSetType( pc, PC_GASM_INTERPOLATE );
        if ( t == "none" ) PCGASMSetType( pc, PC_GASM_NONE );

        int levels = Environment::vm(_name="pc-gasm-overlap",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>();
        PCGASMSetOverlap( pc, levels );
#endif
    }
    if ( std::string(pctype) == "asm" )
    {
        std::string t = Environment::vm(_name="pc-asm-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>();
        if ( t == "restrict" ) PCASMSetType( pc, PC_ASM_RESTRICT );
        if ( t == "basic" ) PCASMSetType( pc, PC_ASM_BASIC );
        if ( t == "interpolate" ) PCASMSetType( pc, PC_ASM_INTERPOLATE );
        if ( t == "none" ) PCASMSetType( pc, PC_ASM_NONE );

        int levels = Environment::vm(_name="pc-asm-overlap",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>();
        PCASMSetOverlap( pc, levels );

    }
    if ( std::string(pctype) == "lu" )
    {
        std::string t = Environment::vm(_name="pc-factor-mat-solver-package-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>();
        VLOG(2) << "mat solver package: " << t << "\n";
        google::FlushLogFiles(google::INFO);
        // set factor package
        PCFactorSetMatSolverPackage( pc, t.c_str() );

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
    if ( std::string(pctype) == "ml" || std::string(pctype) == "mg")
    {
        int nLevels= option(_name="pc-mg-levels",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>();
        //std::vector<MPI_Comm> comms(levels,worldComm.globalComm());
        ierr = PCMGSetLevels( pc, nLevels, /*comms.data()*/ PETSC_NULL);
        CHKERRABORT( worldComm.globalComm(),ierr );

        std::string mgType = option(_name="pc-mg-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>();
        if ( mgType=="multiplicative" ) ierr = PCMGSetType( pc, PC_MG_MULTIPLICATIVE );
        if ( mgType=="additive" ) ierr = PCMGSetType( pc, PC_MG_ADDITIVE );
        if ( mgType=="full" ) ierr = PCMGSetType( pc, PC_MG_FULL );
        if ( mgType=="kaskade" ) ierr = PCMGSetType( pc, PC_MG_KASKADE );
        CHKERRABORT( worldComm.globalComm(),ierr );

        int smoothdown= option(_name="pc-mg-smoothdown",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>();
        ierr = PCMGSetNumberSmoothDown( pc, smoothdown );
        CHKERRABORT( worldComm.globalComm(),ierr );

        if ( std::string(pctype) == "ml" )
        {
            std::string option_pc_ml_maxNlevels = ( boost::format("-pc_ml_maxNlevels %1%") %nLevels ).str();
            ierr = PetscOptionsClearValue( "-pc_ml_maxNlevels" );
            ierr = PetscOptionsInsertString( (option_pc_ml_maxNlevels).c_str() );

            std::string option_pc_ml_reuse_interpolation = "-pc_ml_reuse_interpolation";
            bool mlReuseInterp = option(_name="pc-ml-reuse-interpolation",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<bool>();
            std::string mlReuseInterpStr = (mlReuseInterp)?"true":"false";
            ierr = PetscOptionsClearValue( option_pc_ml_reuse_interpolation.c_str() );
            ierr = PetscOptionsInsertString( (option_pc_ml_reuse_interpolation+" "+mlReuseInterpStr).c_str() );

            std::string option_pc_ml_KeepAggInfo = "-pc_ml_KeepAggInfo";
            bool mlKeepAggInfo = option(_name="pc-ml-keep-agg-info",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<bool>();
            std::string mlKeepAggInfoStr = (mlKeepAggInfo)?"true":"false";
            ierr = PetscOptionsClearValue( option_pc_ml_KeepAggInfo.c_str() );
            ierr = PetscOptionsInsertString( (option_pc_ml_KeepAggInfo+" "+mlKeepAggInfoStr).c_str() );

            std::string option_pc_ml_Reusable = "-pc_ml_Reusable";
            bool mlReusable = option(_name="pc-ml-reusable",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<bool>();
            std::string mlReusableStr = (mlReusable)?"true":"false";
            ierr = PetscOptionsClearValue( option_pc_ml_Reusable.c_str() );
            ierr = PetscOptionsInsertString( (option_pc_ml_Reusable+" "+mlReusableStr).c_str() );

            std::string option_pc_ml_OldHierarchy = "-pc_ml_OldHierarchy";
            bool mlOldHierarchy = option(_name="pc-ml-old-hierarchy",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<bool>();
            std::string mlOldHierarchyStr = (mlOldHierarchy)?"true":"false";
            ierr = PetscOptionsClearValue( option_pc_ml_OldHierarchy.c_str() );
            ierr = PetscOptionsInsertString( (option_pc_ml_OldHierarchy+" "+mlOldHierarchyStr).c_str() );
        }

        std::string name = (!sub.empty())?sub+"-"+prefix:prefix;
        PreconditionerPetsc<double>::setPetscMGCoarsePreconditionerType( pc, worldComm, name );
        PreconditionerPetsc<double>::setPetscMGLevelsPreconditionerType( pc, worldComm, name );

        ierr = PCSetFromOptions( pc );
        CHKERRABORT( worldComm.globalComm(),ierr );

    }


    if ( Environment::vm(_name="pc-view",_sub=sub,_prefix=prefix).as<bool>() )
        PCView( pc, PETSC_VIEWER_STDOUT_SELF );

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
        if ( worldComm.globalSize() == 1 || matSolverPackage_type == MATSOLVER_MUMPS)
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
        break;

    case ML_PRECOND:
        ierr = PCSetType( pc,( char* ) PCML );
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
    // before configurePC else pc-view doesn't work really
    if ( preconditioner_type == FIELDSPLIT_PRECOND )
        setPetscFieldSplitPreconditionerType( pc, worldComm, name );

    configurePC( pc, worldComm, "", name );

    if ( preconditioner_type == ASM_PRECOND ||
         preconditioner_type == GASM_PRECOND ||
         preconditioner_type == BLOCK_JACOBI_PRECOND )
        setPetscSubpreconditionerType( pc, worldComm, name );

}



template <typename T>
#if PETSC_VERSION_LESS_THAN(3,0,0)
void PreconditionerPetsc<T>::setPetscSubpreconditionerType( PC& pc, std::string const& prefix  )
#else
    void PreconditionerPetsc<T>::setPetscSubpreconditionerType( PC& pc, WorldComm const& worldComm, std::string const& prefix )
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
    if ( Environment::numberOfProcessors() > 1 )
    {
        std::string subpctype =  Environment::vm(_name="pc-type",_sub="sub",_prefix=prefix).template as<std::string>();
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
            configurePC( subpc, worldComm, "sub", prefix );
        }
    }
}


template <typename T>
void
PreconditionerPetsc<T>::setPetscFieldSplitPreconditionerType( PC& pc,
                                                              WorldComm const& worldComm,
                                                              std::string const& prefix )
{
    // For catching PETSc error return codes
    int ierr = 0;

    PCCompositeType theFieldSplitType = PC_COMPOSITE_SCHUR;
    std::string t = Environment::vm(_name="fieldsplit-type",_prefix=prefix,_worldcomm=worldComm).template as<std::string>();
    if ( t == "schur" ) theFieldSplitType = PC_COMPOSITE_SCHUR;
    if ( t == "additive" ) theFieldSplitType = PC_COMPOSITE_ADDITIVE;
    if ( t == "multiplicative" ) theFieldSplitType = PC_COMPOSITE_MULTIPLICATIVE;
    if ( t == "symmetric-multiplicative" ) theFieldSplitType = PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE;
    if ( t == "special" ) theFieldSplitType = PC_COMPOSITE_SPECIAL;

    ierr = PCFieldSplitSetType( pc, theFieldSplitType );
    CHKERRABORT( worldComm.globalComm(),ierr );

    if ( t == "schur" )
    {
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
        PCFieldSplitSchurFactType theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_FULL;
        std::string t2 = Environment::vm(_name="fieldsplit-schur-fact-type",_prefix=prefix,_worldcomm=worldComm).template as<std::string>();
        if (t2 == "diag")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_DIAG;
        if (t2 == "lower")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_LOWER;
        if (t2 == "upper")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_UPPER;
        if (t2 == "full")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_FULL;

        ierr = PCFieldSplitSetSchurFactType( pc,theSchurFactType );
        CHKERRABORT( worldComm.globalComm(),ierr );
#endif
    }


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

        std::string subksptype =  Environment::vm(_name="ksp-type",_prefix=prefixSplit).template as<std::string>();
        //std::cout<< " subksptype " << subksptype << std::endl;
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

        VLOG(2) << "configure split " << i << "\n";
        google::FlushLogFiles(google::INFO);

        // Get pointer to sub KSP object's PC
        PC subpc;
        ierr = KSPGetPC( subksps[i], &subpc );
        CHKERRABORT( worldComm.globalComm(),ierr );

        // Set requested type on the sub PC
        std::string subpctype =  Environment::vm(_name="pc-type",_prefix=prefixSplit).template as<std::string>();
        ierr = PCSetType( subpc, subpctype.c_str() );
        CHKERRABORT( worldComm.globalComm(),ierr );

        //LOG(INFO) << "configure split " << i << " (" << prefixSplit << ")" << subpctype <<  "\n";
        //google::FlushLogFiles(google::INFO);

        // configure sub PC
        configurePC( subpc, worldComm, "", prefixSplit );

        // configure maybe sub sub PC
        const char* thesubpctype;
        ierr = PCGetType( subpc, &thesubpctype );
        CHKERRABORT( worldComm.globalComm(),ierr );
        if ( std::string( thesubpctype ) == "block_jacobi" || std::string( thesubpctype ) == "bjacobi" ||
             std::string( thesubpctype ) == "asm" ||
             std::string( thesubpctype ) == "gasm" )
        setPetscSubpreconditionerType( subpc, worldComm, prefixSplit );

    }

}


template <typename T>
void
PreconditionerPetsc<T>::setPetscMGCoarsePreconditionerType( PC& pc,
                                                            WorldComm const& worldComm,
                                                            std::string const& prefix )
{
    // For catching PETSc error return codes
    int ierr = 0;

    KSP coarseksp;
    ierr = PCMGGetCoarseSolve( pc, &coarseksp);
    CHKERRABORT( worldComm.globalComm(),ierr );

    //KSPMonitorSet( coarseksp,KSPMonitorDefault,PETSC_NULL,PETSC_NULL );
    //KSPView(coarseksp,	PETSC_VIEWER_STDOUT_SELF);

    PC coarsepc;
    ierr = KSPGetPC( coarseksp, &coarsepc );
    CHKERRABORT( worldComm.globalComm(),ierr );

    std::string prefixMGCoarse = ( boost::format( "%1%mg-coarse" ) %prefixvm( prefix,"" ) ).str();
    std::string mgCoarsePCtype =  Environment::vm(_name="pc-type",_prefix=prefixMGCoarse).template as<std::string>();
    ierr = PCSetType( coarsepc, mgCoarsePCtype.c_str() );
    CHKERRABORT( worldComm.globalComm(),ierr );

    configurePC( coarsepc, worldComm, "", prefixMGCoarse );
    //ierr = PCSetUp( coarsepc );
}

template <typename T>
void
PreconditionerPetsc<T>::setPetscMGLevelsPreconditionerType( PC& pc,
                                                            WorldComm const& worldComm,
                                                            std::string const& prefix )
{
    // For catching PETSc error return codes
    int ierr = 0;

    // get number of levels
    int nLevels = 1;
    ierr = PCMGGetLevels( pc, &nLevels);
    CHKERRABORT( worldComm.globalComm(),ierr );

    // Loop over ksp objects
    for ( int k=1; k<nLevels; ++k )
    {
        std::string prefixMGLevels = ( boost::format( "%1%mg-levels%2%" ) %prefixvm( prefix,"" ) %k ).str();
        // get ksp
        KSP levelksp;
        ierr = PCMGGetSmoother( pc, k, &levelksp );
        CHKERRABORT( worldComm.globalComm(),ierr );

        // warning : use KSP_NORM_PRECONDITIONED and force convergence
        ierr = KSPSetNormType( levelksp, KSP_NORM_PRECONDITIONED );
        CHKERRABORT( worldComm.globalComm(),ierr );
        void *cctx;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,4,4)
        ierr = KSPConvergedDefaultCreate(&cctx);
        ierr = KSPSetConvergenceTest( levelksp, KSPConvergedDefault, cctx, PETSC_NULL );
#else
        ierr = KSPDefaultConvergedCreate(&cctx);
        ierr = KSPSetConvergenceTest( levelksp, KSPDefaultConverged, cctx, PETSC_NULL );
#endif
        CHKERRABORT( worldComm.globalComm(),ierr );

        // set ksp tolerance
        double levelksprtol = option(_name="ksp-rtol",_prefix=prefixMGLevels).template as<double>();
        size_type levelkspmaxit = option(_name="ksp-maxit",_prefix=prefixMGLevels).template as<size_type>();
        ierr = KSPSetTolerances ( levelksp,levelksprtol,PETSC_DEFAULT,PETSC_DEFAULT,levelkspmaxit );
        CHKERRABORT( worldComm.globalComm(),ierr );

        if (  option(_name="ksp-monitor",_prefix=prefixMGLevels).template as<bool>() )
        {
            ierr = KSPMonitorSet( levelksp,KSPMonitorDefault,PETSC_NULL,PETSC_NULL );
            CHKERRABORT( worldComm.globalComm(),ierr );
        }

        //KSPView(levelksp, PETSC_VIEWER_STDOUT_SELF);

#if 0
        //ierr = KSPSetType(levelksp,KSPRICHARDSON);
        ierr = KSPSetType(levelksp,KSPGMRES);
        CHKERRABORT( worldComm.globalComm(),ierr );

        PC levelpc;
        ierr = KSPGetPC( levelksp, &levelpc );
        CHKERRABORT( worldComm.globalComm(),ierr );
        //std::string prefixMGLevels = ( boost::format( "%1%mg-levels%2%" ) %prefixvm( prefix,"" ) %k ).str();
        std::string mgLevelsPCtype =  Environment::vm(_name="pc-type",_prefix=prefixMGLevels).template as<std::string>();

        ierr = PCSetType( levelpc, mgLevelsPCtype.c_str() );
        CHKERRABORT( worldComm.globalComm(),ierr );

        configurePC( levelpc, worldComm, "", prefixMGLevels );
        //ierr = PCSetUp( levelpc );
        //ierr = PCSetFromOptions( levelpc );
#else

        std::string mgLevelsKSPtype =  Environment::vm(_name="ksp-type",_prefix=prefixMGLevels).template as<std::string>();
        std::string mgLevelsPCtype =  Environment::vm(_name="pc-type",_prefix=prefixMGLevels).template as<std::string>();

        std::string option_ksp_type = ( boost::format( "-mg_levels_%1%_ksp_type") %k ).str();
        std::string option_pc_type = ( boost::format( "-mg_levels_%1%_pc_type") %k ).str();

        ierr = PetscOptionsClearValue( option_ksp_type.c_str() );
        ierr = PetscOptionsInsertString( (option_ksp_type+" "+mgLevelsKSPtype).c_str() );

        ierr = PetscOptionsClearValue( option_pc_type.c_str() );
        ierr = PetscOptionsInsertString( (option_pc_type+" "+mgLevelsPCtype).c_str() );

        if (mgLevelsPCtype=="lu")
        {
            std::string mgLevelsPCFMSPtype =  Environment::vm(_name="pc-factor-mat-solver-package-type",_prefix=prefixMGLevels).template as<std::string>();
            std::string option_pc_factor_mat_solver_package = ( boost::format( "-mg_levels_%1%_pc_factor_mat_solver_package") %k ).str();
            ierr = PetscOptionsClearValue( option_pc_factor_mat_solver_package.c_str() );
            ierr = PetscOptionsInsertString( (option_pc_factor_mat_solver_package+" "+mgLevelsPCFMSPtype).c_str() );
        }
        if (mgLevelsPCtype=="gasm")
        {
            int gasmoverlap = Environment::vm(_name="pc-gasm-overlap",_prefix=prefixMGLevels).template as<int>();

            std::string option_pc_gasm_overlap = ( boost::format( "-mg_levels_%1%_pc_gasm_overlap") %k ).str();
            std::string optionval_pc_gasm_overlap = ( boost::format( "-mg_levels_%1%_pc_gasm_overlap %2%") %k %gasmoverlap ).str();
            ierr = PetscOptionsClearValue( option_pc_gasm_overlap.c_str() );
            ierr = PetscOptionsInsertString( optionval_pc_gasm_overlap.c_str() );

            std::string option_sub_pc_type = ( boost::format( "-mg_levels_%1%_sub_pc_type") %k ).str();
            std::string subpctype =  Environment::vm(_name="pc-type",_sub="sub",_prefix=prefixMGLevels).template as<std::string>();
            ierr = PetscOptionsClearValue( option_sub_pc_type.c_str() );
            ierr = PetscOptionsInsertString( (option_sub_pc_type+" "+subpctype).c_str() );

            if (subpctype=="lu")
            {
                std::string option_sub_pc_factor_mat_solver_package = ( boost::format( "-mg_levels_%1%_sub_pc_factor_mat_solver_package") %k ).str();
                std::string t =  Environment::vm(_name="pc-factor-mat-solver-package-type",_sub="sub",_prefix=prefixMGLevels).template as<std::string>();
                ierr = PetscOptionsClearValue( option_sub_pc_factor_mat_solver_package.c_str() );
                ierr = PetscOptionsInsertString( (option_sub_pc_factor_mat_solver_package+" "+t).c_str() );
            }

        }

#endif
    }
}



template class PreconditionerPetsc<double>;

}
