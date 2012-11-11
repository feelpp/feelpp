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
        const PCType pc_type;
        ierr = PCGetType ( M_pc, &pc_type );
        CHKERRABORT( this->worldComm().globalComm(),ierr );

        LOG(INFO) << "preconditionerpetsc set as "  << pc_type << "\n";
        MatrixPetsc<T> * pmatrix = dynamic_cast<MatrixPetsc<T>*>( this->M_matrix.get() );

        M_mat = pmatrix->mat();
    }

    int ierr;
    ierr = PCSetOperators( M_pc,M_mat,M_mat,( MatStructure )SAME_NONZERO_PATTERN );
    CHKERRABORT( this->worldComm().globalComm(),ierr );

    // Set the PCType.  Note: this used to be done *before* the call to
    // PCSetOperators(), and only when !M_is_initialized, but
    // 1.) Some preconditioners (those employing sub-preconditioners,
    // for example) have to call PCSetUp(), and can only do this after
    // the operators have been set.
    // 2.) It should be safe to call set_petsc_preconditioner_type()
    // multiple times.
    VLOG(2) << "prec : "  << this->M_preconditioner_type << "\n";
    setPetscPreconditionerType( this->M_preconditioner_type,this->M_matSolverPackage_type,M_pc,this->worldComm() );
    VLOG(2) << "mat solver package : "  << this->M_matSolverPackage_type << "("  << Environment::vm()["pc-factor-mat-solver-package-type"].template as<std::string>() << ")\n";
    std::string type =  Environment::vm()["pc-factor-mat-solver-package-type"].template as<std::string>();
    this->setMatSolverPackageType( matSolverPackageEnumType( type ) );







    if ( Environment::vm().count( "pc-view" ) )
        PCView( M_pc, PETSC_VIEWER_STDOUT_SELF );
    this->M_is_initialized = true;
}


template <typename T>
void PreconditionerPetsc<T>::clear ()
{
    if ( this-> M_is_initialized )
    {
        this->M_is_initialized = false;

        PetscTruth is_petsc_initialized;
        PetscInitialized( &is_petsc_initialized );
        if ( is_petsc_initialized )
            PETSc::PCDestroy( M_pc );
    }

}


void
configurePC( PC& pc, WorldComm const& worldComm, std::string sub = "", std::string prefix = "" )
{
    LOG(INFO) << "configuring PC...\n";
    google::FlushLogFiles(google::INFO);
    const char* pctype;
    int ierr = PCGetType ( pc, &pctype );
    CHKERRABORT( worldComm.globalComm(),ierr );
    LOG(INFO) << "configuring PC " << pctype << "\n";
    google::FlushLogFiles(google::INFO);
    if ( std::string(pctype) == "gasm" )
    {
        std::string t = Environment::vm(_name="pc-gasm-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<std::string>();
        if ( t == "restrict" ) PCGASMSetType( pc, PC_GASM_RESTRICT );
        if ( t == "basic" ) PCGASMSetType( pc, PC_GASM_BASIC );
        if ( t == "interpolate" ) PCGASMSetType( pc, PC_GASM_INTERPOLATE );
        if ( t == "none" ) PCGASMSetType( pc, PC_GASM_NONE );

        int levels = Environment::vm(_name="pc-gasm-overlap",_prefix=prefix,_sub=sub,_worldcomm=worldComm).as<int>();
        PCGASMSetOverlap( pc, levels );

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
        LOG(INFO) << "mat solver package: " << t << "\n";
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
    LOG(INFO) << "configuring PC " << pctype << " done\n";
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
            ierr = PCSetType ( pc, ( char* ) PCGASM );
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

    case FIELDSPLIT_PRECOND: {
        ierr = PCSetType( pc,( char* ) PCFIELDSPLIT );
        CHKERRABORT( worldComm.globalComm(),ierr );
        const PCType subpctypes[5] = { PCLU,PCNONE,PCLU,PCLU,PCLU };
        const KSPType subksptypes[5] = { KSPPREONLY,KSPMINRES,KSPMINRES,KSPMINRES,KSPMINRES };
        setPetscFieldSplitPreconditionerType( PC_COMPOSITE_SCHUR, subksptypes, subpctypes, pc, worldComm );
        break; }

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

    configurePC( pc, worldComm, "", name );

    if ( preconditioner_type == ASM_PRECOND ||
         preconditioner_type == GASM_PRECOND ||
         preconditioner_type == BLOCK_JACOBI_PRECOND )
        setPetscSubpreconditionerType( pc, worldComm );
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
    const PCType thepctype;
    ierr = PCGetType( pc, &thepctype );
    CHKERRABORT( worldComm.globalComm(),ierr );
    // To store array of local KSP contexts on this processor
    KSP* subksps;

    // the number of blocks on this processor
    int n_local;

    // The global number of the first block on this processor.
    // This is not used, so we just pass PETSC_NULL instead.
    // int first_local;
    // Fill array of local KSP contexts
    LOG(INFO) << "[setPetscSubpreconditionerType] preconditioner type: " << thepctype << "\n";
    google::FlushLogFiles(google::INFO);
    if ( std::string( thepctype ) == "block_jacobi" )
        ierr = PCBJacobiGetSubKSP( pc, &n_local, PETSC_NULL, &subksps );
    else if ( std::string( thepctype ) == "asm" )
        ierr = PCASMGetSubKSP( pc, &n_local, PETSC_NULL, &subksps );
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
    else if ( std::string( thepctype ) == "gasm" )
        ierr = PCGASMGetSubKSP( pc, &n_local, PETSC_NULL, &subksps );
#endif

    CHKERRABORT( worldComm.globalComm(),ierr );
    std::string subpctype =  Environment::vm(_name="pc-type",_sub="sub",_prefix=prefix).template as<std::string>();
    LOG(INFO) << "subpctype: " << subpctype << "\n";
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
        LOG(INFO) << "pc " << i << "\n";
        google::FlushLogFiles(google::INFO);
        configurePC( subpc, worldComm, "sub", prefix );
    }

}


template <typename T>
void
PreconditionerPetsc<T>::setPetscFieldSplitPreconditionerType( const PCCompositeType type,
                                                              const KSPType * subksptypes,
                                                              const PCType * subpctypes,
                                                              PC& pc,
                                                              WorldComm const& worldComm )
{
    // For catching PETSc error return codes
    int ierr = 0;

    ierr = PCFieldSplitSetType( pc, type );
    CHKERRABORT( worldComm.globalComm(),ierr );

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
        // Set requested type on the sub KSP
        ierr = KSPSetType ( subksps[i], subksptypes[i] );
        CHKERRABORT( worldComm.globalComm(),ierr );

        // Get pointer to sub KSP object's PC
        PC subpc;
        ierr = KSPGetPC( subksps[i], &subpc );
        CHKERRABORT( worldComm.globalComm(),ierr );

        // Set requested type on the sub PC
        ierr = PCSetType( subpc, subpctypes[i] );
        CHKERRABORT( worldComm.globalComm(),ierr );
    }

}




template class PreconditionerPetsc<double>;

}
