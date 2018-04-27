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
#include <feel/feelalg/solverlinearpetsc.hpp>
#include <feel/feelpde/operatorpcdbase.hpp>


extern "C" {

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 ) && PETSC_VERSION_LESS_THAN( 3,6,0 )
#include <petsc-private/pcimpl.h>
#include <petsc-private/kspimpl.h>
#elif PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
#include <petsc/private/pcimpl.h>
#include <petsc/private/kspimpl.h>
#else
#include <private/pcimpl.h>
#include <private/kspimpl.h>
#endif

#include <feel/feelalg/preconditionerpetscpatch.cpp>
#include <feel/feelalg/preconditionerpetsclsc.cpp>
#include <feel/feelalg/preconditionerpetscpmm.cpp>
#include <feel/feelalg/preconditionerpetscpcd.cpp>
#include <feel/feelalg/preconditionerpetscfeelpp.cpp>


PetscErrorCode __feel_petsc_prec_ksp_monitor(KSP ksp,PetscInt it,PetscReal rnorm,void* ctx)
{
#if 0
    // need to store the tree of prec/subsolver else object was already deleted when go here
    Feel::ConfigureKSP *solver  = static_cast<Feel::ConfigureKSP*>( ctx );
    if ( !solver ) return 0;
    if ( solver->worldComm().isMasterRank() )
        std::cout << "   " << it << " " << solver->prefix() << " KSP Residual norm " << std::scientific << rnorm << "\n";
#elif PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,7,0)
    // a temporary solution (not quite good because because create/distroy viewer)
    MPI_Comm comm = PETSC_COMM_WORLD;
    PetscObjectGetComm((PetscObject)ksp,&comm);
    PetscViewerAndFormat *vf;
    PetscViewerAndFormatCreate( (comm == PETSC_COMM_SELF)? PETSC_VIEWER_STDOUT_SELF : PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);
    int ierr = KSPMonitorDefault( ksp,it,rnorm, vf );
    PetscViewerAndFormatDestroy( &vf );
#endif
    return 0;
}

} // extern C


namespace Feel
{
template <typename T>
void PreconditionerPetsc<T>::apply( const Vector<T> & x, Vector<T> & y ) const
{
    //if ( !this->M_is_initialized ) this->init();


    VectorPetsc<T> & x_pvec = dynamic_cast<VectorPetsc<T>&>( const_cast<Vector<T>&>( x ) );
    VectorPetsc<T> & y_pvec = dynamic_cast<VectorPetsc<T>&>( const_cast<Vector<T>&>( y ) );

    Vec x_vec = x_pvec.vec();
    Vec y_vec = y_pvec.vec();

    int ierr = PCApply( M_pc,x_vec,y_vec );
    CHKERRABORT( this->worldComm().globalComm(),ierr );
}
template <typename T>
void PreconditionerPetsc<T>::apply( Vec x, Vec y ) const
{
    int ierr = PCApply( M_pc,x,y );
    CHKERRABORT( this->worldComm().globalComm(),ierr );
}



/*----------------------- inline functions ----------------------------------*/
template <typename T>
PreconditionerPetsc<T>::PreconditionerPetsc ( std::string const& name, WorldComm const& worldComm )
    :
    Preconditioner<T>( name, worldComm ),
    M_indexSplitHasChanged( false )
{
}
template <typename T>
PreconditionerPetsc<T>::PreconditionerPetsc ( PreconditionerPetsc const& p )
    :
    Preconditioner<T>( p ),
    M_indexSplitHasChanged( p.M_indexSplitHasChanged )
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

    static bool petscPCInHouseIsInit =false;
    if ( !petscPCInHouseIsInit )
    {
        check( PCRegister("lsc2",PCCreate_LSC2) );
        check( PCRegister("pmm",PCCreate_PMM_Feelpp) );
        check( PCRegister("pcd",PCCreate_PCD_Feelpp) );
        check( PCRegister("blockns",PCCreate_FEELPP) );
        check( PCRegister("blockms",PCCreate_FEELPP) );
        petscPCInHouseIsInit=true;
    }

    M_indexSplitHasChanged = false;
    std::list<MatNullSpace> nullspList;

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

            bool useComponentsSplit = boption( _prefix=this->name(), _name="fieldsplit-use-components" );
            std::string fieldsDefStr = soption( _prefix=this->name(), _name="fieldsplit-fields" );
            auto fieldsDef = IndexSplit::parseFieldsDef( fieldsDefStr );
            if ( fieldsDef.size() == 0 )
            {
                if ( useComponentsSplit )
                {
                    auto const& isComponents = pmatrix->mapRow().indexSplitWithComponents();
                    pmatrix->updatePCFieldSplit( M_pc, isComponents );
                    for ( int splitId=0;splitId<isComponents->size();++splitId )
                        fieldsDef[splitId].insert(splitId);
                }
                else
                {
                    pmatrix->updatePCFieldSplit( M_pc );
                    for ( int splitId=0;splitId<pmatrix->indexSplit()->size();++splitId )
                        fieldsDef[splitId].insert(splitId);
                }
            }
            else
            {
                //fieldsDef.showMe();
                indexsplit_ptrtype isUsed;
                if ( useComponentsSplit )
                    isUsed = pmatrix->mapRow().indexSplitWithComponents()->applyFieldsDef( fieldsDef  );
                else
                    isUsed = pmatrix->indexSplit()->applyFieldsDef( fieldsDef  );
                //isUsed->showMe();
                pmatrix->updatePCFieldSplit( M_pc,isUsed );
            }
            //fieldsDef.showMe();
            for ( auto const& splitBaseIdsMap : fieldsDef )
            {
                int splitId = splitBaseIdsMap.first;
                std::set<int> splitBaseIds = splitBaseIdsMap.second;

                if ( this->hasNearNullSpace( splitBaseIds ) )
                {
                    IS isToApply;
                    std::string splitIdStr = (boost::format("%1%")%splitId).str();
                    this->check( PCFieldSplitGetIS( M_pc,splitIdStr.c_str(),&isToApply ) );
                    auto const& nearnullspace = this->nearNullSpace( splitBaseIds );
                    int dimNullSpace = nearnullspace->size();
                    std::vector<Vec> petsc_vec(dimNullSpace);
                    for ( int k = 0 ; k<dimNullSpace ; ++k )
                        petsc_vec[k] =  dynamic_cast<const VectorPetsc<T>*>( &nearnullspace->basisVector(k) )->vec();

                    MatNullSpace nullsp;
                    this->check( MatNullSpaceCreate( this->worldComm(), PETSC_FALSE , dimNullSpace, petsc_vec.data(), &nullsp ) );
                    //this->check( PetscObjectCompose((PetscObject) isToApply, "nullspace", (PetscObject) nullsp) );
                    this->check( PetscObjectCompose((PetscObject) isToApply, "nearnullspace", (PetscObject) nullsp) );
                    nullspList.push_back( nullsp );
                }
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
        }
        else if (this->M_preconditioner_type==FEELPP_BLOCKMS_PRECOND)
        {
            check( PCSetType( M_pc, "blockms" ) );
            this->init();
        }
        M_indexSplitHasChanged = true;

        this->M_mat_has_changed = false;
    }

    //check( PCSetOperators( M_pc,M_mat,M_mat, PetscGetMatStructureEnum(MatrixStructure::SAME_NONZERO_PATTERN) ) );
    //check( PCSetOperators( M_pc,M_mat,M_mat, PetscGetMatStructureEnum(MatrixStructure::DIFFERENT_NONZERO_PATTERN) ) );
#if PETSC_VERSION_LESS_THAN(3,5,0)
    check( PCSetOperators( M_pc,M_mat,M_mat, PetscGetMatStructureEnum(this->M_prec_matrix_structure) ) );
#else
    check( PCSetReusePreconditioner(M_pc,(this->M_prec_matrix_structure == Feel::SAME_PRECONDITIONER)? PETSC_TRUE : PETSC_FALSE ) );
    check( PCSetOperators( M_pc,M_mat,M_mat ) );
#endif

    // Set the PCType.  Note: this used to be done *before* the call to
    // PCSetOperators(), and only when !M_is_initialized, but
    // 1.) Some preconditioners (those employing sub-preconditioners,
    // for example) have to call PCSetUp(), and can only do this after
    // the operators have been set.
    // 2.) It should be safe to call set_petsc_preconditioner_type()
    // multiple times.
    VLOG(2) << "prec : "  << this->M_preconditioner_type << "\n";
    setPetscPreconditionerType( this->M_preconditioner_type,this->M_matSolverPackage_type,M_pc,this->worldComm(),this->name() );

    // destroy null space used
    for ( MatNullSpace nullsp : nullspList )
        PETSc::MatNullSpaceDestroy( nullsp );

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
void
PreconditionerPetsc<T>::setPrecMatrixStructure( MatrixStructure mstruct  )
{
    super_type::setPrecMatrixStructure( mstruct );

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,5,0)
    if ( this-> M_is_initialized )
    {
        check( PCSetReusePreconditioner(M_pc,(this->M_prec_matrix_structure == Feel::SAME_PRECONDITIONER)? PETSC_TRUE : PETSC_FALSE ) );
    }
#endif
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
        PetscPCFactorSetMatSolverPackage( pc, matSolverPackage_type );
        break;

    case ICC_PRECOND:
        ierr = PCSetType ( pc, ( char* ) PCICC );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case ILU_PRECOND:
    {
        if ( matSolverPackage_type == MATSOLVER_PETSC )
        {
            ierr = PCSetType ( pc, ( char* ) PCILU );
            CHKERRABORT( worldComm.globalComm(),ierr );
        }
#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 ) //#ifdef FEELPP_HAS_PETSC_HYPRE
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
        else if ( worldComm.globalSize() == 1 )
        {
            ierr = PCSetType ( pc, ( char* ) PCILU );
            CHKERRABORT( worldComm.globalComm(),ierr );
        }
        else
        {
#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 ) //#ifdef FEELPP_HAS_PETSC_HYPRE
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
        if ( worldComm.globalSize() == 1 ||
             matSolverPackage_type == MATSOLVER_MUMPS ||
             matSolverPackage_type == MATSOLVER_PASTIX )
        {
            ierr = PCSetType ( pc, ( char* ) PCLU );
            CHKERRABORT( worldComm.globalComm(),ierr );

            // set factor package
            PetscPCFactorSetMatSolverPackage( pc, matSolverPackage_type );
            //ierr = PCFactorSetMatSolverPackage( pc, ( char* ) matSolverPackage_type );
            //CHKERRABORT( worldComm.globalComm(),ierr );
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

    case LSC2_PRECOND:
        ierr = PCSetType( pc, "lsc2" );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case PMM_PRECOND:
        ierr = PCSetType( pc, "pmm" );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case PCD_PRECOND:
        ierr = PCSetType( pc, "pcd" );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case FEELPP_BLOCKNS_PRECOND:
        ierr = PCSetType( pc, "blockns" );
        CHKERRABORT( worldComm.globalComm(),ierr );
        break;

    case FEELPP_BLOCKMS_PRECOND:
        ierr = PCSetType( pc, "blockms" );
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
#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 ) // #ifdef FEELPP_HAS_PETSC_HYPRE
        ierr = PCSetType( pc,( char* ) PCHYPRE );
        CHKERRABORT( worldComm.globalComm(),ierr );
        ierr = PCHYPRESetType( pc, "boomeramg" );
        CHKERRABORT( worldComm.globalComm(),ierr );
#else
        LOG(ERROR) << "preconditioner boomeramg is available with HYPRE package (petsc >= 3.3.0)";
#endif
        break;
    case AMS_PRECOND:
#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,6,0 )
        ierr = PCSetType( pc,( char* ) PCHYPRE );
        CHKERRABORT( worldComm.globalComm(),ierr );
        ierr = PCHYPRESetType( pc, "ams" );
        CHKERRABORT( worldComm.globalComm(),ierr );
#else
        LOG(ERROR) << "preconditioner ams is available with HYPRE package (petsc >= 3.6.0)";
#endif
        break;

    case REDUNDANT_PRECOND:
        ierr = PCSetType( pc,( char* ) PCREDUNDANT );
        CHKERRABORT( worldComm.globalComm(),ierr );
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
                                                          WorldComm const& worldComm,
                                                          std::string const& name )
{
    int ierr = 0;
    SetPCType( pc, preconditioner_type, matSolverPackage_type, worldComm );

    // configure main preconditioner
    ConfigurePC( pc, this, worldComm, "", name );

    // prepare PC to use
    ierr = PCSetUp( pc );
    CHKERRABORT( worldComm.globalComm(),ierr );

    if ( boption( _prefix=name, _name="pc-view" ) )
    {
        ierr = PCView( pc, PETSC_VIEWER_STDOUT_WORLD );
        CHKERRABORT( worldComm.globalComm(),ierr );
    }
}

/**
 * ConfigurePC
 */
ConfigurePC::ConfigurePC( PC& pc, PreconditionerPetsc<double> * precFeel,
                          WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( precFeel,worldComm, sub, prefix ),
    M_useConfigDefaultPetsc( option(_name="pc-use-config-default-petsc",_prefix=prefix,_sub=sub).as<bool>() ),
    M_factorShiftType( "none"/*option(_name="pc-factor-shift-type",_prefix=prefix,_sub=sub).as<std::string>()*/ )
{
    run( pc );
}
 ConfigurePC::ConfigurePC( PC& pc, PreconditionerPetsc<double> * precFeel,
                          WorldComm const& worldComm, std::string const& sub, std::string const& prefix,
                          std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( precFeel,worldComm, sub, prefix, prefixOverwrite ),
    M_useConfigDefaultPetsc( option(_name="pc-use-config-default-petsc",_prefix=prefix,_sub=sub).as<bool>() ),
    M_factorShiftType( "none"/*option(_name="pc-factor-shift-type",_prefix=prefix,_sub=sub).as<std::string>()*/ )
{
    run( pc );
}
ConfigurePC::ConfigurePC( PC& pc, PreconditionerPetsc<double> * precFeel,
                          WorldComm const& worldComm, std::string const& sub, std::string const& prefix,
                          std::vector<std::string> const& prefixOverwrite,
                          po::variables_map const& vm )
    :
    ConfigurePCBase( precFeel,worldComm, sub, prefix, prefixOverwrite ),
    M_useConfigDefaultPetsc( option(_name="pc-use-config-default-petsc",_prefix=prefix,_sub=sub,_vm=vm).as<bool>() ),
    M_factorShiftType( option(_name="pc-factor-shift-type",_prefix=prefix,_sub=sub,_vm=vm).as<std::string>() )
{
    run( pc );
}

ConfigurePC::ConfigurePC( PreconditionerPetsc<double> * precFeel,
                          WorldComm const& worldComm, std::string const& sub, std::string const& prefix,
                          std::vector<std::string> const& prefixOverwrite,
                          po::variables_map const& vm )
    :
    ConfigurePCBase( precFeel,worldComm, sub, prefix, prefixOverwrite ),
    M_useConfigDefaultPetsc( option(_name="pc-use-config-default-petsc",_prefix=prefix,_sub=sub,_vm=vm).as<bool>() ),
    M_factorShiftType( option(_name="pc-factor-shift-type",_prefix=prefix,_sub=sub,_vm=vm).as<std::string>() )
{}

void
ConfigurePC::run( PC& pc )
{
    VLOG(2) << "configuring PC... (sub: " << this->sub() << ")";
    //PCSetOptionsPrefix( pc, (this->prefix()+"_").c_str());
    google::FlushLogFiles(google::INFO);
    const char* pctype;
    this->check( PCGetType ( pc, &pctype ) );
    VLOG(2) << "configuring PC (" << this->prefix() << "." << this->sub() << ")" << pctype <<  "\n";
    google::FlushLogFiles(google::INFO);

    if ( M_useConfigDefaultPetsc )
        return;

    // init with petsc option if given and not interfaced
    if ( true )
        this->check( PCSetFromOptions( pc ) );

    if ( std::string(pctype) == "lu" || std::string(pctype) == "ilu" )
    {
        if ( M_factorShiftType == "nonzero" )
            this->check( PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO) );
        else if ( M_factorShiftType == "positive_definite" )
            this->check( PCFactorSetShiftType(pc,MAT_SHIFT_POSITIVE_DEFINITE) );
        else if ( M_factorShiftType == "inblocks" )
            this->check( PCFactorSetShiftType(pc,MAT_SHIFT_INBLOCKS) );
    }


    if ( std::string(pctype) == "gasm" )
    {
        ConfigurePCGASM( pc, this->precFeel(), this->worldComm(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( std::string(pctype) == "asm" )
    {
        ConfigurePCASM( pc, this->precFeel(), this->worldComm(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( std::string(pctype) == "bjacobi" || std::string(pctype) == "block_jacobi" )
    {
        ConfigureSubPC( pc, this->precFeel(), this->worldComm().subWorldCommSeq(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( ( std::string(pctype) == "lu" ) ||
              ( std::string(pctype) == "cholesky" ) )
    {
        ConfigurePCLU( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( std::string(pctype) == "ilu" )
    {
        ConfigurePCILU( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( std::string(pctype) == "sor" )
    {
        ConfigurePCSOR( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix(), this->prefixOverwrite() );
    }
    else if ( std::string(pctype) == "gamg" )
    {
        ConfigurePCGAMG( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix() );
    }
    else if ( std::string(pctype) == "ml" )
    {
        ConfigurePCML( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix() );
    }
    else if ( std::string(pctype) == "fieldsplit" )
    {
        ConfigurePCFieldSplit( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix() );
    }
    else if ( std::string(pctype) == "lsc" )
    {
        ConfigurePCLSC( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix() );
    }
    else if ( std::string(pctype) == "lsc2" )
    {
        ConfigurePCLSC( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix(), "in-house" );
    }
    else if ( std::string(pctype) == "pmm" )
    {
        ConfigurePCPMM( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix() );
    }
    else if ( std::string(pctype) == "pcd" )
    {
        ConfigurePCPCD( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix() );
    }
    else if ( std::string(pctype) == "blockns" )
    {
        CHECK( this->precFeel()->hasInHousePreconditioners( "blockns") ) << "blockns in-house prec not attached";
        this->check( PCSetPrecond_FEELPP(pc, this->precFeel()->inHousePreconditioners( "blockns") ) );
    }
    else if ( std::string(pctype) == "blockms" )
    {
        CHECK( this->precFeel()->hasInHousePreconditioners( "blockms") ) << "blockms in-house prec not attached";
        this->check( PCSetPrecond_FEELPP(pc, this->precFeel()->inHousePreconditioners( "blockms") ) );
    }
    else if ( std::string(pctype) == "hypre" )
    {
#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
        const char* hypretype;
        this->check( PCHYPREGetType( pc, &hypretype ) );
        if ( std::string( hypretype ) == "boomeramg" )
            ConfigurePCHYPRE_BOOMERAMG( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix(), this->prefixOverwrite() );
#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_LESS_THAN( 3,6,0 )
        else if ( std::string( hypretype ) == "euclid" )
            ConfigurePCHYPRE_EUCLID( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix(), this->prefixOverwrite() );
#endif
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,6,0 )
        else if ( std::string( hypretype ) == "ams" )
            ConfigurePCHYPRE_AMS( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix(), this->prefixOverwrite() );
#endif
#else
        LOG(ERROR) << "hypre package (euclid, boomeramg, ams) is available from PETSc version >= 3.3  (3.6 for ams)";
#endif

    }
    else if ( std::string(pctype) == "redundant" )
    {
        ConfigurePCRedundant mypc( this->precFeel(), /*pc,*/ this->worldComm(), this->prefix(), this->prefixOverwrite() );
        mypc.setFactorShiftType( M_factorShiftType );
        mypc.run(pc);
    }

    VLOG(2) << "configuring PC " << pctype << " done\n";
    google::FlushLogFiles(google::INFO);
}



void
updateOptionsDescKSP( po::options_description & _options, std::string const& prefix, std::string const& sub, bool useDefaultValue=true,
                   std::string const& kspType = "gmres", double rtol = 1e-8, size_type maxit=1000 )
{
    std::string kspctx = (sub.empty())? "" : sub+"-";

    _options.add_options()
        ( prefixvm( prefix,kspctx+"ksp-type" ).c_str(),
          (useDefaultValue)?Feel::po::value<std::string>()->default_value( kspType ):Feel::po::value<std::string>(),
          "cg, bicgstab, gmres,preonly,..." )
        ( prefixvm( prefix,kspctx+"ksp-view" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( false ):Feel::po::value<bool>(),
          "Prints the KSP data structure" )
        ( prefixvm( prefix,kspctx+"ksp-monitor" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( false ):Feel::po::value<bool>(),
          "monitor ksp" )
        ( prefixvm( prefix,kspctx+"ksp-converged-reason" ).c_str() , "converged reason ksp" )
        ( prefixvm( prefix,kspctx+"ksp-verbose" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( 0 ):Feel::po::value<int>(),
          "(=0,1,2) print solver iterations" )
        ( prefixvm( prefix,kspctx+"ksp-rtol" ).c_str(),
          (useDefaultValue)?Feel::po::value<double>()->default_value( rtol ):Feel::po::value<double>(),
          "relative tolerance" )
        ( prefixvm( prefix,kspctx+"ksp-atol" ).c_str(),
          (useDefaultValue)?Feel::po::value<double>()->default_value( 1e-50 ):Feel::po::value<double>(),
          "absolute tolerance" )
        ( prefixvm( prefix,kspctx+"ksp-dtol" ).c_str(),
          (useDefaultValue)?Feel::po::value<double>()->default_value( 1e5 ):Feel::po::value<double>(),
          "divergence tolerance" )
        ( prefixvm( prefix,kspctx+"ksp-maxit" ).c_str(),
          (useDefaultValue)?Feel::po::value<size_type>()->default_value( maxit ):Feel::po::value<size_type>(),
          "maximum number of iterations" )
        ( prefixvm( prefix,kspctx+"ksp-maxit-reuse" ).c_str(),
          (useDefaultValue)?Feel::po::value<size_type>():Feel::po::value<size_type>(),
          "maximum number of iterations when reuse prec/jac" )
        ( prefixvm( prefix,kspctx+"constant-null-space" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( 0 ):Feel::po::value<bool>(),
          "set the null space to be the constant values" )
        ( prefixvm( prefix,kspctx+"ksp-use-config-default-petsc" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( false ):Feel::po::value<bool>(),
          "configure ksp with default petsc options" )
        ( prefixvm( prefix,kspctx+"gmres-restart" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( 30 ):Feel::po::value<int>(),
          "number of iterations before solver restarts (gmres)" )
        ;
}
po::options_description
getOptionsDescKSP( std::string const& prefix, std::string const& sub, bool useDefaultValue=true,
                   std::string const& kspType = "gmres", double rtol = 1e-8, size_type maxit=1000 )
{
    po::options_description _options( "options KSP",200);
    updateOptionsDescKSP( _options, prefix, sub, useDefaultValue, kspType, rtol, maxit );
    return _options;
}

po::options_description
getOptionsDescKSP( std::string const& prefix, std::string const& sub, std::vector<std::string> const& prefixOverwrite,
                   std::string const& kspType = "gmres", double rtol = 1e-8, size_type maxit=1000 )
{
    po::options_description _options( "options KSP",200);
    updateOptionsDescKSP( _options, prefix, sub, true, kspType, rtol, maxit );
    for ( std::string const& prefixOver : prefixOverwrite )
        updateOptionsDescKSP( _options, prefixOver, sub, false );
    return _options;
}
po::options_description
updateOptionsDescPrecBase( po::options_description & _options, std::string const& prefix, std::string const& sub, bool useDefaultValue=true, std::string const& pcType = "lu" )
{
    std::string pcctx = (sub.empty())? "" : sub+"-";
    //po::options_description _options( "options PC", 200);

    _options.add_options()
        ( prefixvm( prefix,pcctx+"pc-type" ).c_str(),
          (useDefaultValue)?Feel::po::value<std::string>()->default_value( pcType ):Feel::po::value<std::string>(),
          "type of preconditioners (lu, cholesky, icc, ilut, ilutp, diag, id,...)" )
        ( prefixvm( prefix,pcctx+"pc-view" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( false ):Feel::po::value<bool>(),
          "display preconditioner information" )
//#if defined(FEELPP_HAS_MUMPS) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
#if defined(PETSC_HAVE_MUMPS)
        ( prefixvm( prefix,pcctx+"pc-factor-mat-solver-package-type" ).c_str(),
          (useDefaultValue)?Feel::po::value<std::string>()->default_value( "mumps" ):Feel::po::value<std::string>(),
          "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#else
        ( prefixvm( prefix,pcctx+"pc-factor-mat-solver-package-type" ).c_str(),
          (useDefaultValue)?Feel::po::value<std::string>()->default_value( "petsc" ):Feel::po::value<std::string>(),
          "sets the software that is used to perform the factorization (petsc,umfpack, spooles, petsc, superlu, superlu_dist, mumps,...)" )
#endif
        ( prefixvm( prefix,pcctx+"pc-use-config-default-petsc" ).c_str(),
          (useDefaultValue)?Feel::po::value<bool>()->default_value( false ):Feel::po::value<bool>(),
          "configure pc with default petsc options" )
        ( prefixvm( prefix,pcctx+"pc-factor-shift-type" ).c_str(),
          (useDefaultValue)?Feel::po::value<std::string>()->default_value( "none" ):Feel::po::value<std::string>(),
          "adds a particular type of quantity to the diagonal of the matrix during numerical factorization, thus the matrix has nonzero pivots : none, nonzero, positive_definite, inblocks" )
        ;

    return _options;
}
po::options_description
getOptionsDescPrecBase( std::string const& prefix, std::string const& sub, bool useDefaultValue=true, std::string const& pcType = "lu" )
{
    po::options_description _options( "options PC", 200);
    updateOptionsDescPrecBase( _options,prefix,sub,useDefaultValue,pcType );
    //for ( std::string prefixOver : prefixOverwrite )
    //    updateOptionsDescPrecBase( _options,prefixOver,sub );
    return _options;

}
po::options_description
getOptionsDescPrecBase( std::string const& prefix, std::string const& sub, std::vector<std::string> const& prefixOverwrite,
                        std::string const& pcType = "lu" )
{
    po::options_description _options( "options PC", 200);
    updateOptionsDescPrecBase( _options,prefix,sub,true,pcType );
    for ( std::string const& prefixOver : prefixOverwrite )
        updateOptionsDescPrecBase( _options,prefixOver,sub,false );
    return _options;

}
void
updateOptionsDescLU( po::options_description & _options, std::string const& prefix, std::string const& sub )
{
    std::string pcctx = (sub.empty())? "" : sub+"-";

//#if defined(FEELPP_HAS_MUMPS) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
#if defined(PETSC_HAVE_MUMPS)
    for ( int icntl=1 ; icntl<= 33 ; ++icntl )
    {
        std::string mumpsOption = (boost::format("pc-factor-mumps.icntl-%1%")%icntl ).str();
        if( icntl == 7 )
            _options.add_options()
                ( prefixvm( prefix,pcctx+mumpsOption ).c_str(),
                  Feel::po::value<int>()->default_value( 0 ),"configure mumps factorisation (see mumps ICNTL documentation)" );
        else
            _options.add_options()
                ( prefixvm( prefix,pcctx+mumpsOption ).c_str(),
                  Feel::po::value<int>(),"configure mumps factorisation (see mumps ICNTL documentation)" );
    }
#endif
}
po::options_description
getOptionsDescLU( std::string const& prefix, std::string const& sub, std::vector<std::string> const& prefixOverwrite )
{
    po::options_description _options( "options PC LU", 200);
    updateOptionsDescLU( _options,prefix,sub );
    for ( std::string const& prefixOver : prefixOverwrite )
        updateOptionsDescLU( _options,prefixOver,sub );
    return _options;
}

void
updateOptionsDescILU( po::options_description & _options, std::string const& prefix, std::string const& sub, bool useDefaultValue=true )
{
    std::string pcctx = (sub.empty())? "" : sub+"-";
    _options.add_options()
        ( prefixvm( prefix,pcctx+"pc-factor-levels" ).c_str(),
        (useDefaultValue)?Feel::po::value<int>()->default_value( 3 ):Feel::po::value<int>(),
        "Sets the number of levels of fill to use for ilu" )
        ( prefixvm( prefix,pcctx+"pc-factor-fill" ).c_str(),
        (useDefaultValue)?Feel::po::value<double>()->default_value( 6 ):Feel::po::value<double>(),
        "Indicate the amount of fill you expect in the factored matrix, fill = number nonzeros in factor/number nonzeros in original matrix." )
        ;
}
po::options_description
getOptionsDescILU( std::string const& prefix, std::string const& sub, std::vector<std::string> const& prefixOverwrite )
{
    po::options_description _options( "options PC ILU", 200);
    updateOptionsDescILU( _options,prefix,sub,true );
    for ( std::string const& prefixOver : prefixOverwrite )
        updateOptionsDescILU( _options,prefixOver,sub,false );
    return _options;
}

void
updateOptionsDescSOR( po::options_description & _options, std::string const& prefix, std::string const& sub, bool useDefaultValue=true )
{
    std::string pcctx = (sub.empty())? "" : sub+"-";
    _options.add_options()
        ( prefixvm( prefix,pcctx+"pc-sor-omega" ).c_str(),
          (useDefaultValue)?Feel::po::value<double>()->default_value( 1. ):Feel::po::value<double>(),
          "Sets the SOR relaxation coefficient, omega" )
        ( prefixvm( prefix,pcctx+"pc-sor-lits" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( 1 ):Feel::po::value<int>(),
          "number of local iterations, smoothings over just variables on processor" )
        ( prefixvm( prefix,pcctx+"pc-sor-its" ).c_str(),
          (useDefaultValue)?Feel::po::value<int>()->default_value( 1 ):Feel::po::value<int>(),
          "number of parallel iterations to use; each parallel iteration has lits local iterations" )
        ( prefixvm( prefix,pcctx+"pc-sor-type" ).c_str(),
          (useDefaultValue)?Feel::po::value<std::string>()->default_value( "local_symmetric" ):Feel::po::value<std::string>(),
          "(symmetric,forward,backward,local_symmetric,local_forward,local_backward) Sets the SOR preconditioner to use symmetric (SSOR), backward, or forward relaxation. The local variants perform SOR on each processor" )
        ;
}
po::options_description
getOptionsDescSOR( std::string const& prefix, std::string const& sub, std::vector<std::string> const& prefixOverwrite )
{
    po::options_description _options( "options PC SOR", 200);
    updateOptionsDescSOR( _options,prefix,sub,true );
    for ( std::string const& prefixOver : prefixOverwrite )
        updateOptionsDescSOR( _options,prefixOver,sub,false );
    return _options;
}

void
updateOptionsDescGASM( po::options_description & _options, std::string const& prefix, bool useDefaultValue=true )
{
    _options.add_options()
        ( prefixvm( prefix,"pc-gasm-type" ).c_str(), (useDefaultValue)?Feel::po::value<std::string>()->default_value( "restrict" ):Feel::po::value<std::string>(),
          "type of gasm (basic, restrict, interpolate, none)" )
        ( prefixvm( prefix,"pc-gasm-overlap" ).c_str(), (useDefaultValue)?Feel::po::value<int>()->default_value( 1 ):Feel::po::value<int>(),
          "number of overlap levels" )
        ;
}
po::options_description
getOptionsDescGASM( std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
{
    po::options_description _options( "options PC GASM", 100);
    updateOptionsDescGASM( _options,prefix,true);
    for ( std::string const& prefixOver : prefixOverwrite )
        updateOptionsDescGASM( _options,prefixOver,false);
    return _options;
}

void
updateOptionsDescBOOMERAMG( po::options_description & _options, std::string const& prefix, bool useDefaultValue=true )
{
    _options.add_options()
        /// Documentation here : http://computation.llnl.gov/project/linear_solvers/download/hypre-2.10.0b_ref_manual.pdf
        ( prefixvm( prefix,"pc-hypre-boomeramg-max-iter" ).c_str(),Feel::po::value<int>()->default_value( 20 ),"sets the max number of iteration for boomeramg" )
        ( prefixvm( prefix,"pc-hypre-boomeramg-tol").c_str(), Feel::po::value<double>()->default_value( 1e-7 ), "Convergence tolerance PER hypre call (0.0 = use a fixed number of iterations)")
        ( prefixvm( prefix,"pc-hypre-boomeramg-cycle-type").c_str(), Feel::po::value<int>()->default_value( 0 ),"Cycle type (0=V / 11==W")
        ( prefixvm( prefix,"pc-hypre-boomeramg-print-statistics").c_str(),Feel::po::value<bool>()->default_value( false ),"Print statistics")
        ( prefixvm( prefix,"pc-hypre-boomeramg-max-levels").c_str(),Feel::po::value<int>()->default_value( 25 ),"Number of levels (of grids) allowed")
        ( prefixvm( prefix,"pc-hypre-boomeramg-strong-threshold").c_str(),Feel::po::value<double>()->default_value( 0.25 ),"Threshold for being strongly connected")
        ( prefixvm( prefix,"pc-hypre-boomeramg-max-row-sum").c_str(),Feel::po::value<double>()->default_value( 0.9 ),"Maximum row sum")
        ( prefixvm( prefix,"pc-hypre-boomeramg-coarsen-type").c_str(),Feel::po::value<int>()->default_value( 6 ),"Coarsen type")
        ( prefixvm( prefix,"pc-hypre-boomeramg-agg-nl").c_str(),Feel::po::value<int>()->default_value( 1 ),"Number of levels of aggressive coarsening")
        ( prefixvm( prefix,"pc-hypre-boomeramg-relax-type-all").c_str(),Feel::po::value<int>()->default_value( 6 ),"Relax type for the up and down cycles")
        ( prefixvm( prefix,"pc-hypre-boomeramg-interp-type").c_str(),Feel::po::value<int>()->default_value( 6 ),"Interpolation type")
#if 0 // Not interfaced options
        ( prefixvm( prefix,"pc-hypre-boomeramg-max-levels").c_str(),
        Feel::po::value<int>()->default_value( 2 ),
        "Number of levels (of grids) allowed")
        ( prefixvm( prefix,"pc-hypre-boomeramg-tol").c_str(),
        Feel::po::value<double>()->default_value( 1e-7 ),
        "Convergence tolerance PER hypre call (0.0 = use a fixed number of iterations)")
        ( prefixvm( prefix,"pc_hypre_boomeramg_truncfactor").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Truncation factor for interpolation (0=no truncation)")
        ( prefixvm( prefix,"pc_hypre_boomeramg_P_max").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Max elements per row for interpolation operator (0=unlimited)")
        ( prefixvm( prefix,"pc_hypre_boomeramg_agg_nl").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Number of levels of aggressive coarsening")
        ( prefixvm( prefix,"pc_hypre_boomeramg_agg_num_paths").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Number of paths for aggressive coarsening")
        ( prefixvm( prefix,"pc-hypre-boomeramg-strong-threshold").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Threshold for being strongly connected")
        ( prefixvm( prefix,"pc-hypre-boomeramg-max-row-sum").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Maximum row sum")
        ( prefixvm( prefix,"pc_hypre_boomeramg_grid_sweeps_all").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Number of sweeps for the up and down grid levels")
        ( prefixvm( prefix,"pc_hypre_boomeramg_grid_sweeps_down").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Number of sweeps for the down cycles")
        ( prefixvm( prefix,"pc_hypre_boomeramg_grid_sweeps_up").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Number of sweeps for the up cycles")
        ( prefixvm( prefix,"pc_hypre_boomeramg_grid_sweeps_coarse").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Number of sweeps for the coarse level")
        ( prefixvm( prefix,"pc_hypre_boomeramg_relax_type_all").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Relax type for the up and down cycles")
        ( prefixvm( prefix,"pc_hypre_boomeramg_relax_type_down").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Relax type for the down cycles")
        ( prefixvm( prefix,"pc_hypre_boomeramg_relax_type_up").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Relax type for the up cycles")
        ( prefixvm( prefix,"pc_hypre_boomeramg_relax_type_coarse").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Relax type on coarse grid")
        ( prefixvm( prefix,"pc_hypre_boomeramg_relax_weight_all").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Relaxation weight for all levels (0 = hypre estimates, -k = determined with k CG steps)")
        ( prefixvm( prefix,"pc_hypre_boomeramg_relax_weight_level").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Set the relaxation weight for a particular level (weight,level)")
        ( prefixvm( prefix,"pc_hypre_boomeramg_outer_relax_weight_all").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Outer relaxation weight for all levels (-k = determined with k CG steps)")
        ( prefixvm( prefix,"pc_hypre_boomeramg_outer_relax_weight_level").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Set the outer relaxation weight for a particular level (weight,level)")
        ( prefixvm( prefix,"pc_hypre_boomeramg_no_CF").c_str(), 
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Do not use CF-relaxation")
        ( prefixvm( prefix,"pc_hypre_boomeramg_measure_type").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Measure type")
        ( prefixvm( prefix,"pc-hypre-boomeramg-coarsen-type").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Coarsen type"
        ( prefixvm( prefix,"pc_hypre_boomeramg_interp_type").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Interpolation type")
        ( prefixvm( prefix,"pc-hypre-boomeramg-print-statistics").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Print statistics")
        ( prefixvm( prefix,"pc_hypre_boomeramg_print_debug").c_str(),
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Print debug information")
        ( prefixvm( prefix,"pc_hypre_boomeramg_nodal_coarsen").c_str(), 
        Feel::po::value<double>()->default_value( 1e-6 ),
        "HYPRE_BoomerAMGSetNodal()")
        ( prefixvm( prefix,"pc_hypre_boomeramg_nodal_relaxation").c_str(), 
        Feel::po::value<double>()->default_value( 1e-6 ),
        "Nodal relaxation via Schwarz")
#endif
        ;
}
po::options_description
getOptionsDescBOOMERAMG( std::string const& prefix, std::string const& sub, std::vector<std::string> const& prefixOverwrite )
{
    po::options_description _options( "options PC BOOMERAMG", 100);
    updateOptionsDescBOOMERAMG(_options,prefix,true);
    for ( std::string const& prefixOver : prefixOverwrite )
        updateOptionsDescBOOMERAMG( _options,prefixOver,false);
    return _options;
}

void
updateOptionsDescAMS( po::options_description & _options, std::string const& prefix, bool useDefaultValue=true )
{
    _options.add_options()
        /// Documentation here : http://computation.llnl.gov/project/linear_solvers/download/hypre-2.10.0b_ref_manual.pdf p46-53
        ( prefixvm( prefix,"pc-hypre-ams-print-level").c_str(),useDefaultValue?Feel::po::value<int>()->default_value( 1 ):Feel::po::value<int>(),       "Debugging output level for AMS" )
        ( prefixvm( prefix,"pc-hypre-ams-max-iter").c_str(),useDefaultValue?Feel::po::value<int>()->default_value( 1 ):Feel::po::value<int>(),          "Maximum number of AMS multigrid iterations within PCApply" )
        ( prefixvm( prefix,"pc-hypre-ams-cycle-type").c_str(),useDefaultValue?Feel::po::value<int>()->default_value( 13 ):Feel::po::value<int>(),       "Cycle type for AMS multigrid" )
        ( prefixvm( prefix,"pc-hypre-ams-tol").c_str(),useDefaultValue?Feel::po::value<double>()->default_value( 0. ):Feel::po::value<double>(),        "Error tolerance for AMS multigrid" )
        ( prefixvm( prefix,"pc-hypre-ams-relax-type").c_str(),useDefaultValue?Feel::po::value<int>()->default_value( 2 ):Feel::po::value<int>(),        "Relaxation type for AMS smoother" )
        ( prefixvm( prefix,"pc-hypre-ams-relax-times").c_str(),useDefaultValue?Feel::po::value<int>()->default_value( 1 ):Feel::po::value<int>(),       "Number of relaxation steps for AMS smoother" )
        ( prefixvm( prefix,"pc-hypre-ams-relax-weight").c_str(),useDefaultValue?Feel::po::value<double>()->default_value( 1 ):Feel::po::value<double>(),"Relaxation weight for AMS smoother" )
        ( prefixvm( prefix,"pc-hypre-ams-omega").c_str(),useDefaultValue?Feel::po::value<double>()->default_value( 1 ):Feel::po::value<double>(),       "SSOR coefficient for AMS smoother" )
#if 0 // Not interfaced yet
        ( prefixvm( prefix,"pc-hypre-ams-amg-alpha-theta").c_str(),Feel::po::value<double>()->default_value( 20 ),"Threshold for strong coupling of vector Poisson AMG solver" )
        ( prefixvm( prefix,"pc-hypre-ams-amg-alpha-options").c_str(),Feel::po::value<XXX>()->default_value( 20 ), "AMG options for vector Poisson" ) // array of 5 floats
        ( prefixvm( prefix,"pc-hypre-ams-amg-beta-theta").c_str(),Feel::po::value<int>()->default_value( 20 ),    "Threshold for strong coupling of scalar Poisson AMG solver" )
        ( prefixvm( prefix,"pc-hypre-ams-amg-beta-options").c_str(),Feel::po::value<XXX>()->default_value( 20 ),  "AMG options for scalar Poisson solver" ) // array of 5 floats
#endif
        ;
}
po::options_description
getOptionsDescAMS( std::string const& prefix, std::string const& sub, std::vector<std::string> const& prefixOverwrite )
{
    po::options_description _options( "options PC AMS", 100);
    updateOptionsDescAMS(_options,prefix,true);
    for ( std::string const& prefixOver : prefixOverwrite )
        updateOptionsDescAMS( _options,prefixOver,false);
    return _options;
}

void
updateOptionsDescASM( po::options_description & _options, std::string const& prefix, bool useDefaultValue=true )
{
    _options.add_options()
        ( prefixvm( prefix,"pc-asm-type" ).c_str(), (useDefaultValue)?Feel::po::value<std::string>()->default_value( "restrict" ):Feel::po::value<std::string>(),
          "type of asm (basic, restrict, interpolate, none)" )
        ( prefixvm( prefix,"pc-asm-overlap" ).c_str(), (useDefaultValue)?Feel::po::value<int>()->default_value( 1 ):Feel::po::value<int>(),
          "number of overlap levels" )
        ;
}
po::options_description
getOptionsDescASM( std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
{
    po::options_description _options( "options PC ASM", 100);
    updateOptionsDescASM(_options,prefix,true);
    for ( std::string const& prefixOver : prefixOverwrite )
        updateOptionsDescASM( _options,prefixOver,false);
    return _options;
}

po::options_description
getOptionsDescML( std::string const& prefix, std::string const& sub )
{
    std::string pcctx = (sub.empty())? "pc-" : sub+"-pc-";

    po::options_description _options( "options PC ML", 200);
    // multigrid options
    _options.add_options()
        ( prefixvm( prefix,pcctx+"mg-levels" ).c_str(), Feel::po::value<int>()->default_value( 10/*2*/ ),
          "number of levels including finest" )
        ( prefixvm( prefix,pcctx+"mg-type" ).c_str(), Feel::po::value<std::string>()->default_value( "multiplicative" ),
          "Determines the form of multigrid to use: multiplicative, additive, full, kaskade " )
        ( prefixvm( prefix,pcctx+"mg-smoothdown" ).c_str(), Feel::po::value<int>()->default_value( 1 ),
          "number of smoothing steps before applying restriction operator" )
        ;
    // ml options
    _options.add_options()
        ( prefixvm( prefix,pcctx+"ml-reuse-interpolation" ).c_str(), Feel::po::value<bool>()->default_value( false ),
          "Reuse the interpolation operators when possible (cheaper, weaker when matrix entries change a lot)" )
        ( prefixvm( prefix,pcctx+"ml-keep-agg-info" ).c_str(), Feel::po::value<bool>()->default_value( false ),
          "Allows the preconditioner to be reused, or auxilliary matrices to be generated" )
        ( prefixvm( prefix,pcctx+"ml-reusable" ).c_str(), Feel::po::value<bool>()->default_value( false ),
          "Store intermedaiate data structures so that the multilevel hierarchy is reusable" )
        ( prefixvm( prefix,pcctx+"ml-old-hierarchy" ).c_str(), Feel::po::value<bool>()->default_value( false ),
          "Use old routine to generate hierarchy" )
        ;
    // coarse ksp/pc
    std::string mgctx = (sub.empty())? "mg-" : sub+"-mg-";
    std::string prefixMGCoarse = ( boost::format( "%1%%2%coarse" ) %prefixvm(prefix,"") %mgctx ).str();
    //_options.add( getOptionsDescPrecBase(prefixMGCoarse,"",true,"redundant") );
    po::options_description optionsCoarse( "options PC ML Coarse Level", 200);
    updateOptionsDescPrecBase(optionsCoarse,prefixMGCoarse,"",true,"redundant");
    _options.add( optionsCoarse );

    return _options;
}
po::options_description
getOptionsDescGAMG( std::string const& prefix, std::string const& sub )
{
    std::string pcctx = (sub.empty())? "pc-" : sub+"-pc-";

    po::options_description _options( "options PC GAMG", 200);
    // multigrid options
    _options.add_options()
        ( prefixvm( prefix,pcctx+"mg-levels" ).c_str(), Feel::po::value<int>()->default_value( 10/*2*/ ),
          "number of levels including finest" )
        ( prefixvm( prefix,pcctx+"mg-type" ).c_str(), Feel::po::value<std::string>()->default_value( "multiplicative" ),
          "Determines the form of multigrid to use: multiplicative, additive, full, kaskade " )
        ( prefixvm( prefix,pcctx+"mg-smoothdown" ).c_str(), Feel::po::value<int>()->default_value( 1 ),
          "number of smoothing steps before applying restriction operator" )
        ;
    // gamg options
    _options.add_options()
        ( prefixvm( prefix,pcctx+"gamg-type" ).c_str(), Feel::po::value<std::string>()->default_value( "agg" ),
          "type of generalized algebraic multigrid : agg, geo " )
        ( prefixvm( prefix,pcctx+"gamg-proc-eq-lim" ).c_str(), Feel::po::value<int>()->default_value( 50 ),
          "number of equations to aim for on coarse grids via processor reduction" )
        ( prefixvm( prefix,pcctx+"gamg-coarse-eq-lim" ).c_str(), Feel::po::value<int>()->default_value( 800 ),
          "max number of equations on coarse grids" )
        ( prefixvm( prefix,pcctx+"gamg-threshold" ).c_str(), Feel::po::value<double>()->default_value( 0. ),
          "relative threshold to use for dropping edges in aggregation graph" )
        ( prefixvm( prefix,pcctx+"gamg-set-sym-graph" ).c_str(), Feel::po::value<bool>()->default_value( true ),
          "set for asymmetric matrice (if the matrix is sym, put to false)" )
        ( prefixvm( prefix,pcctx+"gamg-reuse-interpolation").c_str(), Feel::po::value<bool>()->default_value( false ),
          "reuse prolongation operator" )
        ( prefixvm( prefix,pcctx+"gamg-verbose").c_str(), Feel::po::value<int>()->default_value( 0 ),
          "verbose internal petsc info for gamg" )
        ( prefixvm( prefix,pcctx+"gamg-nsmooths" ).c_str(), Feel::po::value<int>()->default_value( 1 ),
          "number of smoothing steps" )
        ;

    // coarse ksp/pc
    std::string mgctx = (sub.empty())? "mg-" : sub+"-mg-";
    std::string prefixMGCoarse = ( boost::format( "%1%%2%coarse" ) %prefixvm(prefix,"") %mgctx ).str();
    //_options.add( getOptionsDescPrecBase(prefixMGCoarse,"",true,"redundant") );
    po::options_description optionsCoarse( "options PC GAMG Coarse Level", 200);
    updateOptionsDescPrecBase(optionsCoarse,prefixMGCoarse,"",true,"redundant");
    _options.add( optionsCoarse );

    return _options;
}

po::options_description
getOptionsDescMultiGridLevels( int nLevel,std::string const& prefix, std::string const& sub )
{
    std::string mgctx = (sub.empty())? "mg-" : sub+"-mg-";
    po::options_description _options( "options PC MultiGrid Levels", 100);

    // all levels ksp/pc (not including coarse level) with default values
    std::string prefixMGLevelsGeneric = prefixvm( prefix, mgctx+"levels" );
    po::options_description optionsAllLevel( "options PC MultiGrid Levels (all level)", 200);
    updateOptionsDescPrecBase(optionsAllLevel,prefixMGLevelsGeneric,"",true,"sor");
    _options.add( optionsAllLevel );

    // fine level
    std::string prefixMGFineLevel = prefixvm( prefix, mgctx+"fine-level" );
    po::options_description optionsFineLevel( "options PC MultiGrid Levels (fine level)", 200);
    updateOptionsDescPrecBase(optionsFineLevel,prefixMGFineLevel,"",false);
    _options.add( optionsFineLevel );

    // each levels can be control separately
    po::options_description optionsEachLevel( "options PC MultiGrid Levels (each level)", 200);
    for ( uint16_type i=1; i<=nLevel; ++i )
    {
        std::string prefixMGLevels = ( boost::format( "%1%%2%levels%3%" ) %prefixvm(prefix,"") %mgctx %i ).str();
        updateOptionsDescPrecBase(optionsEachLevel,prefixMGLevels,"",false);
    }
    _options.add( optionsEachLevel );

    return _options;
}

po::options_description
getOptionsDescFieldSplit( std::string const& prefix, std::string const& sub )
{
    std::string pcctx = (sub.empty())? "" : sub+"-";

    po::options_description _options( "options PC FieldSplit", 200,200);

    // field split options
    _options.add_options()
        ( prefixvm( prefix,pcctx+"fieldsplit-type" ).c_str(), Feel::po::value<std::string>()->default_value( "additive" ),
          "type of fieldsplit (additive, multiplicative, symmetric-multiplicative, schur)" )
        //( prefixvm( prefix,pcctx+"fieldsplit-fields" ).c_str(), Feel::po::value<std::string>()->default_value( "" ),
        //  "fields definition (ex: --fieldsplit-fields=0->(0,2),1->(1)" )
        ;
    // schur complement options
    _options.add_options()
        ( prefixvm( prefix,pcctx+"fieldsplit-schur-fact-type" ).c_str(), Feel::po::value<std::string>()->default_value( "full" ),
          "type of schur factorization (diag, lower, upper, full)" )
        ( prefixvm( prefix,pcctx+"fieldsplit-schur-precondition" ).c_str(), Feel::po::value<std::string>()->default_value( "a11" ),
          "self,user,a11" )
        ;

    // inner solver (A^{-1}) of schur complement : S = C-B A^{-1} B^T
    std::string prefixSchurInnerSolver = prefixvm( prefix,pcctx+"fieldsplit-schur-inner-solver" );
    _options.add_options()
        ( prefixvm( prefixSchurInnerSolver,"use-outer-solver" ).c_str(), Feel::po::value<bool>()->default_value( true ), "use-outer-solver" )
        ;
    _options.add( getOptionsDescPrecBase(prefixSchurInnerSolver,"",true,"jacobi") );

    // solver (A^{-1}) used in upper schur preconditioning
    std::string prefixSchurUpperSolver = prefixvm( prefix,pcctx+"fieldsplit-schur-upper-solver" );
    _options.add_options()
        ( prefixvm( prefixSchurUpperSolver,"use-outer-solver" ).c_str(), Feel::po::value<bool>()->default_value( true ), "use-outer-solver" )
        ;
    _options.add( getOptionsDescPrecBase(prefixSchurUpperSolver,"",true,"jacobi") );

    return _options;
}

po::options_description
getOptionsDescSplitInFieldSplit( int nSplit, std::string const& prefix, std::string const& sub )
{
    std::string pcctx = (sub.empty())? "" : sub+"-";
    po::options_description _options( "options PC Split In FieldSplit", 100);

    // ksp/pc options for each split
    for ( uint16_type i=0; i<nSplit; ++i )
    {
        std::string prefixfieldsplit = ( boost::format( "%1%%2%fieldsplit-%3%" ) %prefixvm( prefix,"" ) %pcctx %i ).str();
        _options.add_options()
            ( prefixvm( prefixfieldsplit,"fieldsplit-fields" ).c_str(), Feel::po::value<std::string>()->default_value( "" ),
              "fields definition (ex: --fieldsplit-fields=0->(0,2),1->(1)" )
            ( prefixvm( prefixfieldsplit,"fieldsplit-use-components" ).c_str(), Feel::po::value<bool>()->default_value( false ),"split also with components" )
            ;
        _options.add( getOptionsDescPrecBase(prefixfieldsplit,"",true,(i==0)?"lu":"none") );
    }

    return _options;
}

po::options_description
getOptionsDescLSC( std::string const& prefix, std::string const& sub )
{
    std::string pcctx = (sub.empty())? "" : sub+"-";
    std::string prefixLSC = prefixvm( prefix,pcctx+"lsc");
    po::options_description _options( "options PC LSC", 100);
    _options.add_options()
        ( prefixvm( prefixLSC,"scale-diag" ).c_str(), Feel::po::value<bool>()->default_value( false ), "scale diag" )
        ;
    _options.add( getOptionsDescPrecBase(prefixLSC,"",true,"lu") );

    return _options;
}

po::options_description
getOptionsDescPMM( std::string const& prefix, std::string const& sub )
{
    std::string pcctx = (sub.empty())? "" : sub+"-";
    std::string prefixPMM = prefixvm( prefix,pcctx+"pmm");
    po::options_description _options( "options PC PMM", 100);
    _options.add( getOptionsDescPrecBase(prefixPMM,"",true,"lu") );

    return _options;
}

po::options_description
getOptionsDescPCD( std::string const& prefix, std::string const& sub )
{
    std::string pcctx = (sub.empty())? "" : sub+"-";
    std::string prefixPCD_Ap = prefixvm( prefix,pcctx+"pcd.Ap");
    std::string prefixPCD_Mp = prefixvm( prefix,pcctx+"pcd.Mp");
    po::options_description _options( "options PC PCD", 100);
    _options.add( getOptionsDescPrecBase(prefixPCD_Ap,"",true,"lu") );
    _options.add( getOptionsDescPrecBase(prefixPCD_Mp,"",true,"lu") );
    return _options;
}

/**
 * ConfigureKSP
 */
#if 0
ConfigureKSP::ConfigureKSP( KSP& ksp,WorldComm const& worldComm, std::string const& sub,std::string const& prefix )
    :
    ConfigurePCBase( worldComm,sub,prefix, getOptionsDescKSP( prefix, sub ) ),
    M_type( option(_name="ksp-type",_sub=sub,_prefix=prefix,_vm=this->vm() ).as<std::string>() ),
    M_useConfigDefaultPetsc( option(_name="ksp-use-config-default-petsc",_prefix=prefix,_sub=sub,_vm=this->vm()).as<bool>() ),
    M_rtol( option(_name="ksp-rtol",_sub=sub,_prefix=prefix,_vm=this->vm()).as<double>() ),
    M_maxit( option(_name="ksp-maxit",_sub=sub,_prefix=prefix,_vm=this->vm()).as<size_type>() ),
    M_showMonitor( option(_name="ksp-monitor",_sub=sub,_prefix=prefix,_vm=this->vm()).as<bool>() ),
    M_kspView( option(_name="ksp-view",_sub=sub,_prefix=prefix,_vm=this->vm()).as<bool>() ),
    M_constantNullSpace( option(_name="constant-null-space",_sub=sub,_prefix=prefix,_vm=this->vm()).as<bool>() )
{
    run( ksp );
}
#endif
ConfigureKSP::ConfigureKSP( KSP& ksp, PreconditionerPetsc<double> * precFeel,WorldComm const& worldComm, std::string const& sub,std::string const& prefix,
                            std::vector<std::string> const& prefixOverwrite,
                            std::string const& kspType, double rtol, size_type maxit )
    :
    ConfigurePCBase( precFeel,worldComm,sub,prefix,prefixOverwrite,getOptionsDescKSP( prefix,sub,prefixOverwrite,kspType,rtol,maxit ) ),
    M_type( getOption<std::string>("ksp-type",prefix,sub,prefixOverwrite,this->vm()) ),
    M_useConfigDefaultPetsc( getOption<bool>("ksp-use-config-default-petsc",prefix,sub,prefixOverwrite,this->vm()) ),
    M_rtol( getOption<double>("ksp-rtol",prefix,sub,prefixOverwrite,this->vm()) ),
    M_maxit( getOption<size_type>("ksp-maxit",prefix,sub,prefixOverwrite,this->vm()) ),
    M_showMonitor( getOption<bool>("ksp-monitor",prefix,sub,prefixOverwrite,this->vm()) ),
    M_kspView( getOption<bool>("ksp-view",prefix,sub,prefixOverwrite,this->vm()) ),
    M_constantNullSpace( getOption<bool>("constant-null-space",prefix,sub,prefixOverwrite,this->vm()) ),
    M_nRestartGMRES( getOption<int>("gmres-restart",prefix,sub,prefixOverwrite,this->vm()) )
{
    run( ksp );
}
ConfigureKSP::ConfigureKSP( PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm, std::string const& sub,std::string const& prefix,
                            std::vector<std::string> const& prefixOverwrite,
                            std::string const& kspType, double rtol, size_type maxit )
    :
    ConfigurePCBase( precFeel,worldComm,sub,prefix,prefixOverwrite,getOptionsDescKSP( prefix,sub,prefixOverwrite,kspType,rtol,maxit ) ),
    M_type( getOption<std::string>("ksp-type",prefix,sub,prefixOverwrite,this->vm()) ),
    M_useConfigDefaultPetsc( getOption<bool>("ksp-use-config-default-petsc",prefix,sub,prefixOverwrite,this->vm()) ),
    M_rtol( getOption<double>("ksp-rtol",prefix,sub,prefixOverwrite,this->vm()) ),
    M_maxit( getOption<size_type>("ksp-maxit",prefix,sub,prefixOverwrite,this->vm()) ),
    M_showMonitor( getOption<bool>("ksp-monitor",prefix,sub,prefixOverwrite,this->vm()) ),
    M_kspView( getOption<bool>("ksp-view",prefix,sub,prefixOverwrite,this->vm()) ),
    M_constantNullSpace( getOption<bool>("constant-null-space",prefix,sub,prefixOverwrite,this->vm()) ),
    M_nRestartGMRES( getOption<int>("gmres-restart",prefix,sub,prefixOverwrite,this->vm()) )
{
}

void
ConfigureKSP::run( KSP& ksp ) const
{
    // set ksp type : gmres,cg,preonly,...
    this->check( KSPSetType( ksp, M_type.c_str() ) );

    // init with petsc option if given and not interfaced
    if ( true )
        this->check( KSPSetFromOptions( ksp ) );

    if ( M_useConfigDefaultPetsc )
        return;

    // get ksp type
#if PETSC_VERSION_LESS_THAN(3,4,0) && !PETSC_VERSION_LESS_THAN(3,0,0)
    const KSPType ksp_type;
#else
    KSPType ksp_type;
#endif
    this->check( KSPGetType ( ksp, &ksp_type ) );

    // configure ksp from type
    if ( std::string((char*)ksp_type) == std::string( ( char* )KSPPREONLY ) )
    {
        this->check( KSPSetInitialGuessNonzero ( ksp, PETSC_FALSE ) );
    }
    else if ( std::string((char*)ksp_type) == std::string( ( char* )KSPGMRES ) )
    {
        this->check( KSPGMRESSetRestart( ksp, M_nRestartGMRES ) );
    }

    // Norm that is passed in the Krylov convergence test routines
    //this->check( KSPSetNormType( ksp, KSP_NORM_DEFAULT /*KSP_NORM_PRECONDITIONED*/ ) );

    // set ksp tolerance
    this->check( KSPSetTolerances( ksp,M_rtol,PETSC_DEFAULT,PETSC_DEFAULT,M_maxit ) );
    // monitor
    if ( M_showMonitor )
    {
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,7,0)
        this->check( KSPMonitorSet( ksp,__feel_petsc_prec_ksp_monitor,(void*) this,PETSC_NULL ) );
#else
        this->check( KSPMonitorSet( ksp,KSPMonitorDefault,PETSC_NULL,PETSC_NULL ) );
#endif
    }

    // constant null space
    if ( M_constantNullSpace )
    {
        MatNullSpace nullsp;
        this->check( MatNullSpaceCreate( PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nullsp ) );
#if PETSC_VERSION_LESS_THAN( 3,5,4 )
        this->check( KSPSetNullSpace( ksp, nullsp ) );
#else
        Mat A;
        this->check( KSPGetOperators( ksp, &A, NULL ) );
        this->check( MatSetNullSpace( A, nullsp ) );
#endif
        PETSc::MatNullSpaceDestroy( nullsp );
    }

}



/**
 * ConfigurePCLU
 */
ConfigurePCLU::ConfigurePCLU( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                              std::string const& sub, std::string const& prefix,
                              std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( precFeel,worldComm,sub,prefix,prefixOverwrite, getOptionsDescLU(prefix,sub,prefixOverwrite) ),
    M_matSolverPackage( ""),//getOption<std::string>("pc-factor-mat-solver-package-type",prefix,sub,prefixOverwrite) ),
    M_mumpsParameters( 33, std::make_pair(false,-1) )
{
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3, 9, 0 )    
    MatSolverType ptype;
    this->check( PCFactorGetMatSolverType(pc, &ptype ) );
#else
    const MatSolverPackage ptype;
    this->check( PCFactorGetMatSolverPackage(pc, &ptype ) );
#endif
    M_matSolverPackage = std::string( ( char* ) ptype );

    if ( M_matSolverPackage == "mumps" )
    {
#if defined(PETSC_HAVE_MUMPS)
        for ( int icntl=1 ; icntl<= M_mumpsParameters.size() ; ++icntl )
        {
            std::string mumpsOption = (boost::format("pc-factor-mumps.icntl-%1%")%icntl ).str();
            auto mumpsOptionAsked = getOptionIfAvalaible<int>(mumpsOption,prefix,sub,prefixOverwrite,this->vm());
            if ( mumpsOptionAsked.first )
                M_mumpsParameters[icntl-1] = mumpsOptionAsked;
        }
#endif
    }
    VLOG(2) << "ConfigurePC : LU\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->matSolverPackage : " << M_matSolverPackage << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}
void
ConfigurePCLU::run( PC& pc )
{
    // set factor package
    //this->check( PCFactorSetMatSolverPackage( pc, M_matSolverPackage.c_str() ) );

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
    // allow to tune the factorisation package
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,9,0 )
    this->check( PCFactorSetUpMatSolverType(pc) );
    
#else
    this->check( PCFactorSetUpMatSolverPackage(pc) );
#endif

    // configure mumps
    if ( M_matSolverPackage == "mumps" )
    {
#if defined(PETSC_HAVE_MUMPS)
        Mat F;
        this->check( PCFactorGetMatrix(pc,&F) );
        for ( int icntl=1 ; icntl<= M_mumpsParameters.size() ; ++icntl )
        {
            if ( M_mumpsParameters[icntl-1].first )
            {
                PetscInt ival = M_mumpsParameters[icntl-1].second;
                this->check( MatMumpsSetIcntl(F,icntl,ival) );
            }
        }
#else
        CHECK( false ) << "mumps not installed with PETSc";
#endif
    }
#endif

}

/**
 * ConfigurePCILU
 */
ConfigurePCILU::ConfigurePCILU( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                std::string const& sub, std::string const& prefix,
                                std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix,prefixOverwrite, getOptionsDescILU(prefix,sub,prefixOverwrite) ),
    M_levels( getOption<int>("pc-factor-levels",prefix,sub,prefixOverwrite,this->vm() ) ),
    M_fill( getOption<double>("pc-factor-fill",prefix,sub,prefixOverwrite,this->vm() ) )
{
    VLOG(2) << "ConfigurePC : ILU\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->levels : " << M_levels << "\n"
            << "  |->fill : " << M_fill << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}
void
ConfigurePCILU::run( PC& pc )
{
    // do we need to set the mat solver package for ilu ?
    //PetscPCFactorSetMatSolverPackage( pc, "petsc" );
    this->check( PCFactorSetLevels( pc, M_levels ) );
    this->check( PCFactorSetFill( pc, M_fill ) );
}
 
/**
 * ConfigurePCHYPRE_EUCLID
 */
#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_LESS_THAN( 3,6,0 )
ConfigurePCHYPRE_EUCLID::ConfigurePCHYPRE_EUCLID( PC& pc, PreconditionerPetsc<double> * precFeel,
                                                  WorldComm const& worldComm, std::string const& sub, std::string const& prefix,
                                                  std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix,prefixOverwrite,getOptionsDescILU(prefix,sub,prefixOverwrite) ),
    M_levels( option(_name="pc-factor-levels",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() )
{
    VLOG(2) << "ConfigurePC : HYPRE_EUCLID(ILU)\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->levels : " << M_levels << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}
void
ConfigurePCHYPRE_EUCLID::run( PC& pc )
{
#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
    this->check( PetscImpl::PCHYPRE_EUCLIDSetLevels( pc, M_levels ) );
#endif
}
#endif
/**
 * ConfigurePCHYPRE_BOOMERAMG
 */
ConfigurePCHYPRE_BOOMERAMG::ConfigurePCHYPRE_BOOMERAMG( PC& pc, PreconditionerPetsc<double> * precFeel,
                                                  WorldComm const& worldComm, std::string const& sub, std::string const& prefix,
                                                  std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix,prefixOverwrite,getOptionsDescBOOMERAMG(prefix,sub,prefixOverwrite) ),
    M_max_iter( option(_name="pc-hypre-boomeramg-max-iter",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_tol( option(_name="pc-hypre-boomeramg-tol",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<double>() ),
    M_cycle_type( option(_name="pc-hypre-boomeramg-cycle-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>()+1 ),
    M_max_levels( option(_name="pc-hypre-boomeramg-max-levels",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_coarsen_type( option(_name="pc-hypre-boomeramg-coarsen-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_strong_threshold( option(_name="pc-hypre-boomeramg-strong-threshold",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<double>() ),
    M_agg_nl( option(_name="pc-hypre-boomeramg-agg-nl",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_relax_type_all( option(_name="pc-hypre-boomeramg-relax-type-all",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_interp_type( option(_name="pc-hypre-boomeramg-interp-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() )
{
    VLOG(2) << "ConfigurePC : HYPRE_BOOMERAMG\n"
            << "  |->prefix     : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->max_iter   : " << M_max_iter << "\n"
            << "  |->tol        : " << M_tol << "\n"
            << "  |->cycle_type : " << M_cycle_type << "\n"
            << "  |->max_levels : " << M_max_levels << "\n"
            << "  |->coarsen_type : " << M_coarsen_type << "\n"
            << "  |->strong_threshold : " << M_strong_threshold << "\n"
            << "  |->agg nl : " << M_agg_nl << "\n"
            << "  |->relax_type_all : " << M_relax_type_all << "\n"
            << "  |->interp_type : " << M_interp_type << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}
void
ConfigurePCHYPRE_BOOMERAMG::run( PC& pc )
{
#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
    this->check( PetscImpl::PCHYPRE_BOOMERAMGSetMaxIter( pc, M_max_iter ) );
    this->check( PetscImpl::PCHYPRE_BOOMERAMGSetTol( pc, M_tol ) );
    this->check( PetscImpl::PCHYPRE_BOOMERAMGSetCycleType( pc, M_cycle_type ) );
    this->check( PetscImpl::PCHYPRE_BOOMERAMGSetMaxLevels( pc, M_max_levels ) );
    this->check( PetscImpl::PCHYPRE_BOOMERAMGSetCoarsenType( pc, M_coarsen_type ) );
    this->check( PetscImpl::PCHYPRE_BOOMERAMGSetStrongThreshold( pc, M_strong_threshold ) );
    this->check( PetscImpl::PCHYPRE_BOOMERAMGSetAggNumLevels( pc, M_agg_nl ) );
    this->check( PetscImpl::PCHYPRE_BOOMERAMGSetRelaxTypeAll( pc, M_relax_type_all ) );
    this->check( PetscImpl::PCHYPRE_BOOMERAMGSetInterpolationType( pc, M_interp_type ) );
#if 0
    const char     *petscPrefix;
    this->check( PCGetOptionsPrefix(pc,&petscPrefix ) );

    if( petscPrefix != NULL)
    {
      std::cout << petscPrefix << std::endl;
      std::cout << (boost::format("-%1%pc_hypre_boomeramg_max_iter"        )%std::string(petscPrefix)).str().c_str() << std::endl;
      this->check( PetscOptionsSetValue((boost::format("-%1%pc_hypre_boomeramg_max_iter"        )%std::string(petscPrefix)).str().c_str(), boost::lexical_cast<std::string>(M_max_iter        ).c_str()) );
      this->check( PetscOptionsSetValue((boost::format("-%1%pc_hypre_boomeramg_tol"             )%std::string(petscPrefix)).str().c_str(), boost::lexical_cast<std::string>(M_tol             ).c_str()) );
      this->check( PetscOptionsSetValue((boost::format("-%1%pc_hypre_boomeramg_cycle_type"      )%std::string(petscPrefix)).str().c_str(), boost::lexical_cast<std::string>(M_cycle_type      ).c_str()) );
      this->check( PetscOptionsSetValue((boost::format("-%1%pc_hypre_boomeramg_max_levels"      )%std::string(petscPrefix)).str().c_str(), boost::lexical_cast<std::string>(M_max_levels      ).c_str()) );
      this->check( PetscOptionsSetValue((boost::format("-%1%pc_hypre_boomeramg_coarsen_type"    )%std::string(petscPrefix)).str().c_str(), boost::lexical_cast<std::string>(M_coarsen_type    ).c_str()) );
      this->check( PetscOptionsSetValue((boost::format("-%1%pc_hypre_boomeramg_strong_threshold")%std::string(petscPrefix)).str().c_str(), boost::lexical_cast<std::string>(M_strong_threshold).c_str()) );
      this->check( PetscOptionsSetValue((boost::format("-%1%pc_hypre_boomeramg_agg_nl"          )%std::string(petscPrefix)).str().c_str(), boost::lexical_cast<std::string>(M_agg_nl          ).c_str()) );
      this->check( PetscOptionsSetValue((boost::format("-%1%pc_hypre_boomeramg_relax_type_all"  )%std::string(petscPrefix)).str().c_str(), boost::lexical_cast<std::string>(M_relax_type_all  ).c_str()) );
      this->check( PetscOptionsSetValue((boost::format("-%1%pc_hypre_boomeramg_interp_type"     )%std::string(petscPrefix)).str().c_str(), boost::lexical_cast<std::string>(M_interp_type     ).c_str()) );
    }
    else
    {
      std::cout << "OCucou" << std::endl;
      this->check( PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter"        , boost::lexical_cast<std::string>(M_max_iter        ).c_str()) );
      this->check( PetscOptionsSetValue("-pc_hypre_boomeramg_tol"             , boost::lexical_cast<std::string>(M_tol             ).c_str()) );
      this->check( PetscOptionsSetValue("-pc_hypre_boomeramg_cycle_type"      , boost::lexical_cast<std::string>(M_cycle_type      ).c_str()) );
      this->check( PetscOptionsSetValue("-pc_hypre_boomeramg_max_levels"      , boost::lexical_cast<std::string>(M_max_levels      ).c_str()) );
      this->check( PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type"    , boost::lexical_cast<std::string>(M_coarsen_type    ).c_str()) );
      this->check( PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", boost::lexical_cast<std::string>(M_strong_threshold).c_str()) );
      this->check( PetscOptionsSetValue("-pc_hypre_boomeramg_agg_nl"          , boost::lexical_cast<std::string>(M_agg_nl          ).c_str()) );
      this->check( PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_all"  , boost::lexical_cast<std::string>(M_relax_type_all  ).c_str()) );
      this->check( PetscOptionsSetValue("-pc_hypre_boomeramg_interp_type"     , boost::lexical_cast<std::string>(M_interp_type     ).c_str()) );
    }
#endif
#endif
}

/**
 * ConfigurePCHYPRE_AMS
 */
ConfigurePCHYPRE_AMS::ConfigurePCHYPRE_AMS( PC& pc, PreconditionerPetsc<double> * precFeel,
                                                  WorldComm const& worldComm, std::string const& sub, std::string const& prefix,
                                                  std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix,prefixOverwrite,getOptionsDescAMS(prefix,sub,prefixOverwrite) ),
    M_print_level(option(_name="pc-hypre-ams-print-level",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_max_iter(option(_name="pc-hypre-ams-max-iter",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_cycle_type(option(_name="pc-hypre-ams-cycle-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_tol(option(_name="pc-hypre-ams-tol",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<double>() ),
    M_relax_type(option(_name="pc-hypre-ams-relax-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_relax_times(option(_name="pc-hypre-ams-relax-times",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_relax_weight(option(_name="pc-hypre-ams-relax-weight",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<double>() ),
    M_omega(option(_name="pc-hypre-ams-omega",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<double>() )
    #if 0
    M_amg_alpha_theta(_name="pc-hypre-ams-amg-alpha-theta",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<double>())
    M_amg_alpha_options(_name="pc-hypre-ams-amg-alpha-options",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<double[5]>())
    M_amg_alpha_theta(_name="pc-hypre-ams-amg-beta-theta",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<double>())
    M_amg_alpha_options(_name="pc-hypre-ams-amg-beta-options",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<double[5]>())
    #endif
{
    VLOG(2) << "ConfigurePC : HYPRE_AMS\n"
            << "  |-> prefix       : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |-> print_level  : " << M_print_level << "\n"
            << "  |-> max_iter     : " << M_max_iter << "\n"
            << "  |-> cycle_type   : " << M_cycle_type << "\n"
            << "  |-> tolerance    : " << M_tol << "\n"
            << "  |-> relax_type   : " << M_relax_type << "\n"
            << "  |-> relax_times  : " << M_relax_times << "\n"
            << "  |-> relax_weight : " << M_relax_weight << "\n"
            << "  |-> omega        : " << M_omega;
    #if 0
    << "  |-> amg-alpha-theta        : " << M_amg_alpha_theta << "\n";
    << "  |-> amg-alpha-options      : " << M_amg_alpha_options << "\n";
    << "  |-> amg-beta-theta         : " << M_amg_beta_theta << "\n";
    << "  |-> amg-beta-options       : " << M_amg_beta_options << "\n";
    #endif
    google::FlushLogFiles(google::INFO);
    run( pc );
}
void
ConfigurePCHYPRE_AMS::run( PC& pc )
{
#if defined(PETSC_HAVE_HYPRE) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,6,0 )
    this->check( PetscImpl::PCHYPRE_AMSSetPrintLevel(pc, M_print_level));
    this->check( PetscImpl::PCHYPRE_AMSSetMaxIter(pc, M_max_iter));
    this->check( PetscImpl::PCHYPRE_AMSSetCycleType(pc, M_cycle_type));
    this->check( PetscImpl::PCHYPRE_AMSSetTol(pc, M_tol));
    this->check( PetscImpl::PCHYPRE_AMSSetSmoothingOptions(pc, M_relax_type, M_relax_times, M_relax_weight, M_omega));

#if 0
    this->check( PetscImpl::PCHYPRE_AMSSetAlphaAMGOptions(pc, M_amg_alpha_options[0],       /* AMG coarsen type */
                                                              M_amg_alpha_options[1],       /* AMG agg_levels */
                                                              M_amg_alpha_options[2],       /* AMG relax_type */
                                                              M_amg_alpha_theta,
                                                              M_amg_alpha_options[3],       /* AMG interp_type */
                                                              M_amg_alpha_options[4]);
    this->check( PetscImpl::PCHYPRE_AMSSetBetaAMGOptions(pc, M_amg_beta_options[0],       /* AMG coarsen type */
                                                             M_amg_beta_options[1],       /* AMG agg_levels */
                                                             M_amg_beta_options[2],       /* AMG relax_type */
                                                             M_amg_beta_theta,
                                                             M_amg_beta_options[3],       /* AMG interp_type */
                                                             M_amg_beta_options[4]);
#endif

    if (!pc->setupcalled)
    {
        // from hypre doc, this function should be called before HYPRE AMSSetup(), so must be init once
        if ( this->precFeel()->hasAuxiliarySparseMatrix("G") )
        {
            auto gMat = this->precFeel()->auxiliarySparseMatrix("G");
            CHECK(gMat) << "The pointer gMat is not initialized\n";
            MatrixPetsc<double> * gPetsc   = const_cast<MatrixPetsc<double> *>( dynamic_cast<MatrixPetsc<double> const*>( &(*gMat) ) );
            CHECK( gPetsc && gPetsc->mat() ) << "gPetsc->mat() is not initialized\n";
            this->check( PCHYPRESetDiscreteGradient( pc, gPetsc->mat() ) );
        }
        else
            std::cerr << "G for hypre AMS has not been provided\n";

        if ( this->precFeel()->hasAuxiliaryVector("Px") && this->precFeel()->hasAuxiliaryVector("Py") && this->precFeel()->hasAuxiliaryVector("Pz")  )
        {
            auto pxVec = this->precFeel()->auxiliaryVector("Px");
            auto pyVec = this->precFeel()->auxiliaryVector("Py");
            auto pzVec = this->precFeel()->auxiliaryVector("Pz");
            VectorPetsc<double> * pxPetsc   = const_cast<VectorPetsc<double> *>( dynamic_cast<VectorPetsc<double> const*>( &(*pxVec) ) );
            VectorPetsc<double> * pyPetsc   = const_cast<VectorPetsc<double> *>( dynamic_cast<VectorPetsc<double> const*>( &(*pyVec) ) );
            VectorPetsc<double> * pzPetsc   = const_cast<VectorPetsc<double> *>( dynamic_cast<VectorPetsc<double> const*>( &(*pzVec) ) );
            this->check( PCHYPRESetEdgeConstantVectors(pc, pxPetsc->vec(), pyPetsc->vec(), pzPetsc->vec()) );
        }
        else if ( this->precFeel()->hasAuxiliaryVector("X") && this->precFeel()->hasAuxiliaryVector("Y") && this->precFeel()->hasAuxiliaryVector("Z")  )
        {
            auto pxVec = this->precFeel()->auxiliaryVector("X");
            auto pyVec = this->precFeel()->auxiliaryVector("Y");
            auto pzVec = this->precFeel()->auxiliaryVector("Z");
            VectorPetsc<double> * pxPetsc   = const_cast<VectorPetsc<double> *>( dynamic_cast<VectorPetsc<double> const*>( &(*pxVec) ) );
            VectorPetsc<double> * pyPetsc   = const_cast<VectorPetsc<double> *>( dynamic_cast<VectorPetsc<double> const*>( &(*pyVec) ) );
            VectorPetsc<double> * pzPetsc   = const_cast<VectorPetsc<double> *>( dynamic_cast<VectorPetsc<double> const*>( &(*pzVec) ) );
            this->check( PetscImpl::PCHYPRE_AMSSetCoordinateVectors(pc, pxPetsc->vec(), pyPetsc->vec(), pzPetsc->vec()));
        }
        else
            std::cerr << "Nor (Px, Py, Pz), nor (X, Y, Z) has been provided\n";
    } // !pc->setupcalled


    if ( this->precFeel()->hasAuxiliarySparseMatrix("a_alpha") )
    {
        auto gMat = this->precFeel()->auxiliarySparseMatrix("a_alpha");
        MatrixPetsc<double> * gPetsc   = const_cast<MatrixPetsc<double> *>( dynamic_cast<MatrixPetsc<double> const*>( &(*gMat) ) );
        this->check( PetscImpl::PCHYPRE_AMSSetAlphaPoissonMatrix_HYPRE(pc, gPetsc->mat()));
    }
    if ( this->precFeel()->hasAuxiliarySparseMatrix("a_beta") )
    {
        auto gMat = this->precFeel()->auxiliarySparseMatrix("a_beta");
        if(!gMat)
        {
          this->check( PetscImpl::PCHYPRE_AMSSetBetaPoissonMatrix_HYPRE(pc, NULL));
        }
        else
        {
          MatrixPetsc<double> * gPetsc   = const_cast<MatrixPetsc<double> *>( dynamic_cast<MatrixPetsc<double> const*>( &(*gMat) ) );
          this->check( PetscImpl::PCHYPRE_AMSSetBetaPoissonMatrix_HYPRE(pc, gPetsc->mat()));
        }
    }
#endif
}


/**
 * ConfigurePCSOR
 */
ConfigurePCSOR::ConfigurePCSOR( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                std::string const& sub, std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix,prefixOverwrite,getOptionsDescSOR(prefix,sub,prefixOverwrite) ),
    M_type( getOption<std::string>("pc-sor-type",prefix,sub,prefixOverwrite,this->vm() ) ),
    M_omega( getOption<double>("pc-sor-omega",prefix,sub,prefixOverwrite,this->vm()) ),
    M_nIteration( getOption<int>("pc-sor-its",prefix,sub,prefixOverwrite,this->vm()) ),
    M_nLocalIteration( getOption<int>("pc-sor-lits",prefix,sub,prefixOverwrite,this->vm()) )
{
    VLOG(2) << "ConfigurePC : ILU\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->omega : " << M_omega << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}
void
ConfigurePCSOR::run( PC& pc )
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
ConfigurePCGASM::ConfigurePCGASM( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                  std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( precFeel, worldComm,"",prefix,prefixOverwrite, getOptionsDescGASM( prefix,prefixOverwrite ) ),
    M_type( getOption<std::string>("pc-gasm-type",prefix,"",prefixOverwrite, this->vm() ) ),
    M_overlap( getOption<int>("pc-gasm-overlap",prefix,"",prefixOverwrite, this->vm()) )
{
    VLOG(2) << "ConfigurePC : GASM\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->type : " << M_type  << "\n"
            << "  |->overlap : " << M_overlap << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}
void
ConfigurePCGASM::run( PC& pc )
{
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
    /**/ if ( M_type == "restrict" )    this->check( PCGASMSetType( pc, PC_GASM_RESTRICT ) );
    else if ( M_type == "basic" )       this->check( PCGASMSetType( pc, PC_GASM_BASIC ) );
    else if ( M_type == "interpolate" ) this->check( PCGASMSetType( pc, PC_GASM_INTERPOLATE ) );
    else if ( M_type == "none" )        this->check( PCGASMSetType( pc, PC_GASM_NONE ) );
    else                                CHECK( false ) << "invalid gasm type : " << M_type << "\n";
    this->check( PCGASMSetOverlap( pc, M_overlap ) );
#endif
    ConfigureSubPC( pc,this->precFeel(),this->worldComm().subWorldCommSeq(),this->prefix(),this->prefixOverwrite() );
}

/**
 * ConfigurePCGASM
 */
ConfigurePCASM::ConfigurePCASM( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( precFeel, worldComm,"",prefix,prefixOverwrite, getOptionsDescASM( prefix,prefixOverwrite ) ),
    M_type( getOption<std::string>("pc-asm-type",prefix,"",prefixOverwrite,this->vm() ) ),
    M_overlap( getOption<int>("pc-asm-overlap",prefix,"",prefixOverwrite,this->vm() ) )
{
    VLOG(2) << "ConfigurePC : ASM\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->type : " << M_type  << "\n"
            << "  |->overlap : " << M_overlap << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}
void
ConfigurePCASM::run( PC& pc )
{
    /**/ if ( M_type == "restrict" )    this->check( PCASMSetType( pc, PC_ASM_RESTRICT ) );
    else if ( M_type == "basic" )       this->check( PCASMSetType( pc, PC_ASM_BASIC ) );
    else if ( M_type == "interpolate" ) this->check( PCASMSetType( pc, PC_ASM_INTERPOLATE ) );
    else if ( M_type == "none" )        this->check( PCASMSetType( pc, PC_ASM_NONE ) );
    else                                CHECK( false ) << "invalid asm type : " << M_type << "\n";
    this->check( PCASMSetOverlap( pc, M_overlap ) );

    ConfigureSubPC( pc,this->precFeel(),this->worldComm().subWorldCommSeq(),this->prefix(),this->prefixOverwrite() );
}

/**
 * ConfigureSubPC
 */
ConfigureSubPC::ConfigureSubPC( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( precFeel, worldComm,"",prefix, prefixOverwrite, getOptionsDescPrecBase(prefix,"sub", prefixOverwrite ) ),
    M_subPCtype( getOption<std::string>("pc-type",prefix,"sub",prefixOverwrite,this->vm() ) ),
    M_subMatSolverPackage( getOption<std::string>("pc-factor-mat-solver-package-type",prefix,"sub",prefixOverwrite,this->vm() ) ),
    M_subPCview( getOption<bool>("pc-view",prefix,"sub",prefixOverwrite,this->vm() ) ),
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

    ConfigureKSP kspConf( this->precFeel(),this->worldComm(), "sub", this->prefix(), this->prefixOverwrite() );

    for ( int i=0; i<M_nBlock; ++i )
    {
        run( subksps[i],kspConf );
    }
}
void
ConfigureSubPC::run( KSP& ksp, ConfigureKSP const& kspConf )
{
    // configure coarse ksp
    //ConfigureKSP kspConf( ksp, this->worldComm(), "sub", this->prefix(), this->prefixOverwrite() );
    kspConf.run( ksp );
    this->check( KSPSetUp( ksp ) );

    PC subpc;
    // Get pointer to sub KSP object's PC
    this->check( KSPGetPC( ksp, &subpc ) );

    // configure sub-pc
    SetPCType( subpc, pcTypeConvertStrToEnum( M_subPCtype ),
               matSolverPackageConvertStrToEnum( M_subMatSolverPackage ),
               this->worldComm() );
    ConfigurePC( subpc, this->precFeel(), this->worldComm(), "sub", this->prefix(), this->prefixOverwrite(), this->vm() );
    this->check( PCSetUp( subpc ) );

    if ( kspConf.kspView() )
        this->check( KSPView( ksp, PETSC_VIEWER_STDOUT_SELF ) );
    else if ( M_subPCview )
        this->check( PCView( subpc, PETSC_VIEWER_STDOUT_SELF ) );
}

/**
 * ConfigurePCML
 */
ConfigurePCML::ConfigurePCML( PC& pc, PreconditionerPetsc<double> * precFeel,
                              WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix, getOptionsDescML(prefix,sub) ),
    M_mgType( option(_name="pc-mg-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<std::string>() ),
    M_nLevels( option(_name="pc-mg-levels",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_mlReuseInterp( option(_name="pc-ml-reuse-interpolation",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<bool>() ),
    M_mlKeepAggInfo( option(_name="pc-ml-keep-agg-info",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<bool>() ),
    M_mlReusable( option(_name="pc-ml-reusable",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<bool>() ),
    M_mlOldHierarchy( option(_name="pc-ml-old-hierarchy",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<bool>() ),
    M_prefixMGCoarse( (boost::format( "%1%%2%mg-coarse" ) %prefixvm( prefix,"" ) %std::string((sub.empty())?"":sub+"-")  ).str() ),
    M_coarsePCtype( option(_name="pc-type",_prefix=M_prefixMGCoarse,_vm=this->vm()).as<std::string>() ),
    M_coarsePCMatSolverPackage( option(_name="pc-factor-mat-solver-package-type",_prefix=M_prefixMGCoarse,_vm=this->vm()).as<std::string>() ),
    M_coarsePCview( option(_name="pc-view",_prefix=M_prefixMGCoarse,_vm=this->vm()).as<bool>() )
{
    VLOG(2) << "ConfigurePC : ML\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->mgType : " << M_mgType << "\n"
            << "  |->maxLevels : " << M_nLevels << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}
void
ConfigurePCML::run( PC& pc )
{
#if 0
    // Sets the number of levels to use with MG.
    // Must be called before any other MG routine
    this->check( PCMGSetLevels( pc, M_nLevels, PETSC_NULL) );
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
    this->check( PetscImpl::PCMLSetMaxNlevels( pc, M_nLevels ) );
    this->check( PetscImpl::PCMLSetReuseInterpolation( pc, static_cast<PetscBool>( M_mlReuseInterp ) ) );
    this->check( PetscImpl::PCMLSetKeepAggInfo( pc, static_cast<PetscBool>( M_mlKeepAggInfo ) ) );
    this->check( PetscImpl::PCMLSetReusable( pc, static_cast<PetscBool>( M_mlReusable ) ) );
    this->check( PetscImpl::PCMLSetOldHierarchy( pc, static_cast<PetscBool>( M_mlOldHierarchy ) ) );
#endif

    // setup ml pc
    this->check( PCSetUp( pc ) );
    // configure coarse solver
    configurePCMLCoarse( pc );
    // configure fine solvers
    ConfigurePCMGLevels( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix() );

    // must be called after setup ml pc because this one PCMGSetLevels and reset mg prec associated
    if ( M_mgType=="multiplicative" ) this->check( PCMGSetType( pc, PC_MG_MULTIPLICATIVE ) );
    else if ( M_mgType=="additive" ) this->check( PCMGSetType( pc, PC_MG_ADDITIVE ) );
    else if ( M_mgType=="full" ) this->check( PCMGSetType( pc, PC_MG_FULL ) );
    else if ( M_mgType=="kaskade" ) this->check( PCMGSetType( pc, PC_MG_KASKADE ) );
    else CHECK( false ) << "invalid mgType :" << M_mgType << "\n";

}
void
ConfigurePCML::configurePCMLCoarse( PC& pc )
{
    std::vector<std::string> prefixOverwrite;

    // get coarse-ksp
    KSP coarseksp;
    this->check( PCMGGetCoarseSolve( pc, &coarseksp) );

    // get coarse pc
    PC coarsepc;
    this->check( KSPGetPC( coarseksp, &coarsepc ) );

    // in order to setup our ksp config, call PCSetType (with != name) reset the prec
    if ( coarsepc->setupcalled )
        this->check( PCSetType(coarsepc, ( char* )PCNONE) );
    // configure coarse pc
    SetPCType( coarsepc, pcTypeConvertStrToEnum( M_coarsePCtype ),
               matSolverPackageConvertStrToEnum( M_coarsePCMatSolverPackage ),
               this->worldComm() );
    ConfigurePC coarsepcConf( this->precFeel(), /*coarsepc, is,*/ this->worldComm(), "", M_prefixMGCoarse, prefixOverwrite, this->vm() );
    coarsepcConf.setFactorShiftType( "inblocks" );
    coarsepcConf.run( coarsepc );
    // setup coarse pc
    this->check( PCSetUp( coarsepc ) );

    // configure coarse ksp
    ConfigureKSP kspConf( coarseksp, this->precFeel(), this->worldComm(), "", M_prefixMGCoarse, prefixOverwrite, "preonly",1e-5,50 );
    // setup coarse ksp
    this->check( KSPSetUp( coarseksp ) );

    PetscViewer viewer = (this->sub().empty())? PETSC_VIEWER_STDOUT_WORLD : PETSC_VIEWER_STDOUT_SELF;
    if ( kspConf.kspView() )
        this->check( KSPView( coarseksp, viewer ) );
    else if ( M_coarsePCview )
        this->check( PCView( coarsepc, viewer ) );
}

/**
 * ConfigurePCGAMG
 */
ConfigurePCGAMG::ConfigurePCGAMG( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                  std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix, getOptionsDescGAMG(prefix,sub) ),
    M_mgType( option(_name="pc-mg-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<std::string>() ),
    M_gamgType( option(_name="pc-gamg-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<std::string>() ),
    M_nLevels( option(_name="pc-mg-levels",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_procEqLim( option(_name="pc-gamg-proc-eq-lim",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_coarseEqLim(option(_name="pc-gamg-coarse-eq-lim",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_threshold( option(_name="pc-gamg-threshold",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<double>() ),
    M_setSymGraph( option(_name="pc-gamg-set-sym-graph",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<bool>() ),
    M_reuseInterpolation( option(_name="pc-gamg-reuse-interpolation",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<bool>() ),
    M_gamgVerbose( option(_name="pc-gamg-verbose",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_nSmooths( option(_name="pc-gamg-nsmooths",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<int>() ),
    M_prefixMGCoarse( (boost::format( "%1%%2%mg-coarse" ) %prefixvm( prefix,"" ) %std::string((sub.empty())?"":sub+"-")  ).str() ),
    M_coarsePCtype( option(_name="pc-type",_prefix=M_prefixMGCoarse,_vm=this->vm()).as<std::string>() ),
    M_coarsePCMatSolverPackage( option(_name="pc-factor-mat-solver-package-type",_prefix=M_prefixMGCoarse,_vm=this->vm()).as<std::string>() ),
    M_coarsePCview( option(_name="pc-view",_prefix=M_prefixMGCoarse,_vm=this->vm()).as<bool>() )
{
    VLOG(2) << "ConfigurePC : GAMG\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->mgType : " << M_mgType << "\n"
            << "  |->maxLevels : " << M_nLevels << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}

void
ConfigurePCGAMG::run( PC& pc )
{
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
    if ( !pc->setupcalled )
    {
        // set type of multigrid (agg only supported)
        this->check( PCGAMGSetType( pc, M_gamgType.c_str() ) );
        // Set for asymmetric matrices
        //this->check( PetscOptionsSetValue("-pc_gamg_sym_graph", boost::lexical_cast<std::string>(M_setSymGraph).c_str()) );
        //this->check( PetscOptionsSetValue("-fieldsplit_0_pc_gamg_sym_graph", boost::lexical_cast<std::string>(M_setSymGraph).c_str()) );
        const char     *petscPrefix;
        this->check( PCGetOptionsPrefix(pc,&petscPrefix ) );
        if ( petscPrefix != NULL )
        {
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,7,0 )
            this->check( PetscOptionsSetValue( NULL, (boost::format("-%1%pc_gamg_sym_graph")%std::string(petscPrefix)).str().c_str(),
                                               boost::lexical_cast<std::string>(M_setSymGraph).c_str()) );
            this->check( PetscOptionsSetValue( NULL, (boost::format("-%1%pc_gamg_verbose")%std::string(petscPrefix)).str().c_str(),
                                               boost::lexical_cast<std::string>(M_gamgVerbose).c_str()) );
            this->check( PetscOptionsSetValue( NULL, (boost::format("-%1%pc_gamg_agg_nsmooths")%std::string(petscPrefix)).str().c_str(),
                                               boost::lexical_cast<std::string>(M_nSmooths).c_str()) );
#else
            this->check( PetscOptionsSetValue( (boost::format("-%1%pc_gamg_sym_graph")%std::string(petscPrefix)).str().c_str(),
                                               boost::lexical_cast<std::string>(M_setSymGraph).c_str()) );
            this->check( PetscOptionsSetValue( (boost::format("-%1%pc_gamg_verbose")%std::string(petscPrefix)).str().c_str(),
                                               boost::lexical_cast<std::string>(M_gamgVerbose).c_str()) );
            this->check( PetscOptionsSetValue( (boost::format("-%1%pc_gamg_agg_nsmooths")%std::string(petscPrefix)).str().c_str(),
                                               boost::lexical_cast<std::string>(M_nSmooths).c_str()) );
#endif
        }
        else
        {
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,7,0 )
            this->check( PetscOptionsSetValue(NULL, "-pc_gamg_sym_graph", boost::lexical_cast<std::string>(M_setSymGraph).c_str()) );
            this->check( PetscOptionsSetValue(NULL, "-pc_gamg_verbose", boost::lexical_cast<std::string>(M_gamgVerbose).c_str()) );
            this->check( PetscOptionsSetValue(NULL, "-pc_gamg_agg_nsmooths", boost::lexical_cast<std::string>(M_nSmooths).c_str()) );
#else
            this->check( PetscOptionsSetValue("-pc_gamg_sym_graph", boost::lexical_cast<std::string>(M_setSymGraph).c_str()) );
            this->check( PetscOptionsSetValue("-pc_gamg_verbose", boost::lexical_cast<std::string>(M_gamgVerbose).c_str()) );
            this->check( PetscOptionsSetValue("-pc_gamg_agg_nsmooths", boost::lexical_cast<std::string>(M_nSmooths).c_str()) );
#endif
        }

        // PCSetFromOptions is called here because PCGAMGSetType destroy all unless the type_name
        this->check( PCSetFromOptions( pc ) );

        //
        this->check( PCGAMGSetNlevels( pc,M_nLevels ) );
        // Set number of equations to aim for on coarse grids via processor reduction
        this->check( PCGAMGSetProcEqLim( pc, M_procEqLim ) );
        // Set max number of equations on coarse grids
        this->check( PCGAMGSetCoarseEqLim( pc, M_coarseEqLim ) );
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,8,0 )
        this->check( PCGAMGSetThresholdScale( pc, M_threshold ) );
#else
        // Relative threshold to use for dropping edges in aggregation graph
        this->check( PCGAMGSetThreshold( pc, M_threshold ) );
#endif
        // not works!!(seems to be missing PetscObjectComposeFunction with this function)
        //this->check( PCGAMGSetSymGraph( pc, ( M_setSymGraph )?PETSC_TRUE : PETSC_FALSE ) );
        // Reuse prolongation operator
#if PETSC_VERSION_LESS_THAN( 3,5,5 )
        this->check( PCGAMGSetReuseProl( pc, ( M_reuseInterpolation )?PETSC_TRUE : PETSC_FALSE ) );
#else
        this->check( PCGAMGSetReuseInterpolation( pc, ( M_reuseInterpolation )?PETSC_TRUE : PETSC_FALSE ) );
#endif
        // not work also
        //this->check( PCGAMGSetNSmooths( pc, 30 ) );

    }

    // setup sub-pc
    this->check( PCSetUp( pc ) );
    //this->check( PCView( pc, PETSC_VIEWER_STDOUT_WORLD ) );

    // configure coarse pc
    configurePCGAMGCoarse( pc );
    // configure level pc
    ConfigurePCMGLevels( pc, this->precFeel(), this->worldComm(), this->sub(), this->prefix() );

    // must be called after setup gamg pc because this one call PCMGSetLevels and reset mg prec associated
    if ( M_mgType=="multiplicative" ) this->check( PCMGSetType( pc, PC_MG_MULTIPLICATIVE ) );
    else if ( M_mgType=="additive" ) this->check( PCMGSetType( pc, PC_MG_ADDITIVE ) );
    else if ( M_mgType=="full" ) this->check( PCMGSetType( pc, PC_MG_FULL ) );
    else if ( M_mgType=="kaskade" ) this->check( PCMGSetType( pc, PC_MG_KASKADE ) );
    else CHECK( false ) << "invalid mgType :" << M_mgType << "\n";
#else // petsc version >= 3.3
    CHECK( false ) << "gamg supported only from petsc 3.3";
#endif
}

void
ConfigurePCGAMG::configurePCGAMGCoarse( PC& pc )
{
    std::vector<std::string> prefixOverwrite;

    // get coarse-ksp
    KSP coarseksp;
    this->check( PCMGGetCoarseSolve( pc, &coarseksp) );

    // get coarse pc
    PC coarsepc;
    this->check( KSPGetPC( coarseksp, &coarsepc ) );

    // in order to setup our ksp config, call PCSetType (with != name) reset the prec
    if ( coarsepc->setupcalled )
        this->check( PCSetType(coarsepc, ( char* )PCNONE) );
    // configure coarse pc
    SetPCType( coarsepc, pcTypeConvertStrToEnum( M_coarsePCtype ),
               matSolverPackageConvertStrToEnum( M_coarsePCMatSolverPackage ),
               this->worldComm() );
    ConfigurePC coarsepcConf( this->precFeel(), /*coarsepc, is,*/ this->worldComm(), "", M_prefixMGCoarse, prefixOverwrite, this->vm() );
    coarsepcConf.setFactorShiftType( "inblocks" );
    coarsepcConf.run( coarsepc );
    // setup coarse pc
    this->check( PCSetUp( coarsepc ) );

    // configure coarse ksp
    ConfigureKSP kspConf( coarseksp, this->precFeel(), this->worldComm(), "", M_prefixMGCoarse, prefixOverwrite, "preonly",1e-5,50 );
    // setup coarse ksp
    this->check( KSPSetUp( coarseksp ) );


    PetscViewer viewer = (this->sub().empty())? PETSC_VIEWER_STDOUT_WORLD : PETSC_VIEWER_STDOUT_SELF;
    if ( kspConf.kspView() )
        this->check( KSPView( coarseksp, viewer ) );
    else if ( M_coarsePCview )
        this->check( PCView( coarsepc, viewer ) );
}

/**
 * ConfigurePCMGLevels
 */
ConfigurePCMGLevels::ConfigurePCMGLevels( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                          std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( precFeel,worldComm,sub,prefix )//, getOptionsDescMultiGridLevels( prefix,sub) )
{
    //this->check( PCSetUp( pc ) );
    this->check( PCMGGetLevels( pc, &M_nLevels) );
    //std::cout << "M_nLevels " << M_nLevels << std::endl;
    this->initVariableMap( getOptionsDescMultiGridLevels( M_nLevels-1,prefix,sub) );


    M_prefixMGLevels.resize( M_nLevels-1 );
    M_mgLevelsPCtype.resize( M_nLevels-1 );
    M_mgLevelsPCview.resize( M_nLevels-1 );
    //M_mgLevelsKSPview.resize( M_nLevels-1 );
    M_mgLevelsMatSolverPackage.resize( M_nLevels-1 );
    std::string mgctx = (sub.empty())? "mg-" : sub+"-mg-";

    // get generic option for all levels
    std::string prefixAllLevel = ( boost::format( "%1%%2%levels" ) %prefixvm( this->prefix(),"" ) %mgctx ).str();
    //M_mgLevelsKSPview[0] = option(_name="ksp-view",_prefix=prefixAllLevel).as<bool>();
    M_mgLevelsPCtype[0] = option(_name="pc-type",_prefix=prefixAllLevel,_vm=this->vm()).as<std::string>();
    M_mgLevelsPCview[0] = option(_name="pc-view",_prefix=prefixAllLevel,_vm=this->vm()).as<bool>();
    M_mgLevelsMatSolverPackage[0] = option(_name="pc-factor-mat-solver-package-type",_prefix=prefixAllLevel,_vm=this->vm()).as<std::string>();
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
    for ( int level=1; level<M_nLevels/*std::min(M_nLevels,6)*/; ++level )
    {
        std::string prefixCurLevel = ( boost::format( "%1%%2%levels%3%" ) %prefixvm( this->prefix(),"" ) %mgctx %level ).str();
        M_prefixMGLevels[level-1] = prefixCurLevel;
        //if ( this->vm().count( prefixvm(prefixCurLevel,"ksp-view") ) )
        //    M_mgLevelsKSPview[level-1] = option(_name="ksp-view",_prefix=prefixCurLevel,_vm=this->vm()).as<bool>();
        if ( this->vm().count( prefixvm(prefixCurLevel,"pc-type") ) )
            M_mgLevelsPCtype[level-1] = option(_name="pc-type",_prefix=prefixCurLevel,_vm=this->vm()).as<std::string>();
        if ( this->vm().count( prefixvm(prefixCurLevel,"pc-view") ) )
            M_mgLevelsPCview[level-1] = option(_name="pc-view",_prefix=prefixCurLevel,_vm=this->vm()).as<bool>();
        if ( this->vm().count( prefixvm(prefixCurLevel,"pc-factor-mat-solver-package-type") ) )
            M_mgLevelsMatSolverPackage[level-1] = option(_name="pc-factor-mat-solver-package-type",_prefix=prefixCurLevel,_vm=this->vm()).as<std::string>();
    }
    // overwrite options for fine level
    std::string prefixFineLevel = ( boost::format( "%1%%2%fine-level" ) %prefixvm( this->prefix(),"" ) %mgctx ).str();
    //M_prefixMGLevels[M_nLevel-2] = prefixFineLevel;
    //if ( this->vm().count( prefixvm(prefixFineLevel,"ksp-view") ) )
    //M_mgLevelsKSPview[M_nLevels-2] = option(_name="ksp-view",_prefix=prefixFineLevel,_vm=this->vm()).as<bool>();
    if ( /*Environment::*/this->vm().count( prefixvm(prefixFineLevel,"pc-type") ) )
        M_mgLevelsPCtype[M_nLevels-2] = option(_name="pc-type",_prefix=prefixFineLevel,_vm=this->vm()).as<std::string>();
    if ( /*Environment::*/this->vm().count( prefixvm(prefixFineLevel,"pc-view") ) )
        M_mgLevelsPCview[M_nLevels-2] = option(_name="pc-view",_prefix=prefixFineLevel,_vm=this->vm()).as<bool>();
    if ( /*Environment::*/this->vm().count( prefixvm(prefixFineLevel,"pc-factor-mat-solver-package-type") ) )
        M_mgLevelsMatSolverPackage[M_nLevels-2] = option(_name="pc-factor-mat-solver-package-type",_prefix=prefixFineLevel,_vm=this->vm()).as<std::string>();

    // configure each levels
    for ( int level=1; level<M_nLevels; ++level )
        run( pc,level );
}

void
ConfigurePCMGLevels::run( PC& pc, int level )
{
    std::string prefixCurLevel = M_prefixMGLevels[level-1];
    std::string mgctx = (this->sub().empty())? "mg-" : this->sub()+"-mg-";
    std::string prefixAllLevel = ( boost::format( "%1%%2%levels" ) %prefixvm( this->prefix(),"" ) %mgctx ).str();
    std::string prefixFineLevel = ( boost::format( "%1%%2%fine-level" ) %prefixvm( this->prefix(),"" ) %mgctx ).str();

    std::vector<std::string> prefixLevelOverwrite;
    //if ( level<6 )
    prefixLevelOverwrite.push_back( prefixCurLevel );
    if ( level == (M_nLevels-1) )
        prefixLevelOverwrite.push_back( prefixFineLevel );
    //-------------------------------------------------------------------//
    KSP levelksp;
    // get ksp
    this->check( PCMGGetSmoother( pc, level, &levelksp ) );
    // get level pc
    PC levelpc;
    this->check( KSPGetPC( levelksp, &levelpc ) );

    //-------------------------------------------------------------------//
    //if ( levelpc->setupcalled )
    //    this->check( PCSetType(levelpc, ( char* )PCNONE) );
    // configure level pc
    SetPCType( levelpc, pcTypeConvertStrToEnum( M_mgLevelsPCtype[level-1] ),
               matSolverPackageConvertStrToEnum( M_mgLevelsMatSolverPackage[level-1] ),
               this->worldComm() );
    ConfigurePC( levelpc, this->precFeel(), this->worldComm(), "", prefixAllLevel , prefixLevelOverwrite, this->vm() );
    // setup level pc
    this->check( PCSetUp( levelpc ) );

    //-------------------------------------------------------------------//
    // configure ksp
    ConfigureKSP kspConf( levelksp,this->precFeel(), this->worldComm(), "", prefixAllLevel, prefixLevelOverwrite,"richardson", 1e-5, 2 );
#if 0
    // warning : use KSP_NORM_PRECONDITIONED and force convergence
    this->check( KSPSetNormType( levelksp, KSP_NORM_PRECONDITIONED ) );
    void *cctx;
    this->check( KSPDefaultConvergedCreate(&cctx) );
    this->check( KSPSetConvergenceTest( levelksp, KSPDefaultConverged, cctx, PETSC_NULL ) );
#endif
    // setup coarse ksp
    this->check( KSPSetUp( levelksp ) );

    PetscViewer viewer = (this->sub().empty())? PETSC_VIEWER_STDOUT_WORLD : PETSC_VIEWER_STDOUT_SELF;
    if ( kspConf.kspView() )
        this->check( KSPView( levelksp, viewer ) );
    else if ( M_mgLevelsPCview[level-1] )
        this->check( PCView( levelpc, viewer ) );
}

/**
 * ConfigurePCFieldSplit
 */
ConfigurePCFieldSplit::ConfigurePCFieldSplit( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                              std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix, getOptionsDescFieldSplit( prefix,sub ) ),
    M_type( option(_name="fieldsplit-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<std::string>() ),
    M_schurFactType( option(_name="fieldsplit-schur-fact-type",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<std::string>() ),
    M_schurPrecond( option(_name="fieldsplit-schur-precondition",_prefix=prefix,_sub=sub,_worldcomm=worldComm,_vm=this->vm()).as<std::string>() )
{
    VLOG(2) << "ConfigurePC : FieldSplit\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->type : " << M_type << "\n";
    run( pc );
}
void
ConfigurePCFieldSplit::run( PC& pc )
{
    /*const PetscInt ufields[] = {0,2},pfields[] = {1};
      this->check( PCFieldSplitSetFields( pc , NULL, 2, ufields,ufields) );
      this->check( PCFieldSplitSetFields( pc , NULL, 1, pfields,pfields) );*/

    PCCompositeType theFieldSplitType = PC_COMPOSITE_SCHUR;
    /**/ if ( M_type == "schur" ) theFieldSplitType = PC_COMPOSITE_SCHUR;
    else if ( M_type == "additive" ) theFieldSplitType = PC_COMPOSITE_ADDITIVE;
    else if ( M_type == "multiplicative" ) theFieldSplitType = PC_COMPOSITE_MULTIPLICATIVE;
    else if ( M_type == "symmetric-multiplicative" ) theFieldSplitType = PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE;
    else if ( M_type == "special" ) theFieldSplitType = PC_COMPOSITE_SPECIAL;
    this->check( PCFieldSplitSetType( pc, theFieldSplitType ) );

#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,3,0 )
    if ( M_type == "schur" )
    {
        PCFieldSplitSchurFactType theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_FULL;
        /**/ if ( M_schurFactType == "diag")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_DIAG;
        else if ( M_schurFactType == "lower")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_LOWER;
        else if ( M_schurFactType == "upper")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_UPPER;
        else if ( M_schurFactType == "full")  theSchurFactType = PC_FIELDSPLIT_SCHUR_FACT_FULL;
        this->check( PCFieldSplitSetSchurFactType( pc,theSchurFactType ) );

        PCFieldSplitSchurPreType theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_SELF;
        /**/ if ( M_schurPrecond == "self")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_SELF;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,6,0 )
        else if ( M_schurPrecond == "selfp")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_SELFP;
#endif
        else if ( M_schurPrecond == "user")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_USER;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,4,0 )
        else if ( M_schurPrecond == "a11")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_A11;
#else
        else if ( M_schurPrecond == "a11")  theSchurPrecond = PC_FIELDSPLIT_SCHUR_PRE_DIAG;
#endif

        Mat schurMatPrecond = NULL;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,4,0 )
        if ( M_schurPrecond == "user" )
        {
            this->check( PCSetUp( pc ) );

            Mat schur = NULL, A, B, C, D;
            this->check( PetscImpl::PCFieldSplit_GetMatSchurComplement( pc, schur ) );
#if PETSC_VERSION_LESS_THAN(3,5,0)
            this->check( MatSchurComplementGetSubmatrices( schur,&A,NULL,&B,&C,&D ) );
#else
            this->check( MatSchurComplementGetSubMatrices( schur,&A,NULL,&B,&C,&D ) );
#endif

            Mat Bcopy;
            this->check( MatDuplicate(B,MAT_COPY_VALUES,&Bcopy) );

            Vec scaleDiag;
#if PETSC_VERSION_LESS_THAN(3,6,0)
            this->check( MatGetVecs(A,&scaleDiag,NULL) );
#else
            this->check( MatCreateVecs(A,&scaleDiag,NULL) );
#endif
            if ( false ) //this->precFeel()->hasAuxiliarySparseMatrix("mass-matrix") )
            {
                //std::cout << "hasAuxiliarySparseMatrix\n";
                auto massMat = this->precFeel()->auxiliarySparseMatrix("mass-matrix");
                MatrixPetsc<double> * massMatPetsc   = const_cast<MatrixPetsc<double> *>( dynamic_cast<MatrixPetsc<double> const*>( &(*massMat) ) );
                this->check( MatGetDiagonal(massMatPetsc->mat(),scaleDiag) );
            }
            else
            {
                this->check( MatGetDiagonal(A,scaleDiag) );
            }
            this->check( VecReciprocal(scaleDiag) );
            this->check( MatDiagonalScale( Bcopy, scaleDiag ,NULL) );

            //std::cout << "rebuild schur prec\n";

            //MatMatMultSymbolic(C,B,PETSC_DEFAULT,SchurMat);

            this->check( MatMatMult(C,Bcopy,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&schurMatPrecond) );
            //MatView(schurMatPrecond,PETSC_VIEWER_STDOUT_WORLD);

            if ( D != NULL )
            {
#if 0
                this->check( MatAYPX( schurMatPrecond,-1,D,MatStructure::DIFFERENT_NONZERO_PATTERN ) );
#else
                this->check( MatScale( schurMatPrecond, -1 ) );
                this->check( MatAYPX( schurMatPrecond,1,D,MatStructure::DIFFERENT_NONZERO_PATTERN ) );
#endif
            }
            else
              this->check( MatScale( schurMatPrecond, -1 ) );

            // update KSPSetOperators of jac->kspschur
            this->check( PetscImpl::PCFieldSplit_UpdateMatPrecondSchurComplement( pc, schurMatPrecond ) );
            // clean temporary mat and vec
            this->check( MatDestroy( &Bcopy ) );
            this->check( VecDestroy( &scaleDiag ) );
        }
#endif
#if PETSC_VERSION_LESS_THAN(3,5,0)
        this->check( PCFieldSplitSchurPrecondition( pc, theSchurPrecond, schurMatPrecond/*NULL*/ ) );
#else
        this->check( PCFieldSplitSetSchurPre( pc, theSchurPrecond, schurMatPrecond/*NULL*/ ) );
#endif
        // need to call MatDestroy because PCFieldSplitSchurPrecondition call PetscObjectReference ( which increase the object counter)
        // if we not call this  MatDestroy, we have a memory leak
        this->check( MatDestroy( &schurMatPrecond ) );

    }
#endif
    // call necessary before next seting
    this->check( PCSetUp( pc ) );


    // To store array of local KSP contexts on this processor
    KSP* subksps;
    int nSplit;
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,4,0 )
    if ( M_type == "schur" )
        this->check( PetscImpl::PCFieldSplitGetSubKSP_FieldSplit_Schur(pc,&nSplit,&subksps) );
    else
        this->check( PCFieldSplitGetSubKSP(pc,&nSplit,&subksps ) );
#else
        this->check( PCFieldSplitGetSubKSP(pc,&nSplit,&subksps ) );
#endif


    if ( M_type == "schur" )
    {
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,4,0 )
        std::string pcctx = (this->sub().empty())? "" : this->sub()+"-";
        std::string prefixSchurInnerSolver = prefixvm( this->prefix(),pcctx+"fieldsplit-schur-inner-solver" );
        bool noBuildInnerSolver = option(_name="use-outer-solver",_prefix=prefixSchurInnerSolver,_vm=this->vm()).as<bool>();
        if ( !noBuildInnerSolver )
        {
            std::vector<std::string> prefixSchurInnerSolverOverwrite;

            KSP kspInner;
            this->check( PetscImpl::PCFieldSplit_GetKSPInnerSchur( pc, kspInner ) );
            ConfigureKSP kspConf( kspInner, this->precFeel(), this->worldComm(), "", prefixSchurInnerSolver,prefixSchurInnerSolverOverwrite,"preonly", 1e-5,  10 );

            // setup sub-ksp
            this->check( KSPSetUp( kspInner ) );
            //-----------------------------------------------------------//
            // get sub-pc
            PC pcInner;
            this->check( KSPGetPC( kspInner, &pcInner ) );
            // configure sub-pc
            std::string M_innerSchurPCType = option(_name="pc-type",_prefix=prefixSchurInnerSolver,_vm=this->vm()).as<std::string>();
            std::string M_innerSchurPCFactMatSolverPackage = option(_name="pc-factor-mat-solver-package-type",_prefix=prefixSchurInnerSolver,_vm=this->vm()).as<std::string>();
            SetPCType( pcInner, pcTypeConvertStrToEnum( M_innerSchurPCType ),
                       matSolverPackageConvertStrToEnum( M_innerSchurPCFactMatSolverPackage ),
                       this->worldComm() );
            ConfigurePC( pcInner, this->precFeel(), this->worldComm(), "", prefixSchurInnerSolver, prefixSchurInnerSolverOverwrite, this->vm() );
            // setup sub-pc
            this->check( PCSetUp( pcInner ) );
            //-----------------------------------------------------------//
#if 0
            // ksp and pc view
            if ( kspConf.kspView() )
                this->check( KSPView( kspInner, PETSC_VIEWER_STDOUT_WORLD ) );
            else if ( M_subPCview )
                this->check( PCView( pcInner, PETSC_VIEWER_STDOUT_WORLD ) );
#endif
        }

        if ( M_schurFactType == "full" )
        {
            //std::string pcctx = (this->sub().empty())? "" : this->sub()+"-";
            std::string prefixSchurUpperSolver = prefixvm( this->prefix(),pcctx+"fieldsplit-schur-upper-solver" );
            bool noBuildUpperSolver = option(_name="use-outer-solver",_prefix=prefixSchurUpperSolver,_vm=this->vm()).as<bool>();
            if ( !noBuildUpperSolver )
            {
                std::vector<std::string> prefixSchurUpperSolverOverwrite;

                KSP kspUpper;
                this->check( PetscImpl::PCFieldSplit_GetKSPUpperSchur( pc, kspUpper ) );
                ConfigureKSP kspConf( kspUpper, this->precFeel(), this->worldComm(), "", prefixSchurUpperSolver, prefixSchurUpperSolverOverwrite, "preonly", 1e-5,  10 );
                // setup sub-ksp
                this->check( KSPSetUp( kspUpper ) );
                //-----------------------------------------------------------//
                // get sub-pc
                PC pcUpper;
                this->check( KSPGetPC( kspUpper, &pcUpper ) );
                // configure sub-pc
                std::string M_upperSchurPCType = option(_name="pc-type",_prefix=prefixSchurUpperSolver,_vm=this->vm()).as<std::string>();
                std::string M_upperSchurPCFactMatSolverPackage = option(_name="pc-factor-mat-solver-package-type",_prefix=prefixSchurUpperSolver,_vm=this->vm()).as<std::string>();
                SetPCType( pcUpper, pcTypeConvertStrToEnum( M_upperSchurPCType ),
                           matSolverPackageConvertStrToEnum( M_upperSchurPCFactMatSolverPackage ),
                           this->worldComm() );
                ConfigurePC( pcUpper, this->precFeel(), this->worldComm(), "", prefixSchurUpperSolver, prefixSchurUpperSolverOverwrite, this->vm() );
                // setup sub-pc
                this->check( PCSetUp( pcUpper ) );
                //-----------------------------------------------------------//
#if 0
                // ksp and pc view
                if ( kspConf.kspView() )
                    this->check( KSPView( kspUpper, PETSC_VIEWER_STDOUT_WORLD ) );
                else if ( M_subPCview )
                    this->check( PCView( pcUpper, PETSC_VIEWER_STDOUT_WORLD ) );
#endif
            }
        }
#endif // PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,4,0 )
    }

    // config sub ksp/pc for each split
    ConfigurePCFieldSplit::ConfigureSubKSP( &subksps/*pc*/,nSplit,this->precFeel(), M_type,this->worldComm(),this->sub(),this->prefix() );
}

/**
 * ConfigurePCFieldSplitSubKSP
 */
ConfigurePCFieldSplit::ConfigureSubKSP::ConfigureSubKSP( KSP ** subksps/*PC& pc*/, int nSplit,PreconditionerPetsc<double> * precFeel,
                                                         std::string const& typeFieldSplit, WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix, getOptionsDescSplitInFieldSplit(nSplit,prefix,sub) ),
    M_nSplit( nSplit ),
    M_typeFieldSplit( typeFieldSplit )
{
#if 0
    // call necessary before PCFieldSplitGetSubKSP
    this->check( PCSetUp( pc ) );

    // To store array of local KSP contexts on this processor
    KSP* subksps;
    this->check( PCFieldSplitGetSubKSP(pc,&M_nSplit,&subksps ) );
#endif
    M_prefixSplit.resize(M_nSplit);
    M_subPCview.resize(M_nSplit);
    M_subPCtype.resize(M_nSplit);
    M_subMatSolverPackage.resize(M_nSplit);
    for ( int i=0; i<M_nSplit; ++i )
    {
        std::string prefixSplit = prefixvm(this->prefix() , (boost::format( "fieldsplit-%1%" )  %i ).str() );
        M_prefixSplit[i] = prefixSplit;
        M_subPCview[i] = option(_name="pc-view",_prefix=prefixSplit,_vm=this->vm()).as<bool>();
        M_subPCtype[i] = option(_name="pc-type",_prefix=prefixSplit,_vm=this->vm()).as<std::string>();
        M_subMatSolverPackage[i] = option(_name="pc-factor-mat-solver-package-type",_prefix=prefixSplit,_vm=this->vm()).as<std::string>();
    }

    // Loop over sub-ksp objects
    for ( int splitId=0; splitId<M_nSplit; ++splitId )
    {
        VLOG(2) << "configure split " << splitId << " with prefix "<< M_prefixSplit[splitId] << "\n";
        google::FlushLogFiles(google::INFO);

        run( (*subksps)[splitId], splitId );
    }
}

void
ConfigurePCFieldSplit::ConfigureSubKSP::run(KSP& ksp, int splitId )
{
    std::string prefixSplit = M_prefixSplit[splitId];
    std::vector<std::string> prefixSplitOverwrite;
#if 0
    Mat A00;Mat A01;Mat A10; Mat A11;
    this->check( PCFieldSplitGetSchurBlocks(pc,&A00,&A01,&A10, &A11) );
    if (i==0)
        this->check( KSPSetOperators( ksp, A00, A00,
                                      PetscGetMatStructureEnum(MatrixStructure::SAME_PRECONDITIONER)) );
#endif
    ConfigureKSP kspConf( ksp, this->precFeel(), this->worldComm(), "", prefixSplit,prefixSplitOverwrite,"preonly", 1e-5, 50 );
    /*int ierr = KSPSetInitialGuessNonzero ( ksp, PETSC_TRUE );
    CHKERRABORT( this->worldComm().globalComm(),ierr );
    this->check( KSPSetNormType( ksp, KSP_NORM_NONE ) );*/
    // setup ksp
#if PETSC_VERSION_LESS_THAN(3,5,0) || PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,6,0 )
    this->check( KSPSetUp( ksp ) );
#else
    if( M_typeFieldSplit != "schur" || splitId == 0 )
        this->check( KSPSetUp( ksp ) );
#endif
    PC subpc;
    // get sub-pc
    this->check( KSPGetPC( ksp, &subpc ) );
    // configure sub PC
    SetPCType( subpc, pcTypeConvertStrToEnum( M_subPCtype[splitId] ),
               matSolverPackageConvertStrToEnum( M_subMatSolverPackage[splitId] ),
               this->worldComm() );

    // in case of fieldsplit in fieldsplit, need to pass the index split corresponding
    if ( M_subPCtype[splitId] == "fieldsplit" )
    {
#if (PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR >= 2)
        IS isTest;
        this->check( PCFieldSplitGetIS( subpc,(boost::format("%1%")%splitId).str().c_str(),&isTest ) );
        if ( isTest == NULL || this->precFeel()->indexSplitHasChanged() )
        {
            std::string fieldsDefStr = option(_name="fieldsplit-fields",_prefix=prefixSplit,_vm=this->vm()).as<std::string>();
            bool useComponentsSplit = option( _name="fieldsplit-use-components",_prefix=prefixSplit,_vm=this->vm() ).as<bool>();
            auto fieldsDef = IndexSplit::parseFieldsDef( fieldsDefStr );
            PreconditionerPetsc<double>::indexsplit_ptrtype isUsed;
            if ( useComponentsSplit )
                isUsed = this->precFeel()->matrix()->mapRow().indexSplitWithComponents()->applyFieldsDef( IndexSplit::FieldsDef( fieldsDef ) );
            else
                isUsed = this->precFeel()->matrix()->mapRow().indexSplit()->applyFieldsDef( IndexSplit::FieldsDef( fieldsDef ) );
            //isUsed->showMe();

            std::vector<IS> isPetsc;
            PetscConvertIndexSplit( isPetsc ,*isUsed,this->worldComm());
            for ( int i = 0 ; i < isPetsc.size(); ++i )
            {
                this->check( PCFieldSplitSetIS( subpc,(boost::format("%1%")%i).str().c_str(),isPetsc[i] ) );
                this->check( ISDestroy(&isPetsc[i]) );
            }
        }
#else
        CHECK(false) << "use fieldsplit requiert a petsc version greater or equal to 3.2\n";
#endif
    }

    if ( M_subPCtype[splitId] == "lsc" )
        CHECK( splitId==1 ) << "lsc must be use with only field 1, not " << splitId << "\n";

    // force to setup the schur preconditioner
    if( M_typeFieldSplit == "schur" && splitId == 1 )
    {
        if ( subpc->setupcalled )
            PetscObjectStateIncrease((PetscObject)subpc->pmat);
    }

    // configure sub-pc
    ConfigurePC( subpc, this->precFeel(), this->worldComm(), "", prefixSplit,prefixSplitOverwrite,this->vm() );

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
ConfigurePCLSC::ConfigurePCLSC( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                std::string const& sub, std::string const& prefix, std::string lscVersion )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix,getOptionsDescLSC(prefix,sub) ),
    M_version( lscVersion ),
    M_prefixLSC( prefixvm(this->prefix(),"lsc") ),
    M_scaleDiag( option(_name="scale-diag",_prefix=M_prefixLSC,_vm=this->vm()).as<bool>() ),
    M_subPCtype( option(_name="pc-type",_prefix=M_prefixLSC,_vm=this->vm()).as<std::string>() ),
    M_subMatSolverPackage( option(_name="pc-factor-mat-solver-package-type",_prefix=M_prefixLSC,_vm=this->vm()).as<std::string>() ),
    M_subPCview( option(_name="pc-view",_prefix=M_prefixLSC,_vm=this->vm()).as<bool>() )
{
    if ( M_version == "petsc" )
        this->check( PetscImpl::PCLSCSetScaleDiag( pc, static_cast<PetscBool>( M_scaleDiag ) ) );
    else
        this->check( PCLSC2SetScaleDiag( pc, static_cast<PetscBool>( M_scaleDiag ) ) );

    if ( M_version != "petsc" )
    {
        if ( this->precFeel()->hasAuxiliarySparseMatrix("mass-matrix") )
        {
            //std::cout << "LSC hasAuxiliarySparseMatrix\n";
            auto massMat = this->precFeel()->auxiliarySparseMatrix("mass-matrix");
            MatrixPetsc<double> * massMatPetsc   = const_cast<MatrixPetsc<double> *>( dynamic_cast<MatrixPetsc<double> const*>( &(*massMat) ) );
            this->check( PCSetMassMatrix_LSC2( pc, massMatPetsc->mat()) );
        }
    }


    VLOG(2) << "ConfigurePC : LSC\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->prefixLSC : " << M_prefixLSC  << "\n"
            << "  |->subPCtype : " << M_subPCtype << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}
void
ConfigurePCLSC::run( PC& pc )
{
    std::vector<std::string> prefixOverwrite;

    // setup sub-pc
    this->check( PCSetUp( pc ) );
    //-----------------------------------------------------------//
    // get sub-ksp
    KSP subksp;
    if ( M_version == "petsc" )
        this->check( PetscImpl::PCLSCGetKSP( pc, subksp ) );
    else
        this->check( PCLSC2GetKSP( pc, subksp ) );
    // configure sub-ksp
    ConfigureKSP kspConf( subksp, this->precFeel(), this->worldComm(), "", M_prefixLSC,prefixOverwrite,"preonly", 1e-5, 50 );
    // setup sub-ksp
    this->check( KSPSetUp( subksp ) );
    //-----------------------------------------------------------//
    // get sub-pc
    PC subpc;
    this->check( KSPGetPC( subksp, &subpc ) );
    // configure sub-pc
    SetPCType( subpc, pcTypeConvertStrToEnum( M_subPCtype ),
               matSolverPackageConvertStrToEnum( M_subMatSolverPackage ),
               this->worldComm() );
    ConfigurePC( subpc, this->precFeel(), this->worldComm(), "", M_prefixLSC,prefixOverwrite,this->vm() );
    // setup sub-pc
    this->check( PCSetUp( subpc ) );
    //-----------------------------------------------------------//
    // ksp and pc view
    if ( kspConf.kspView() )
        this->check( KSPView( subksp, PETSC_VIEWER_STDOUT_WORLD ) );
    else if ( M_subPCview )
        this->check( PCView( subpc, PETSC_VIEWER_STDOUT_WORLD ) );
}

/**
 * ConfigurePCPMM
 */
ConfigurePCPMM::ConfigurePCPMM( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix,getOptionsDescPMM(prefix,sub) ),
    M_prefixPMM( prefixvm(this->prefix(),"pmm") ),
    M_subPCtype( option(_name="pc-type",_prefix=M_prefixPMM,_vm=this->vm()).as<std::string>() ),
    M_subMatSolverPackage( option(_name="pc-factor-mat-solver-package-type",_prefix=M_prefixPMM,_vm=this->vm()).as<std::string>() ),
    M_subPCview( option(_name="pc-view",_prefix=M_prefixPMM,_vm=this->vm()).as<bool>() )
{
    CHECK( this->precFeel()->hasAuxiliarySparseMatrix("pmm") ) << "pmm is not given";

    auto massMat = this->precFeel()->auxiliarySparseMatrix("pmm");
    MatrixPetsc<double> * massMatPetsc   = const_cast<MatrixPetsc<double> *>( dynamic_cast<MatrixPetsc<double> const*>( &(*massMat) ) );
    this->check( PCSetPressureMassMatrix_PMM_Feelpp( pc, massMatPetsc->mat()) );


    VLOG(2) << "ConfigurePC : PMM\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->prefixPMM : " << M_prefixPMM  << "\n"
            << "  |->subPCtype : " << M_subPCtype << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}

void
ConfigurePCPMM::run( PC& pc )
{
    std::vector<std::string> prefixOverwrite;

    // setup sub-pc
    this->check( PCSetUp( pc ) );

    //-----------------------------------------------------------//
    // get sub-ksp
    KSP subksp;
    this->check( PCGetKSP_PMM_Feelpp( pc, subksp ) );
    // configure sub-ksp
    ConfigureKSP kspConf( subksp, this->precFeel(), this->worldComm(), "", M_prefixPMM,prefixOverwrite,"preonly", 1e-5, 50 );
    // setup sub-ksp
    this->check( KSPSetUp( subksp ) );
    //-----------------------------------------------------------//
    // get sub-pc
    PC subpc;
    this->check( KSPGetPC( subksp, &subpc ) );
    // configure sub-pc
    SetPCType( subpc, pcTypeConvertStrToEnum( M_subPCtype ),
               matSolverPackageConvertStrToEnum( M_subMatSolverPackage ),
               this->worldComm() );
    ConfigurePC( subpc, this->precFeel(), this->worldComm(), "", M_prefixPMM,prefixOverwrite,this->vm() );
    // setup sub-pc
    this->check( PCSetUp( subpc ) );
    //-----------------------------------------------------------//
    // ksp and pc view
    if ( kspConf.kspView() )
        this->check( KSPView( subksp, PETSC_VIEWER_STDOUT_WORLD ) );
    else if ( M_subPCview )
        this->check( PCView( subpc, PETSC_VIEWER_STDOUT_WORLD ) );
}


/**
 * ConfigurePCPCD
 */
ConfigurePCPCD::ConfigurePCPCD( PC& pc, PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                std::string const& sub, std::string const& prefix )
    :
    ConfigurePCBase( precFeel, worldComm,sub,prefix,getOptionsDescPCD(prefix,sub) ),
    M_prefixPCD_Ap( prefixvm(this->prefix(),"pcd.Ap") ),
    M_prefixPCD_Mp( prefixvm(this->prefix(),"pcd.Mp") ),
    M_subPCtype_Ap( option(_name="pc-type",_prefix=M_prefixPCD_Ap,_vm=this->vm()).as<std::string>() ),
    M_subPCtype_Mp( option(_name="pc-type",_prefix=M_prefixPCD_Mp,_vm=this->vm()).as<std::string>() ),
    M_subMatSolverPackage_Ap( option(_name="pc-factor-mat-solver-package-type",_prefix=M_prefixPCD_Ap,_vm=this->vm()).as<std::string>() ),
    M_subMatSolverPackage_Mp( option(_name="pc-factor-mat-solver-package-type",_prefix=M_prefixPCD_Mp,_vm=this->vm()).as<std::string>() ),
    M_subPCview_Ap( option(_name="pc-view",_prefix=M_prefixPCD_Ap,_vm=this->vm()).as<bool>() ),
    M_subPCview_Mp( option(_name="pc-view",_prefix=M_prefixPCD_Mp,_vm=this->vm()).as<bool>() )
{
    CHECK( this->precFeel()->hasOperatorPCD("pcd") ) << "operator pcd is not given";

    auto opPCD = this->precFeel()->operatorPCD("pcd");
    MatrixPetsc<double> * pmMatPetsc   = const_cast<MatrixPetsc<double> *>( dynamic_cast<MatrixPetsc<double> const*>( &(*(opPCD->pressureMassMatrix())) ) );
    MatrixPetsc<double> * pdcMatPetsc   = const_cast<MatrixPetsc<double> *>( dynamic_cast<MatrixPetsc<double> const*>( &(*(opPCD->pressureDiffusionConvectionMatrix())) ) );
    this->check( PCSetMatMp_PCD_Feelpp( pc, pmMatPetsc->mat()) );
    this->check( PCSetMatFp_PCD_Feelpp( pc, pdcMatPetsc->mat()) );
    this->check( PCSetOrder_PCD_Feelpp( pc, opPCD->pcdOrder() ) );
    this->check( PCSetApType_PCD_Feelpp( pc, opPCD->pcdDiffusionType().c_str() ) );
    if ( opPCD->pressureLaplacianMatrix() )
    {
        MatrixPetsc<double> * plMatPetsc   = const_cast<MatrixPetsc<double> *>( dynamic_cast<MatrixPetsc<double> const*>( &(*(opPCD->pressureLaplacianMatrix())) ) );
        this->check( PCSetMatApLaplacian_PCD_Feelpp( pc, plMatPetsc->mat()) );
    }
    if ( opPCD->velocityMassMatrix() )
    {
        MatrixPetsc<double> * vmMatPetsc   = const_cast<MatrixPetsc<double> *>( dynamic_cast<MatrixPetsc<double> const*>( &(*(opPCD->velocityMassMatrix())) ) );
        this->check( PCSetMatMv_PCD_Feelpp( pc, vmMatPetsc->mat()) );
    }

    VLOG(2) << "ConfigurePC : PCD\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->prefixPCD_Ap : " << M_prefixPCD_Ap  << "\n"
            << "  |->subPCtype_Ap : " << M_subPCtype_Ap << "\n"
            << "  |->prefixPCD_Mp : " << M_prefixPCD_Mp  << "\n"
            << "  |->subPCtype_Mp : " << M_subPCtype_Mp << "\n";
    google::FlushLogFiles(google::INFO);
    run( pc );
}

void
ConfigurePCPCD::run( PC& pc )
{
    std::vector<std::string> prefixOverwrite;

    // setup sub-pc
    this->check( PCSetUp( pc ) );

    //-----------------------------------------------------------//
    // get sub-ksp
    KSP subksp_Ap, subksp_Mp;
    this->check( PCGetKSP_Ap_PCD_Feelpp( pc, subksp_Ap ) );
    this->check( PCGetKSP_Mp_PCD_Feelpp( pc, subksp_Mp ) );
    // configure sub-ksp
    ConfigureKSP kspConf_Ap( subksp_Ap, this->precFeel(), this->worldComm(), "", M_prefixPCD_Ap,prefixOverwrite,"preonly", 1e-5, 50 );
    ConfigureKSP kspConf_Mp( subksp_Mp, this->precFeel(), this->worldComm(), "", M_prefixPCD_Mp,prefixOverwrite,"preonly", 1e-5, 50 );
    // setup sub-ksp
    this->check( KSPSetUp( subksp_Ap ) );
    this->check( KSPSetUp( subksp_Mp ) );
    //-----------------------------------------------------------//
    // get sub-pc
    PC subpc_Ap, subpc_Mp;
    this->check( KSPGetPC( subksp_Ap, &subpc_Ap ) );
    this->check( KSPGetPC( subksp_Mp, &subpc_Mp ) );
    // configure sub-pc
    SetPCType( subpc_Ap, pcTypeConvertStrToEnum( M_subPCtype_Ap ),
               matSolverPackageConvertStrToEnum( M_subMatSolverPackage_Ap ),
               this->worldComm() );
    SetPCType( subpc_Mp, pcTypeConvertStrToEnum( M_subPCtype_Mp ),
               matSolverPackageConvertStrToEnum( M_subMatSolverPackage_Mp ),
               this->worldComm() );
    ConfigurePC( subpc_Ap, this->precFeel(), this->worldComm(), "", M_prefixPCD_Ap,prefixOverwrite,this->vm() );
    ConfigurePC( subpc_Mp, this->precFeel(), this->worldComm(), "", M_prefixPCD_Mp,prefixOverwrite,this->vm() );
    // setup sub-pc
    this->check( PCSetUp( subpc_Ap ) );
    this->check( PCSetUp( subpc_Mp ) );
    //-----------------------------------------------------------//
    // ksp and pc view
    if ( kspConf_Ap.kspView() )
        this->check( KSPView( subksp_Ap, PETSC_VIEWER_STDOUT_WORLD ) );
    else if ( M_subPCview_Ap )
        this->check( PCView( subpc_Ap, PETSC_VIEWER_STDOUT_WORLD ) );
    if ( kspConf_Mp.kspView() )
        this->check( KSPView( subksp_Mp, PETSC_VIEWER_STDOUT_WORLD ) );
    else if ( M_subPCview_Mp )
        this->check( PCView( subpc_Mp, PETSC_VIEWER_STDOUT_WORLD ) );
}


/**
 * ConfigurePCRedundant
 */
ConfigurePCRedundant::ConfigurePCRedundant( PreconditionerPetsc<double> * precFeel, WorldComm const& worldComm,
                                            std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
    :
    ConfigurePCBase( precFeel, worldComm,"",prefix,prefixOverwrite, getOptionsDescPrecBase( prefixvm(prefix,"redundant"),"" ) ),
    M_innerPCtype( option(_name="pc-type",_prefix= prefixvm(prefix,"redundant"),_vm=this->vm()).as<std::string>() ),
    M_innerPCMatSolverPackage(option(_name="pc-factor-mat-solver-package-type",_prefix= prefixvm(prefix,"redundant"),_vm=this->vm()).as<std::string>() ),
    M_innerPCview( option(_name="pc-view",_prefix= prefixvm(prefix,"redundant"),_vm=this->vm()).as<bool>() ),
    M_factorShiftType( "none" )
{
    VLOG(2) << "ConfigurePC : Redundant\n"
            << "  |->prefix    : " << this->prefix() << std::string((this->sub().empty())? "" : " -sub="+this->sub()) << "\n"
            << "  |->innerPCtype : " << M_innerPCtype << "\n"
            << "  |->innerPCMatSolverPackage : " << M_innerPCMatSolverPackage << "\n";
    google::FlushLogFiles(google::INFO);
}
void
ConfigurePCRedundant::run( PC& pc )
{
#if PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,5,0 )
    // redifine PCSetUp for PCREDUNDANT because originaly KSPSetUp for innerksp is called in this function
    this->check( PetscImpl::PCRedundantChangeSetup(pc) );
    // build operators
    this->check( PCSetUp( pc ) );

    std::vector<std::string> prefixOverwrite;

    KSP innerksp;
    this->check( PCRedundantGetKSP(pc,&innerksp) );

    PC innerpc;
    this->check( KSPGetPC(innerksp,&innerpc) );

    // configure sequential pc
    SetPCType( innerpc, pcTypeConvertStrToEnum( M_innerPCtype ),
               matSolverPackageConvertStrToEnum( M_innerPCMatSolverPackage ),
               this->worldComm().subWorldCommSeq() );
    ConfigurePC mypc( this->precFeel(), /*innerpc,*/ this->worldComm(), "", prefixvm(this->prefix(),"redundant"), prefixOverwrite, this->vm() );
    mypc.setFactorShiftType( M_factorShiftType );
    mypc.run( innerpc );
    // setup inner pc
    this->check( PCSetUp( innerpc ) );

    // configure inner ksp
    ConfigureKSP kspConf( innerksp, this->precFeel(), this->worldComm(), "", prefixvm(this->prefix(),"redundant"), prefixOverwrite, "preonly",1e-8,50 );
    // setup inner ksp
    this->check( KSPSetUp( innerksp ) );

    PetscViewer viewer = (this->sub().empty())? PETSC_VIEWER_STDOUT_WORLD : PETSC_VIEWER_STDOUT_SELF;
    if ( kspConf.kspView() )
        this->check( KSPView( innerksp, viewer ) );
    else if ( M_innerPCview )
        this->check( PCView( innerpc, viewer ) );
#else
    CHECK( false ) << "redundant not supported for this PETSc version ( must be >= 3.5 )";
#endif
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

        if ((subpctype=="lu") || (subpctype=="cholesky"))
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
