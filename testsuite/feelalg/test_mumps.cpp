/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#define BOOST_TEST_MODULE mumps testsuite

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelconfig.h>

#if defined(FEELPP_HAS_PETSC_H)
#include <petscksp.h>
#endif

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{

#if defined(FEELPP_HAS_PETSC_H)
void
requirePetscSuccess( PetscErrorCode ierr, char const* where )
{
    if ( ierr )
    {
        char const* text = nullptr;
        PetscErrorMessage( ierr, &text, nullptr );
        BOOST_FAIL( where << " failed with PETSc error " << ierr << ": " << ( text ? text : "unknown error" ) );
    }
}
#endif

} // namespace

BOOST_AUTO_TEST_SUITE( mumpssuite )

BOOST_AUTO_TEST_CASE( distributed_lu_mumps_solve )
{
#if defined(FEELPP_HAS_PETSC_H) && defined(PETSC_HAVE_MUMPS)
    MPI_Comm comm = PETSC_COMM_WORLD;
    int commSize = 1;
    MPI_Comm_size( comm, &commSize );

    PetscInt n = 40 * commSize;
    Mat A = nullptr;
    Vec exact = nullptr;
    Vec rhs = nullptr;
    Vec sol = nullptr;
    KSP ksp = nullptr;
    PC pc = nullptr;

    requirePetscSuccess( PetscPushErrorHandler( PetscReturnErrorHandler, nullptr ), "PetscPushErrorHandler" );

    BOOST_TEST_CHECKPOINT( "creating distributed AIJ matrix" );
    requirePetscSuccess( MatCreateAIJ( comm, PETSC_DECIDE, PETSC_DECIDE, n, n, 3, nullptr, 3, nullptr, &A ), "MatCreateAIJ" );

    PetscInt rowStart = 0;
    PetscInt rowEnd = 0;
    requirePetscSuccess( MatGetOwnershipRange( A, &rowStart, &rowEnd ), "MatGetOwnershipRange" );

    BOOST_TEST_CHECKPOINT( "filling distributed tridiagonal matrix" );
    for ( PetscInt row = rowStart; row < rowEnd; ++row )
    {
        PetscInt cols[3];
        PetscScalar vals[3];
        PetscInt nnz = 0;

        if ( row == 0 || row == n - 1 )
        {
            cols[nnz] = row;
            vals[nnz] = 1.0;
            ++nnz;
        }
        else
        {
            cols[nnz] = row - 1;
            vals[nnz] = -1.0;
            ++nnz;

            cols[nnz] = row;
            vals[nnz] = 2.0;
            ++nnz;

            cols[nnz] = row + 1;
            vals[nnz] = -1.0;
            ++nnz;
        }

        requirePetscSuccess( MatSetValues( A, 1, &row, nnz, cols, vals, INSERT_VALUES ), "MatSetValues" );
    }

    BOOST_TEST_CHECKPOINT( "assembling distributed matrix" );
    requirePetscSuccess( MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY ), "MatAssemblyBegin" );
    requirePetscSuccess( MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY ), "MatAssemblyEnd" );

    BOOST_TEST_CHECKPOINT( "creating vectors" );
    requirePetscSuccess( VecCreateMPI( comm, PETSC_DECIDE, n, &exact ), "VecCreateMPI" );
    requirePetscSuccess( VecDuplicate( exact, &rhs ), "VecDuplicate(rhs)" );
    requirePetscSuccess( VecDuplicate( exact, &sol ), "VecDuplicate(sol)" );
    requirePetscSuccess( VecSet( exact, 1.0 ), "VecSet(exact)" );
    requirePetscSuccess( VecSet( sol, 0.0 ), "VecSet(sol)" );
    requirePetscSuccess( MatMult( A, exact, rhs ), "MatMult" );

    BOOST_TEST_CHECKPOINT( "configuring KSP/PC MUMPS LU" );
    requirePetscSuccess( KSPCreate( comm, &ksp ), "KSPCreate" );
    requirePetscSuccess( KSPSetOperators( ksp, A, A ), "KSPSetOperators" );
    requirePetscSuccess( KSPSetType( ksp, KSPPREONLY ), "KSPSetType" );
    requirePetscSuccess( KSPGetPC( ksp, &pc ), "KSPGetPC" );
    requirePetscSuccess( PCSetType( pc, PCLU ), "PCSetType" );
    requirePetscSuccess( PCFactorSetMatSolverType( pc, MATSOLVERMUMPS ), "PCFactorSetMatSolverType" );

    // Keep MUMPS failures as PETSc errors when possible, instead of aborting MPI.
    PetscBool hasMumpsIcntl21 = PETSC_FALSE;
    requirePetscSuccess( PetscOptionsHasName( nullptr, nullptr, "-mat_mumps_icntl_21", &hasMumpsIcntl21 ), "PetscOptionsHasName(-mat_mumps_icntl_21)" );
    if ( !hasMumpsIcntl21 )
        requirePetscSuccess( PetscOptionsSetValue( nullptr, "-mat_mumps_icntl_21", "0" ), "PetscOptionsSetValue(-mat_mumps_icntl_21)" );
    requirePetscSuccess( KSPSetFromOptions( ksp ), "KSPSetFromOptions" );

    BOOST_TEST_CHECKPOINT( "setting up MUMPS LU factorization" );
    requirePetscSuccess( KSPSetUp( ksp ), "KSPSetUp" );

    BOOST_TEST_CHECKPOINT( "solving with MUMPS LU factorization" );
    PetscErrorCode solveErr = KSPSolve( ksp, rhs, sol );
    requirePetscSuccess( solveErr, "KSPSolve" );

    BOOST_TEST_CHECKPOINT( "checking MUMPS solve result" );
    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    requirePetscSuccess( KSPGetConvergedReason( ksp, &reason ), "KSPGetConvergedReason" );
    BOOST_REQUIRE_GT( reason, 0 );

    requirePetscSuccess( VecAXPY( sol, -1.0, exact ), "VecAXPY" );
    PetscReal err = 0;
    requirePetscSuccess( VecNorm( sol, NORM_INFINITY, &err ), "VecNorm" );
    BOOST_CHECK_SMALL( static_cast<double>( err ), 1e-10 );

    requirePetscSuccess( KSPDestroy( &ksp ), "KSPDestroy" );
    requirePetscSuccess( VecDestroy( &sol ), "VecDestroy(sol)" );
    requirePetscSuccess( VecDestroy( &rhs ), "VecDestroy(rhs)" );
    requirePetscSuccess( VecDestroy( &exact ), "VecDestroy(exact)" );
    requirePetscSuccess( MatDestroy( &A ), "MatDestroy" );
    requirePetscSuccess( PetscPopErrorHandler(), "PetscPopErrorHandler" );
#else
    BOOST_TEST_MESSAGE( "PETSc/MUMPS support is not available; skipping MUMPS reproducer" );
#endif
}

BOOST_AUTO_TEST_SUITE_END()
