/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-27

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef FEELPP_ALG_ENUMS_HPP
#define FEELPP_ALG_ENUMS_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelpetsc.hpp>


namespace Feel
{
/**
 * Context for 'on' operation on sparse matrices
 */
struct ContextOn
{
    enum Options
    {
        PENALISATION                = 0x0,  /**< penalisation */
        ELIMINATION                 = 0x1, /**< elimination */
        KEEP_DIAGONAL   = 0x2, /**< enables elimination and keep diagonal entry(ie don't put 1), modify rhs accordingly */
        SYMMETRIC       = 0x4,  /**< enables elimination and make a symmetric elimination */
        CHECK           = 0x6
    };
};
extern std::map<std::string, size_type> ContextOnMap;


enum   MatrixProperties
{
    HERMITIAN          = 0x1,   /**< hermitian : \f$A^* = A\f$ */
    NON_HERMITIAN      = 0x2,   /**< non hermitian : \f$A^* != A\f$ */
    POSITIVE_DEFINITE  = 0x4,   /**< positive definite matrix : \f$v^* A v > 0 \f$ for all non-zero v */
    SINGULAR           = 0x8,    /**< singular matrix : \f$det(A)=0\f$ and 0 is an eigenvalue */
    DENSE              = 0x10    /**< dense matrix */
};

enum MatrixTranspose
{
    MATRIX_TRANSPOSE_ASSEMBLED   = 0x1,
    MATRIX_TRANSPOSE_UNASSEMBLED = 0x2,
    MATRIX_TRANSPOSE_CHECK       = 0x4
};
        

/**
 * Backend types
 *
 * At the moment, we support GMM(serial), PETSC and TRILINOS(serial and parallel)
 */
enum BackendType
{
    BACKEND_GMM = 0,
    BACKEND_PETSC,
    BACKEND_TRILINOS,
    BACKEND_EIGEN,
    BACKEND_EIGEN_DENSE
};

/**
 * Defines an \p enum for iterative solver types
 */
enum SolverType {CG=0,
                 CGN,
                 CGS,
                 CR,
                 QMR,
                 TCQMR,
                 TFQMR,
                 BICG,
                 BICGSTAB,
                 MINRES,
                 GMRES,
                 FGMRES,
                 LSQR,
                 JACOBI,
                 SOR_FORWARD,
                 SOR_BACKWARD,
                 SSOR,
                 RICHARDSON,
                 CHEBYSHEV,
                 PREONLY,
                 GCR,
                 INVALID_SOLVER
                };

/**
 * Defines an \p enum for preconditioner types
 */
enum PreconditionerType {IDENTITY_PRECOND =0,
                         JACOBI_PRECOND,
                         BLOCK_JACOBI_PRECOND,
                         SOR_PRECOND,
                         SSOR_PRECOND,
                         EISENSTAT_PRECOND,
                         ASM_PRECOND,
                         GASM_PRECOND,
                         CHOLESKY_PRECOND,
                         ICC_PRECOND,
                         ILU_PRECOND,
                         LU_PRECOND,
                         AMG_PRECOND,
                         USER_PRECOND,
                         SHELL_PRECOND,
                         FIELDSPLIT_PRECOND,
                         LSC_PRECOND,
                         ML_PRECOND,
                         GAMG_PRECOND,
                         BOOMERAMG_PRECOND,
                         REDUNDANT_PRECOND,
                         NONE_PRECOND,
                         INVALID_PRECONDITIONER
                        };


/**
 * Defines an \p enum for field split types
 */
enum FieldSplitType {ADDITIVE=0,
                     MULTIPLICATIVE,
                     SCHUR
                    };

/**
 * indicates the structure of the matrix versus preconditioner
 */
enum MatrixStructure
{
    SAME_NONZERO_PATTERN,
    DIFFERENT_NONZERO_PATTERN,
    SAME_PRECONDITIONER,
    SUBSET_NONZERO_PATTERN,
    INVALID_STRUCTURE
};
/**
 * Defines an \p enum for iterative eigenproblem solver types
 */
enum EigenSolverType
{
    POWER=0,
    LAPACK,
    SUBSPACE,
    ARNOLDI,
    LANCZOS,
    KRYLOVSCHUR,
    // SLEPc optional packages
    ARPACK,
    // EPSBLZPACK,
    // EPSPLANSO,
    // EPSTRLAN,

    INVALID_EIGENSOLVER
}; // EigenSolverType

/**
 * Defines an \p enum for eigenproblem types.  This can be Hermitian (HEP),
 * generalized Hermitian (GHEP), non-Hermitian (NHEP), generalized non-Hermitian
 * (GNHEP) and Generalized Non-Hermitian GNHEP with positive (semi-)definite B
 */
enum EigenProblemType {NHEP=0,
                       HEP,
                       GNHEP,
                       GHEP,
                       PGNHEP,

                       INVALID_EIGENPROBLEMTYPE
                      };



/**
 * Defines an \p enum for the position of
 * the spectrum, i.e. the eigenvalues to be computed.
 */
enum PositionOfSpectrum {LARGEST_MAGNITUDE=0,
                         SMALLEST_MAGNITUDE,
                         LARGEST_REAL,
                         SMALLEST_REAL,
                         LARGEST_IMAGINARY,
                         SMALLEST_IMAGINARY,

                         INVALID_Postion_of_Spectrum
                        };

/**
 * Spectral transform type
 */
enum SpectralTransformType {SHIFT=0,
                            SINVERT,
                            FOLD,
                            CAYLEY
                           };

extern std::map<std::string, size_type> EigenMap;

/**
 * Defines an \p enum for various linear solver packages.  This
 * allows for run-time switching between solver packages
 *
 */
enum SolverPackage
{
    SOLVERS_FEEL=0,
    SOLVERS_GMM,
    SOLVERS_EIGEN,
    SOLVERS_PETSC,
    SOLVERS_TRILINOS,
    SOLVERS_SLEPC,
    SOLVER_INVALID_PACKAGE
};

/**
 * Define an \p enum for non linear solver type
 * if SELECT_IN_ARGLIST the choice is done by the arguments in the line command
 */
enum SolverNonLinearType
{
    SELECT_IN_ARGLIST=0,
    LINE_SEARCH,
    TRUST_REGION,
    NRICHARDSON,
    NKSPONLY,
    VINEWTONRSLS,
    VINEWTONRSTR,
    NGMRES,
    QN,
    NSHELL,
    GS,
    NCG,
    FAS,
    MS,
    NASM,
    ANDERSON,
    ASPIN
};

/**
 *
 */
enum AitkenType
{
    AITKEN_STANDARD=0,
    AITKEN_METHOD_1=1,
    FIXED_RELAXATION_METHOD=2
};
/**
 *
 */
enum DirichletType
{
    STRONG=0,
    WEAK=1
};

enum ProjectorType
{
    NODAL=-1,
    L2=0,
    H1=1,
    DIFF=2,
    HDIV=3,
    HCURL=4,
    LIFT=5,
    CIP=6
};

enum MatSolverPackageType
{
    MATSOLVER_SPOOLES=0,
    MATSOLVER_SUPERLU,
    MATSOLVER_SUPERLU_DIST,
    MATSOLVER_UMFPACK,
    MATSOLVER_ESSL,
    MATSOLVER_LUSOL,
    MATSOLVER_MUMPS,
    MATSOLVER_MKL_PARDISO,
    MATSOLVER_PASTIX,
    MATSOLVER_DSCPACK,
    MATSOLVER_MATLAB,
    MATSOLVER_PETSC,
    MATSOLVER_PLAPACK,
    MATSOLVER_BAS,
    MATSOLVER_BOOMERAMG,
    MATSOLVER_EUCLID,
    MATSOLVER_PILUT,

};
#if defined(FEELPP_HAS_MUMPS) && PETSC_VERSION_GREATER_OR_EQUAL_THAN( 3,2,0 )
const auto MATSOLVER_DEFAULT = MATSOLVER_MUMPS;
#else
const auto MATSOLVER_DEFAULT = MATSOLVER_PETSC;
#endif

PreconditionerType
pcTypeConvertStrToEnum( std::string const& type );

SolverType
kspTypeConvertStrToEnum( std::string const& type );

SolverNonLinearType
snesTypeConvertStrToEnum( std::string const& type );
std::string
snesTypeConvertEnumToStr( SolverNonLinearType type );

MatSolverPackageType
matSolverPackageConvertStrToEnum( std::string const& type );

FieldSplitType
fieldsplitTypeConvertStrToEnum( std::string const& type );



} // Feel
#endif /* FEELPP_ALG_ENUMS_HPP */
