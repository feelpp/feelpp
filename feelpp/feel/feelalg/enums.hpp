/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

    Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
    Date: 2005-11-27 (updated 2025-01-18)

    Copyright (C) 2005,2006 EPFL
    Copyright (C) 2006,2007 Universite Joseph Fourier (Grenoble I)
    Copyright (C) 2011-2025 Feel++ Consortium

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

#if defined(FEELPP_HAS_PETSC)
#include <feel/feelcore/feelpetsc.hpp>
#endif

namespace Feel
{

namespace solve {

enum class strategy
{
    //! monolithic solve
    monolithic=0,
    // static condensation solve strategy
    static_condensation=1,
    //! local solve
    local=2
};

} // namespace solve


/**
 * Context for 'on' operation on sparse matrices
 */
struct ContextOn
{
    enum Options
    {
        PENALISATION                = 0x0,  /**< penalisation */
        ELIMINATION                 = 0x1, /**< elimination */
        KEEP_DIAGONAL               = 0x2, /**< enables elimination and keep diagonal entry(ie don't put 1), modify rhs accordingly */
        SYMMETRIC                   = 0x4,  /**< enables elimination and make a symmetric elimination */
        CHECK                       = 0x6
    };
};
extern std::map<std::string, size_type> ContextOnMap;


enum MatrixProperties
{
    SYMMETRIC          = 1 << 1, /**< symmetric : \f$A^T = A */
    STRUCTURALLY_SYMMETRIC = 1 << 2, 
    HERMITIAN          = 1 << 3, /**< hermitian : \f$A^* = A\f$ */
    NON_HERMITIAN      = 1 << 4, /**< non hermitian : \f$A^* != A\f$ */
    POSITIVE_DEFINITE  = 1 << 5, /**< positive definite matrix : \f$v^* A v > 0 \f$ for all non-zero v */
    NEGATIVE_DEFINITE  = 1 << 6, /**< negative definite matrix : \f$v^* A v < 0 \f$ for all non-zero v */
    INDEFINITE         = 1 << 7, /**< negative and positive eigenvalues */
    SPD                = SYMMETRIC | POSITIVE_DEFINITE,
    SINGULAR           = 1 << 11,    /**< singular matrix : \f$det(A)=0\f$ and 0 is an eigenvalue */
    DENSE              = 1 << 20,    /**< dense matrix */
};

enum MatrixTranspose
{
    MATRIX_TRANSPOSE_ASSEMBLED   = 0x1,
    MATRIX_TRANSPOSE_UNASSEMBLED = 0x2,
    MATRIX_TRANSPOSE_CHECK       = 0x4
};

/**
 * Backend types
 */
enum BackendType
{
    BACKEND_NONE = -1,
    BACKEND_GMM,
#if FEELPP_HAS_PETSC
    BACKEND_PETSC,
#endif
    BACKEND_TRILINOS,
    BACKEND_EIGEN,
    BACKEND_EIGEN_DENSE
};

#if FEELPP_HAS_PETSC
const BackendType BACKEND_DEFAULT = BACKEND_PETSC;
#else
const BackendType BACKEND_DEFAULT = BACKEND_EIGEN;
#endif

/**
 * Defines an enum for iterative solver types (KSP)
 *
 * Updated to include newer PETSc methods such as pipelined variants,
 * deflated GMRES, and others.  Add as needed.
 */
enum SolverType
{
    // Existing solvers in Feel++ (kept in the same numeric order)
    CG=0,            // KSPCG               -> "cg"
    CGN,             // KSPCGNE             -> "cgne"
    CGS,             // KSPCGS              -> "cgs"
    CR,              // KSPCR               -> "cr"
    QMR,             // (Legacy; PETSc has KSPQMRCGS => "qmrcgs" not exactly the same)
    TCQMR,           // KSPTCQMR            -> "tcqmr"
    TFQMR,           // KSPTFQMR            -> "tfqmr"
    BICG,            // KSPBICG             -> "bicg"
    BICGSTAB,        // KSPBCGS             -> "bcgs"
    MINRES,          // KSPMINRES           -> "minres"
    GMRES,           // KSPGMRES            -> "gmres"
    FGMRES,          // KSPFGMRES           -> "fgmres"
    LSQR,            // KSPLSQR             -> "lsqr"

    // These four are not actually "KSP" methods in PETSc but historically 
    // included in Feel++ for solver selection; in PETSc, Jacobi & SOR are PC types.
    JACOBI,
    SOR_FORWARD,
    SOR_BACKWARD,
    SSOR,

    RICHARDSON,      // KSPRICHARDSON       -> "richardson"
    CHEBYSHEV,       // KSPCHEBYSHEV        -> "chebyshev"
    PREONLY,         // KSPPREONLY          -> "preonly"
    GCR,             // KSPGCR              -> "gcr"

    // Existing “newer” pipeline & deflated methods already in Feel++
    DGMRES,          // KSPDGMRES           -> "dgmres"
    LGMRES,          // KSPLGMRES           -> "lgmres"
    PGMRES,          // KSPPGMRES           -> "pgmres"
    PIPEGMRES,       // (Historically: could map to "pipefgmres" or "pgmres" in PETSc)
    PIPECG,          // KSPPIPECG           -> "pipecg"
    PIPECR,          // KSPPIPECR           -> "pipecr"
    BCGSL,           // KSPBCGSL            -> "bcgsl"
    FBCGS,           // KSPFBCGS            -> "fbcgs"
    IBCGS,           // KSPIBCGS            -> "ibcgs"
    FCG,             // KSPFCG              -> "fcg"

    // Newly added from PETSc table (appended at the end)
    FBCGSR,    // KSPFBCGSR    -> "fbcgsr"
    GROPPCG,   // KSPGROPPCG   -> "groppcg"
    PIPECGRR,  // KSPPIPECGRR  -> "pipecgrr"
    PIPEFCG,   // KSPPIPEFCG   -> "pipefcg"
    CGLS,      // KSPCGLS      -> "cgls"
    NASH,      // KSPNASH      -> "nash"
    STCG,      // KSPSTCG      -> "stcg"
    GLTR,      // KSPGLTR      -> "gltr"
    QCG,       // KSPQCG       -> "qcg"
    FETIDP,    // KSPFETIDP    -> "fetidp"
    TSIRM,     // KSPTSIRM     -> "tsirm"
    SYMMLQ,    // KSPSYMMLQ    -> "symmlq"
    PYTHON,    // KSPPYTHON    -> "python"
    NONE,      // KSPNONE      -> "none"

    INVALID_SOLVER
};

/**
 * Defines an enum for preconditioner types (PC)
 *
 * Updated to include PCBDDC, PCKSP, PCPYTHON, etc.
 */
enum PreconditionerType
{
    // Existing Feel++ preconditioners (kept in the same numeric order for backward compatibility)
    IDENTITY_PRECOND = 0,   // Not a PETSc PCType, historically used in Feel++
    JACOBI_PRECOND,         // PCJACOBI ("jacobi")
    BLOCK_JACOBI_PRECOND,   // PCBJACOBI ("bjacobi")
    SOR_PRECOND,            // PCSOR ("sor")
    SSOR_PRECOND,           // Not a distinct PETSc KSP type, but historically in Feel++
    EISENSTAT_PRECOND,      // PCEISENSTAT ("eisenstat")
    ASM_PRECOND,            // PCASM ("asm")
    GASM_PRECOND,           // PCGASM ("gasm")
    CHOLESKY_PRECOND,       // PCCHOLESKY ("cholesky")
    ICC_PRECOND,            // PCICC ("icc")
    ILU_PRECOND,            // PCILU ("ilu")
    LU_PRECOND,             // PCLU ("lu")
    AMG_PRECOND,            // Not directly a single PETSc PCType (historically used)
    USER_PRECOND,           // Typically a user-defined or external approach
    SHELL_PRECOND,          // PCSHELL ("shell")
    FIELDSPLIT_PRECOND,     // PCFIELDSPLIT ("fieldsplit")
    LSC_PRECOND,            // PCLSC ("lsc")
    LSC2_PRECOND,           // Another LSC variant in Feel++, not a separate PETSc PC
    PMM_PRECOND,            // Historical Feel++ usage
    PCD_PRECOND,            // A specialized block approach in Feel++ (PCD)
    FEELPP_BLOCKNS_PRECOND, // Another custom block solver in Feel++
    FEELPP_BLOCKMS_PRECOND, // Another custom block solver in Feel++
    ML_PRECOND,             // PCML ("ml")
    GAMG_PRECOND,           // PCGAMG ("gamg")
    BOOMERAMG_PRECOND,      // PCHYPRE with BoomerAMG, or direct invocation
    AMS_PRECOND,            // PCHYPRE with AMS or "pc_hypre_type ams"
    REDUNDANT_PRECOND,      // PCREDUNDANT ("redundant")
    NONE_PRECOND,           // PCNONE ("none")

    // Recently added in previous updates
    BDDC_PRECOND,           // PCBDDC ("bddc")
    KSP_PRECOND,            // PCKSP  ("ksp")
    PYTHON_PRECOND,         // PCPYTHON ("python")

    // Append new PETSc PC types (not previously in Feel++). 
    // The numeric IDs come *after* the above to avoid renumbering existing entries.

    AMGX_PRECOND,            // PCAMGX ("amgx")
    QR_PRECOND,              // PCQR ("qr")
    NN_PRECOND,              // PCNN ("nn")
    SPAI_PRECOND,            // PCSPAI ("spai")
    MAT_PRECOND,             // PCMAT ("mat")
    HYPRE_PRECOND,           // PCHYPRE ("hypre")
    PARMS_PRECOND,           // PCPARMS ("parms")
    TFS_PRECOND,             // PCTFS ("tfs")
    GALERKIN_PRECOND,        // PCGALERKIN ("galerkin")
    EXOTIC_PRECOND,          // PCEXOTIC ("exotic")
    CP_PRECOND,              // PCCP ("cp")
    BFBT_PRECOND,            // PCBFBT ("bfbt")
    PFMG_PRECOND,            // PCPFMG ("pfmg")
    SMG_PRECOND,             // PCSMG ("smg")
    SYSPFMG_PRECOND,         // PCSYSPFMG ("syspfmg")
    REDISTRIBUTE_PRECOND,    // PCREDISTRIBUTE ("redistribute")
    SVD_PRECOND,             // PCSVD ("svd")
    CHOWILUVIENNACL_PRECOND, // PCCHOWILUVIENNACL ("chowiluviennacl")
    ROWSCALINGVIENNACL_PRECOND, // PCROWSCALINGVIENNACL ("rowscalingviennacl")
    SAVIENNACL_PRECOND,      // PCSAVIENNACL ("saviennacl")
    KACZMARZ_PRECOND,        // PCKACZMARZ ("kaczmarz")
    TELESCOPE_PRECOND,       // PCTELESCOPE ("telescope")
    PATCH_PRECOND,           // PCPATCH ("patch")
    LMVM_PRECOND,            // PCLMVM ("lmvm")
    HMG_PRECOND,             // PCHMG ("hmg")
    DEFLATION_PRECOND,       // PCDEFLATION ("deflation")
    HPDDM_PRECOND,           // PCHPDDM ("hpddm")
    H2OPUS_PRECOND,          // PCH2OPUS ("h2opus")
    MPI_PRECOND,             // PCMPI ("mpi")

    INVALID_PRECONDITIONER
};

/**
 * Defines an enum for field split type
 */
enum FieldSplitType
{
    ADDITIVE=0,
    MULTIPLICATIVE,
    SCHUR
};

/**
 * indicates the structure of the matrix vs. preconditioner
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
 * Defines an enum for iterative eigenproblem solver types
 * (SLEPc or otherwise).
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
};

/**
 * Defines an enum for eigenproblem types
 */
enum EigenProblemType
{
    NHEP=0,  // Non-Hermitian
    HEP,     // Hermitian
    GNHEP,   // Generalized non-Hermitian
    GHEP,    // Generalized Hermitian
    PGNHEP,  // Positive-definite B (GNHEP)
    INVALID_EIGENPROBLEMTYPE
};

/**
 * Defines an enum for the position of the spectrum
 */
enum PositionOfSpectrum
{
    LARGEST_MAGNITUDE=0,
    SMALLEST_MAGNITUDE,
    LARGEST_REAL,
    SMALLEST_REAL,
    LARGEST_IMAGINARY,
    SMALLEST_IMAGINARY,
#if (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR >= 9)
    TARGET_MAGNITUDE,
    TARGET_REAL,
    TARGET_IMAGINARY,
#endif
    INVALID_Postion_of_Spectrum
};

/**
 * Spectral transform type
 */
enum SpectralTransformType
{
    SHIFT=0,
    SINVERT,
    FOLD,
    CAYLEY
};

extern std::map<std::string, size_type> EigenMap;

/**
 * Defines an enum for solver packages (Feel, PETSc, etc.)
 */
enum SolverPackage
{
    SOLVERS_FEEL=0,
    SOLVERS_GMM,
    SOLVERS_EIGEN,
#if FEELPP_HAS_PETSC
    SOLVERS_PETSC,
#endif
    SOLVERS_TRILINOS,
    SOLVERS_SLEPC,
    SOLVER_INVALID_PACKAGE
};

/**
 * Define an enum for non linear solver types (SNES)
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
 * Enum for SNES line search types
 */
enum class SolverNonLinearLineSearchType
{
    BT = 0,
    NLEQERR,
    BASIC,
    L2,
    CP
};

/**
 * Aitken type
 */
enum AitkenType
{
    AITKEN_STANDARD=0,
    AITKEN_METHOD_1=1,
    FIXED_RELAXATION_METHOD=2
};

/**
 * Dirichlet enforcement type
 */
enum DirichletType
{
    STRONG=0,
    WEAK=1
};

/**
 * Projector type
 */
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
    MATSOLVER_NONE=-1,
    MATSOLVER_SPOOLES,
    MATSOLVER_SUPERLU,
    MATSOLVER_SUPERLU_DIST,
    MATSOLVER_UMFPACK,
    MATSOLVER_ESSL,
    MATSOLVER_LUSOL,
    MATSOLVER_MUMPS,
    MATSOLVER_MKL_PARDISO,
    MATSOLVER_MKL_CPARDISO,
    MATSOLVER_PASTIX,
    MATSOLVER_DSCPACK,
    MATSOLVER_MATLAB,
#if FEELPP_HAS_PETSC
    MATSOLVER_PETSC,
#endif
    MATSOLVER_PLAPACK,
    MATSOLVER_BAS,
    MATSOLVER_BOOMERAMG,
    MATSOLVER_AMS,
    MATSOLVER_EUCLID,
    MATSOLVER_PILUT
};

#if FEELPP_HAS_PETSC
#if defined(FEELPP_HAS_MUMPS) && PETSC_VERSION_GREATER_OR_EQUAL_THAN(3,2,0)
    const auto MATSOLVER_DEFAULT = MATSOLVER_MUMPS;
#else
    const auto MATSOLVER_DEFAULT = MATSOLVER_PETSC;
#endif
#else
    const auto MATSOLVER_DEFAULT = MATSOLVER_NONE;
#endif


// Forward declarations of conversion functions
PreconditionerType
pcTypeConvertStrToEnum( std::string const& type );

#if FEELPP_HAS_PETSC
KSPNormType
kspNormTypeConvertStrToEnum( std::string const& type );
#endif

SolverType
kspTypeConvertStrToEnum( std::string const& type );

SolverNonLinearType
snesTypeConvertStrToEnum( std::string const& type );
std::string
snesTypeConvertEnumToStr( SolverNonLinearType type );

SolverNonLinearLineSearchType
snesLineSearchTypeConvertStrToEnum( std::string const& type );

MatSolverPackageType
matSolverPackageConvertStrToEnum( std::string const& type );

FieldSplitType
fieldsplitTypeConvertStrToEnum( std::string const& type );

} // namespace Feel

#endif /* FEELPP_ALG_ENUMS_HPP */