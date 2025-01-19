/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-02-06

  Copyright (C) 2014-2016 Feel++ Consortium

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
#include <feel/feelalg/enums.hpp>

namespace Feel {

std::map<std::string, size_type> ContextOnMap = {
    {"elimination", ContextOn::ELIMINATION},
    {"penalisation", ContextOn::PENALISATION},
    {"elimination_keep_diagonal", ContextOn::ELIMINATION|ContextOn::KEEP_DIAGONAL},
    {"elimination_symmetric", ContextOn::ELIMINATION|ContextOn::SYMMETRIC},
    {"elimination_symmetric_keep_diagonal", ContextOn::ELIMINATION|ContextOn::SYMMETRIC|ContextOn::KEEP_DIAGONAL}

};


std::map<std::string, size_type> EigenMap = {
    // solver type
    {"power", POWER},
    {"laplack", LAPACK},
    {"subspace", SUBSPACE},
    {"arnoldi", ARNOLDI},
    {"krylovschur", KRYLOVSCHUR },
    {"arpack", ARPACK },
    // problem type
    { "nhep",  NHEP },
    { "hep",  HEP },
    { "gnhep",  GNHEP },
    { "ghep",  GHEP },
    { "pgnhep",  PGNHEP },
    // spectrum type
    {"largest_magnitude", LARGEST_MAGNITUDE },
    {"smallest_magnitude", SMALLEST_MAGNITUDE },
    {"largest_real", LARGEST_REAL },
    {"smallest_real", SMALLEST_REAL },
    {"largest_imaginary", LARGEST_IMAGINARY },
    {"smallest_imaginary", SMALLEST_IMAGINARY },
#if (SLEPC_VERSION_MAJOR == 3) && (SLEPC_VERSION_MINOR >= 9)
    {"target_magnitude", TARGET_MAGNITUDE },
    {"target_real", TARGET_REAL },
    {"target_imaginary", TARGET_IMAGINARY },
#endif
    // spectral transform type
    {"shift", SHIFT },
    {"shift_invert", SINVERT },
    {"fold", FOLD },
    {"cayley", CAYLEY }
};


PreconditionerType
pcTypeConvertStrToEnum(std::string const& type)
{
    // Keep the existing matches
    if      ( type == "lu" )            return PreconditionerType::LU_PRECOND;
    else if ( type == "ilu" )           return PreconditionerType::ILU_PRECOND;
    else if ( type == "id" )            return PreconditionerType::IDENTITY_PRECOND;
    else if ( type == "cholesky" )      return PreconditionerType::CHOLESKY_PRECOND;
    else if ( type == "icc" )           return PreconditionerType::ICC_PRECOND;
    else if ( type == "asm" )           return PreconditionerType::ASM_PRECOND;
    else if ( type == "gasm" )          return PreconditionerType::GASM_PRECOND;
    else if ( type == "jacobi" )        return PreconditionerType::JACOBI_PRECOND;
    else if ( type == "block_jacobi" )  return PreconditionerType::BLOCK_JACOBI_PRECOND;
    else if ( type == "bjacobi" )       return PreconditionerType::BLOCK_JACOBI_PRECOND;
    else if ( type == "sor" )           return PreconditionerType::SOR_PRECOND;
    else if ( type == "ssor" )          return PreconditionerType::SSOR_PRECOND;
    else if ( type == "eisenstat" )     return PreconditionerType::EISENSTAT_PRECOND;
    else if ( type == "shell" )         return PreconditionerType::SHELL_PRECOND;
    else if ( type == "fieldsplit" )    return PreconditionerType::FIELDSPLIT_PRECOND;
    else if ( type == "lsc" )           return PreconditionerType::LSC_PRECOND;
    else if ( type == "lsc2" )          return PreconditionerType::LSC2_PRECOND;
    else if ( type == "pmm" )           return PreconditionerType::PMM_PRECOND;
    else if ( type == "pcd" )           return PreconditionerType::PCD_PRECOND;
    else if ( type == "blockns" )       return PreconditionerType::FEELPP_BLOCKNS_PRECOND;
    else if ( type == "blockms" )       return PreconditionerType::FEELPP_BLOCKMS_PRECOND;
    else if ( type == "ml" )            return PreconditionerType::ML_PRECOND;
    else if ( type == "gamg" )          return PreconditionerType::GAMG_PRECOND;
    else if ( type == "boomeramg" )     return PreconditionerType::BOOMERAMG_PRECOND;
    else if ( type == "ams" )           return PreconditionerType::AMS_PRECOND;
    else if ( type == "redundant" )     return PreconditionerType::REDUNDANT_PRECOND;
    else if ( type == "none" )          return PreconditionerType::NONE_PRECOND;

    // Newly added in previous iteration ( Feel++ side ) 
    else if ( type == "bddc" )          return PreconditionerType::BDDC_PRECOND;
    else if ( type == "ksp" )           return PreconditionerType::KSP_PRECOND;
    else if ( type == "python" )        return PreconditionerType::PYTHON_PRECOND;

    // Additional PETSc PCs from the extended enum
    else if ( type == "amgx" )          return PreconditionerType::AMGX_PRECOND;
    else if ( type == "qr" )            return PreconditionerType::QR_PRECOND;
    else if ( type == "nn" )            return PreconditionerType::NN_PRECOND;
    else if ( type == "spai" )          return PreconditionerType::SPAI_PRECOND;
    else if ( type == "mat" )           return PreconditionerType::MAT_PRECOND;
    else if ( type == "hypre" )         return PreconditionerType::HYPRE_PRECOND;
    else if ( type == "parms" )         return PreconditionerType::PARMS_PRECOND;
    else if ( type == "tfs" )           return PreconditionerType::TFS_PRECOND;
    else if ( type == "galerkin" )      return PreconditionerType::GALERKIN_PRECOND;
    else if ( type == "exotic" )        return PreconditionerType::EXOTIC_PRECOND;
    else if ( type == "cp" )            return PreconditionerType::CP_PRECOND;
    else if ( type == "bfbt" )          return PreconditionerType::BFBT_PRECOND;
    else if ( type == "pfmg" )          return PreconditionerType::PFMG_PRECOND;
    else if ( type == "smg" )           return PreconditionerType::SMG_PRECOND;
    else if ( type == "syspfmg" )       return PreconditionerType::SYSPFMG_PRECOND;
    else if ( type == "redistribute" )  return PreconditionerType::REDISTRIBUTE_PRECOND;
    else if ( type == "svd" )           return PreconditionerType::SVD_PRECOND;
    else if ( type == "chowiluviennacl" ) return PreconditionerType::CHOWILUVIENNACL_PRECOND;
    else if ( type == "rowscalingviennacl" ) return PreconditionerType::ROWSCALINGVIENNACL_PRECOND;
    else if ( type == "saviennacl" )    return PreconditionerType::SAVIENNACL_PRECOND;
    else if ( type == "kaczmarz" )      return PreconditionerType::KACZMARZ_PRECOND;
    else if ( type == "telescope" )     return PreconditionerType::TELESCOPE_PRECOND;
    else if ( type == "patch" )         return PreconditionerType::PATCH_PRECOND;
    else if ( type == "lmvm" )          return PreconditionerType::LMVM_PRECOND;
    else if ( type == "hmg" )           return PreconditionerType::HMG_PRECOND;
    else if ( type == "deflation" )     return PreconditionerType::DEFLATION_PRECOND;
    else if ( type == "hpddm" )         return PreconditionerType::HPDDM_PRECOND;
    else if ( type == "h2opus" )        return PreconditionerType::H2OPUS_PRECOND;
    else if ( type == "mpi" )           return PreconditionerType::MPI_PRECOND;

    // Optionally map "amg" -> AMG_PRECOND or "user" -> USER_PRECOND if you want:
    else if ( type == "amg" )          return PreconditionerType::AMG_PRECOND;
    else if ( type == "user" )         return PreconditionerType::USER_PRECOND;

    // If nothing matches, default to LU
    else
        return PreconditionerType::LU_PRECOND;
}

#if FEELPP_HAS_PETSC
KSPNormType
kspNormTypeConvertStrToEnum( std::string const& type )
{
    /**/ if ( type=="default" )          return KSP_NORM_DEFAULT;
    else if ( type=="none" )             return KSP_NORM_NONE;
    else if ( type=="preconditioned" )   return KSP_NORM_PRECONDITIONED;
    else if ( type=="unpreconditioned" ) return KSP_NORM_UNPRECONDITIONED;
    else if ( type=="natural" )          return KSP_NORM_NATURAL;
    else                                 return KSP_NORM_DEFAULT;
}
#endif


SolverType
kspTypeConvertStrToEnum(std::string const& type)
{
    // --- Existing mappings from the old code ---
    if      ( type == "cg" )         return SolverType::CG;
    else if ( type == "cr" )         return SolverType::CR;
    else if ( type == "cgs" )        return SolverType::CGS;
    else if ( type == "bicg" )       return SolverType::BICG;
    else if ( type == "tcqmr" )      return SolverType::TCQMR;
    else if ( type == "tfqmr" )      return SolverType::TFQMR;
    else if ( type == "lsqr" )       return SolverType::LSQR;
    else if ( type == "bicgstab" )   return SolverType::BICGSTAB;
    else if ( type == "minres" )     return SolverType::MINRES;
    else if ( type == "gmres" )      return SolverType::GMRES;
    else if ( type == "fgmres" )     return SolverType::FGMRES;
    else if ( type == "richardson" ) return SolverType::RICHARDSON;
    else if ( type == "chebyshev" )  return SolverType::CHEBYSHEV;
    else if ( type == "preonly" )    return SolverType::PREONLY;
    else if ( type == "gcr" )        return SolverType::GCR;

    // --- Other Feel++ or PETSc solver names we haven't yet mapped ---

    else if ( type == "cgn" )        return SolverType::CGN;      // Conjugate Gradient on Normal Eqns
    else if ( type == "qmr" )        return SolverType::QMR;      // (Legacy QMR, note PETSc has KSPQMRCGS)
    else if ( type == "jacobi" )     return SolverType::JACOBI;   // (Feel++ historically used as KSP solver)
    else if ( type == "sor_forward" )return SolverType::SOR_FORWARD; 
    else if ( type == "sor_backward")return SolverType::SOR_BACKWARD;
    else if ( type == "ssor" )       return SolverType::SSOR;

    // Pipelined, deflated, flexible, etc.
    else if ( type == "dgmres" )     return SolverType::DGMRES;
    else if ( type == "lgmres" )     return SolverType::LGMRES;
    else if ( type == "pgmres" )     return SolverType::PGMRES;    // Pipelined GMRES
    else if ( type == "pipegmres" )  return SolverType::PIPEGMRES; // could map to "pipefgmres" or "pgmres" in PETSc
    else if ( type == "pipecg" )     return SolverType::PIPECG;
    else if ( type == "pipecr" )     return SolverType::PIPECR;
    else if ( type == "bcgsl" )      return SolverType::BCGSL;
    else if ( type == "fbcgs" )      return SolverType::FBCGS;
    else if ( type == "ibcgs" )      return SolverType::IBCGS;
    else if ( type == "fcg" )        return SolverType::FCG;

    // New additions from the PETSc KSP table
    else if ( type == "fbcgsr" )     return SolverType::FBCGSR;
    else if ( type == "groppcg" )    return SolverType::GROPPCG;
    else if ( type == "pipecgrr" )   return SolverType::PIPECGRR;
    else if ( type == "pipefcg" )    return SolverType::PIPEFCG;
    else if ( type == "cgls" )       return SolverType::CGLS;
    else if ( type == "nash" )       return SolverType::NASH;
    else if ( type == "stcg" )       return SolverType::STCG;
    else if ( type == "gltr" )       return SolverType::GLTR;
    else if ( type == "qcg" )        return SolverType::QCG;
    else if ( type == "fetidp" )     return SolverType::FETIDP;
    else if ( type == "tsirm" )      return SolverType::TSIRM;
    else if ( type == "symmlq" )     return SolverType::SYMMLQ;
    else if ( type == "python" )     return SolverType::PYTHON;
    else if ( type == "none" )       return SolverType::NONE;

    // If nothing matches, default to GMRES
    else
        return SolverType::GMRES;
}

SolverNonLinearType
snesTypeConvertStrToEnum( std::string const& type )
{
    /**/ if ( type == "ls" || type == "newtonls" )  return SolverNonLinearType::LINE_SEARCH;
    else if ( type == "tr" || type == "newtontr" )  return SolverNonLinearType::TRUST_REGION;
    else if ( type == "nrichardson" )               return SolverNonLinearType::NRICHARDSON;
    else if ( type == "ksponly" )                   return SolverNonLinearType::NKSPONLY;
    else if ( type == "vinewtonrsls" )              return SolverNonLinearType::VINEWTONRSLS;
    else if ( type == "vinewtonssls" )              return SolverNonLinearType::VINEWTONRSTR;
    else if ( type == "ngmres" )                    return SolverNonLinearType::NGMRES;
    else if ( type == "qn" )                        return SolverNonLinearType::QN;
    else if ( type == "shell" )                     return SolverNonLinearType::NSHELL;
    else if ( type == "gs" )                        return SolverNonLinearType::GS;
    else if ( type == "ncg" )                       return SolverNonLinearType::NCG;
    else if ( type == "fas" )                       return SolverNonLinearType::FAS;
    else if ( type == "ms" )                        return SolverNonLinearType::MS;
    else if ( type == "nasm" )                      return SolverNonLinearType::NASM;
    else if ( type == "anderson" )                  return SolverNonLinearType::ANDERSON;
    else if ( type == "aspin" )                     return SolverNonLinearType::ASPIN;
    else                                            return SolverNonLinearType::LINE_SEARCH;
}
std::string
snesTypeConvertEnumToStr( SolverNonLinearType type )
{
    /**/ if ( type == SolverNonLinearType::LINE_SEARCH )  return std::string("newtonls");
    else if ( type == SolverNonLinearType::TRUST_REGION ) return std::string("newtontr");
    else if ( type == SolverNonLinearType::NRICHARDSON )  return std::string("nrichardson");
    else if ( type == SolverNonLinearType::NKSPONLY )     return std::string("ksponly");
    else if ( type == SolverNonLinearType::VINEWTONRSLS ) return std::string("vinewtonrsls");
    else if ( type == SolverNonLinearType::VINEWTONRSTR ) return std::string("vinewtonssls");
    else if ( type == SolverNonLinearType::NGMRES )       return std::string("ngmres");
    else if ( type == SolverNonLinearType::QN )           return std::string("qn");
    else if ( type == SolverNonLinearType::NSHELL )       return std::string("shell");
    else if ( type == SolverNonLinearType::GS )           return std::string("gs");
    else if ( type == SolverNonLinearType::NCG )          return std::string("ncg");
    else if ( type == SolverNonLinearType::FAS )          return std::string("fas" );
    else if ( type == SolverNonLinearType::MS )           return std::string("ms" );
    else if ( type == SolverNonLinearType::NASM )         return std::string("nasm");
    else if ( type == SolverNonLinearType::ANDERSON )     return std::string("anderson");
    else if ( type == SolverNonLinearType::ASPIN )        return std::string("aspin" );
    else                                                  return std::string("newtonls");
}

SolverNonLinearLineSearchType
snesLineSearchTypeConvertStrToEnum( std::string const& type )
{
    /**/ if ( type == "bt" ) return SolverNonLinearLineSearchType::BT;
    else if ( type == "nleqerr" ) return SolverNonLinearLineSearchType::NLEQERR;
    else if ( type == "basic" ) return SolverNonLinearLineSearchType::BASIC;
    else if ( type == "l2" ) return SolverNonLinearLineSearchType::L2;
    else if ( type == "cp" ) return SolverNonLinearLineSearchType::CP;
    else
    {
        LOG(WARNING) << "line search type " << type << " is unknown, switch to bt";
        return SolverNonLinearLineSearchType::BT;
    }
}

MatSolverPackageType
matSolverPackageConvertStrToEnum( std::string const& type )
{
    /**/ if ( type =="spooles" )     return MatSolverPackageType::MATSOLVER_SPOOLES;
    else if ( type=="superlu" )      return MatSolverPackageType::MATSOLVER_SUPERLU;
    else if ( type=="superlu-dist" ) return MatSolverPackageType::MATSOLVER_SUPERLU_DIST;
    else if ( type=="umfpack" )      return MatSolverPackageType::MATSOLVER_UMFPACK;
    else if ( type=="essl" )         return MatSolverPackageType::MATSOLVER_ESSL;
    else if ( type=="lusol" )        return MatSolverPackageType::MATSOLVER_LUSOL;
    else if ( type=="mumps" )        return MatSolverPackageType::MATSOLVER_MUMPS;
    else if ( type=="mkl_pardiso" )  return MatSolverPackageType::MATSOLVER_MKL_PARDISO;
    else if ( type=="mkl_cpardiso" )  return MatSolverPackageType::MATSOLVER_MKL_CPARDISO;
    else if ( type=="pastix" )       return MatSolverPackageType::MATSOLVER_PASTIX;
    else if ( type=="dscpack" )      return MatSolverPackageType::MATSOLVER_DSCPACK;
    else if ( type=="matlab" )       return MatSolverPackageType::MATSOLVER_MATLAB;
#if FEELPP_HAS_PETSC
    else if ( type=="petsc" )        return MatSolverPackageType::MATSOLVER_PETSC;
#endif
    else if ( type=="plapack" )      return MatSolverPackageType::MATSOLVER_PLAPACK;
    else if ( type=="bas" )          return MatSolverPackageType::MATSOLVER_BAS;
    else if ( type=="boomeramg" )    return MatSolverPackageType::MATSOLVER_BOOMERAMG;
    else if ( type=="ams" )          return MatSolverPackageType::MATSOLVER_AMS;
    else if ( type=="euclid" )       return MatSolverPackageType::MATSOLVER_EUCLID;
    else if ( type=="pilut" )        return MatSolverPackageType::MATSOLVER_PILUT;
    else                             return MatSolverPackageType::MATSOLVER_NONE;
}

FieldSplitType
fieldsplitTypeConvertStrToEnum( std::string const& type )
{
    /**/ if ( type=="additive" )        return FieldSplitType::ADDITIVE;
    else if ( type=="multiplicative" )  return FieldSplitType::MULTIPLICATIVE;
    else if ( type=="schur" )           return FieldSplitType::SCHUR;
    else                                return FieldSplitType::ADDITIVE;
}


}
