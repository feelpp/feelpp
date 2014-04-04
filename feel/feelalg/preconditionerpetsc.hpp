/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

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
   \file preconditionerpetsc.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-01-16
 */
#ifndef __PreconditionerPetsc_H
#define __PreconditionerPetsc_H 1

#include <feel/feelcore/feelpetsc.hpp>
#include <feel/feelalg/preconditioner.hpp>

namespace Feel
{
/**
 * \class PreconditionerPetsc
 * \brief Petsc preconditioner class
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename T>
class PreconditionerPetsc
    : public Preconditioner<T>
{
public:

    typedef typename MatrixSparse<T>::indexsplit_type indexsplit_type;
    typedef typename MatrixSparse<T>::indexsplit_ptrtype indexsplit_ptrtype;


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    PreconditionerPetsc( std::string const& name, WorldComm const& worldComm=Environment::worldComm() );
    //! copy constructor
    PreconditionerPetsc( PreconditionerPetsc const & );
    //! destructor
    virtual ~PreconditionerPetsc();

    /**
     * Release all memory and clear data structures.
     */
    virtual void clear ();

    /**
     * Initialize data structures if not done so already.
     */
    virtual void init ();

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    PreconditionerPetsc& operator=( PreconditionerPetsc const & o )
    {
        if ( this != &o )
        {
            M_pc = o.M_pc;
            M_mat = o.M_mat;
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{
    /**
     * Returns the actual Petsc PC struct.  Useful for more advanced
     * purposes
     */
    PC pc()
    {
        return M_pc;
    }



    //@}

    /** @name  Mutators
     */
    //@{


    /**
     * Tells PETSC to use the user-specified preconditioner
     */
    static void setPetscPreconditionerType ( const PreconditionerType & preconditioner_type,
                                             const MatSolverPackageType & matSolverPackage_type,
                                             PC & pc,
                                             WorldComm const& worldComm=Environment::worldComm(),
                                             std::string const& prefix="");

    static void setPetscPreconditionerType ( const PreconditionerType & preconditioner_type,
                                             const MatSolverPackageType & matSolverPackage_type,
                                             PC & pc,
                                             indexsplit_ptrtype const& is,
                                             WorldComm const& worldComm=Environment::worldComm(),
                                             std::string const& prefix="");


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Computes the preconditioned vector "y" based on input "x".
     * Usually by solving Py=x to get the action of P^-1 x.
     */
    virtual void apply( const Vector<T> & x, Vector<T> & y );


    //@}

    /**
     * Preconditioner context
     */
    PC M_pc;

    /**
     * Petsc Matrix that's been pulled out of the _matrix object.
     * This happens during init...
     */
    Mat M_mat;

    static void setPetscSubpreconditionerType( PC& pc,
                                               indexsplit_ptrtype const& is,
                                               WorldComm const& worldComm=Environment::worldComm(),
                                               std::string const& prefix="" );

    static void setPetscFieldSplitPreconditionerType( PC& pc,
                                                      indexsplit_ptrtype const& is,
                                                      WorldComm const& worldComm=Environment::worldComm(),
                                                      std::string const& prefix="" );

    static void setPetscLSCPreconditionerType( PC& pc,
                                               indexsplit_ptrtype const& is,
                                               WorldComm const& worldComm=Environment::worldComm(),
                                               std::string const& prefix="" );


protected:
    void check( int err ) { CHKERRABORT( this->worldComm().globalComm(), err ); }
private:
    /**
     * Some PETSc preconditioners (ILU, LU) don't work in parallel.  This function
     * is called from setPetscPreconditionerType() to set additional options
     * for those so-called sub-preconditioners.  This method ends up being static
     * so that it can be called from set_petsc_preconditioner_type().  Not sure
     * why setPetscPreconditionerType() needs to be static though...
     */
    //#if PETSC_VERSION_LESS_THAN(3,0,0)
    ////    // In Petsc 2.3.3, PCType was #define'd as const char*
    //static void setPetscSubpreconditionerType(PCType type, PC& pc);
    //#else
    // In later versions, PCType is #define'd as char*, so we need the const

    //#endif
};


/**
 * ConfigurePCBase
 */
struct ConfigurePCBase
{
public :
    ConfigurePCBase( WorldComm const& worldComm, std::string const& sub, std::string const& prefix )
        :
        M_worldComm( worldComm ),
        M_sub( sub ),
        M_prefix( prefix )
    {}

    WorldComm const& worldComm() const { return M_worldComm; }
    std::string const& prefix() const { return M_prefix; }
    std::string const& sub() const { return M_sub; }

    void check( int err ) { CHKERRABORT( this->worldComm().globalComm(), err ); }

private :

    WorldComm const& M_worldComm;
    std::string M_sub, M_prefix;
};


struct ConfigureKSP : public ConfigurePCBase
{
public :
    ConfigureKSP( KSP& ksp,WorldComm const& worldComm, std::string const& sub, std::string const& prefix );

    bool kspView() const { return M_kspView; }

private :

    std::string M_type;
    double M_rtol;
    size_type M_maxit;
    bool M_showMonitor,M_kspView;

private :
    void runConfigureKSP( KSP& ksp );
};


/**
 * ConfigurePCML
 */
class ConfigurePCML : public ConfigurePCBase
{
public :
    ConfigurePCML( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                   WorldComm const& worldComm, std::string const& sub, std::string const& prefix );

private :
    void runConfigurePCML( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is );
    void configurePCMLCoarse( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is );

private :
    std::string M_mgType;
    int M_nLevels;
    bool M_mlReuseInterp, M_mlKeepAggInfo, M_mlReusable, M_mlOldHierarchy;

    std::string M_prefixMGCoarse;
    std::string M_coarsePCtype;//, M_coarseKSPtype;
    bool M_coarsePCview;
};

/**
 * ConfigurePCGAMG
 */
class ConfigurePCGAMG : public ConfigurePCBase
{
public :
    ConfigurePCGAMG( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                     WorldComm const& worldComm, std::string const& sub, std::string const& prefix );

private :
    void runConfigurePCGAMG( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is );
    void configurePCGAMGCoarse( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is );

private :
    std::string M_mgType;
    std::string M_gamgType;
    int M_nLevels;
    int M_procEqLim, M_coarseEqLim;
    double M_threshold;

    std::string M_prefixMGCoarse;
    std::string M_coarsePCtype;//, M_coarseKSPtype;
    bool M_coarsePCview;
};

/**
 * ConfigurePCMGLevels
 */
class ConfigurePCMGLevels : public ConfigurePCBase
{
public :
    ConfigurePCMGLevels( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
                         WorldComm const& worldComm, std::string const& sub, std::string const& prefix );

private :
    void runConfigurePCMGLevels( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is, int level );

private :
    std::vector<std::string> M_prefixMGLevels;

    // get number of levels (including coarse)
    int M_nLevels;
    // ksp parameters on each levels (not including coarse)
    //std::vector<double> M_levelksprtol;
    //std::vector<size_type> M_levelkspmaxit;
    //std::vector<bool> M_levelksponitor;
    //std::vector<std::string> M_mgLevelsKSPtype;
    std::vector<bool> M_mgLevelsKSPview;
    std::vector<std::string> M_mgLevelsPCtype;
    std::vector<bool> M_mgLevelsPCview;

};



void
configurePC( PC& pc, PreconditionerPetsc<double>::indexsplit_ptrtype const& is,
             WorldComm const& worldComm, std::string const& sub = "", std::string const& prefix = "" );




} // Feel
#endif /* __PreconditionerPetsc_H */
