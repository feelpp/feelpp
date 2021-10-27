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
#include <feel/feelalg/backend.hpp>

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
    typedef Preconditioner<T> super_type;
public:
    using super = super_type;
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
    PreconditionerPetsc( std::string const& name, worldcomm_ptr_t const& worldComm=Environment::worldCommPtr() );
    //! copy constructor
    PreconditionerPetsc( PreconditionerPetsc const & );
    //! destructor
    ~PreconditionerPetsc() override;

    /**
     * Release all memory and clear data structures.
     */
    void clear () override;

    /**
     * Initialize data structures if not done so already.
     */
    void init () override;

    /**
     * View preconditioner context
     */
    void view() const override;

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
    /*static*/ void setPetscPreconditionerType ( const PreconditionerType & preconditioner_type,
                                             const MatSolverPackageType & matSolverPackage_type,
                                             PC & pc,
                                             worldcomm_ptr_t const& worldComm=Environment::worldCommPtr(),
                                             std::string const& prefix="");

    void setPrecMatrixStructure( MatrixStructure mstruct  ) override;

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Computes the preconditioned vector "y" based on input "x".
     * Usually by solving Py=x to get the action of P^-1 x.
     */
    void apply( const Vector<T> & x, Vector<T> & y ) const override;

    void apply( Vec x, Vec y ) const;


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

    /**
     * indicate in subpreconditioner if index split has changed
     */
    bool indexSplitHasChanged() const { return M_indexSplitHasChanged; }


protected:
    void check( int err ) const { CHKERRABORT( this->worldComm().globalComm(), err ); }
private:

    bool M_indexSplitHasChanged;
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






template <typename T>
T
getOption( std::string const& name, std::string const& prefix, std::string const& sub, std::vector<std::string> const& prefixOverwrite,
           po::variables_map vm = Environment::vm() )
{
    T res = option(_name=name,_prefix=prefix,_sub=sub,_vm=vm).template as<T>();
    std::string optctx = (sub.empty())? "": sub+"-";
    for ( std::string const& prefixAdded : prefixOverwrite )
        if ( /*Environment::vm()*/vm.count( prefixvm(prefixAdded,optctx+name) ) )
            res = option(_name=name,_prefix=prefixAdded,_sub=sub,_vm=vm).template as<T>();

    return res;
}
template <typename T>
struct getOptionIfAvalaible
{
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION(
        ( std::optional<T> ), static apply, tag,
        ( required
          ( name,( std::string ) ) )
        ( optional
          ( sub,( std::string ),std::string() )
          ( prefix,( std::string ),std::string() )
          ( prefix_overload, *, std::vector<std::string>() )
          ( vm, ( po::variables_map const& ), Environment::vm() )
          ) )
#endif
    template <typename ... Ts>
    static std::optional<T> apply( Ts && ... v )
        {
            auto args = NA::make_arguments( std::forward<Ts>(v)... );
            std::string const& name = args.get(_name);
            std::string const& sub = args.get_else(_sub, "" );
            std::string const& prefix = args.get_else(_prefix, "" );
            std::vector<std::string> const& prefix_overload = args.get_else(_prefix_overload, std::vector<std::string>{} );
            po::variables_map const& vm = args.get_else(_vm, Environment::vm() );

            std::optional<T> res;
            std::string optctx = (sub.empty())? "": sub+"-";
            if ( vm.count( prefixvm(prefix,optctx+name) ) )
            {
                res = option(_name=name,_prefix=prefix,_sub=sub,_vm=vm).template as<T>();
            }
            for ( std::string const& prefixAdded : prefix_overload/*prefixOverwrite*/ )
            {
                if ( vm.count( prefixvm(prefixAdded,optctx+name) ) )
                {
                    res = option(_name=name,_prefix=prefixAdded,_sub=sub,_vm=vm).template as<T>();
                }
            }
            return res;
        }

};

template <typename T>
std::pair<bool,T>
getOptionIfAvalaibleOLD( std::string const& name, std::string const& prefix, std::string const& sub, std::vector<std::string> const& prefixOverwrite,
                         po::variables_map vm = Environment::vm() )
{
    bool hasOption=false;
    T res{};
    std::string optctx = (sub.empty())? "": sub+"-";
    if ( vm.count( prefixvm(prefix,optctx+name) ) )
    {
        hasOption = true;
        res = option(_name=name,_prefix=prefix,_sub=sub,_vm=vm).template as<T>();
    }
    for ( std::string const& prefixAdded : prefixOverwrite )
        if ( vm.count( prefixvm(prefixAdded,optctx+name) ) )
        {
            hasOption = true;
            res = option(_name=name,_prefix=prefixAdded,_sub=sub,_vm=vm).template as<T>();
        }

    return std::make_pair(hasOption,res);
}



/**
 * ConfigurePCBase
 */
struct ConfigurePCBase : public CommObject
{
public :
    using super = CommObject;
    ConfigurePCBase( PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                     std::string const& sub, std::string const& prefix )
        :
        super( worldComm ),
        M_precFeel( precFeel ),
        M_sub( sub ),
        M_prefix( prefix )
    {}

    ConfigurePCBase( PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm, std::string const& sub,
                     std::string const& prefix, std::vector<std::string> const& prefixOverwrite )
        :
        super( worldComm ),
        M_precFeel( precFeel ),
        M_sub( sub ),
        M_prefix( prefix ),
        M_prefixOverwrite( prefixOverwrite )
    {}

    ConfigurePCBase( PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm, std::string const& sub,
                     std::string const& prefix, po::options_description const& _options )
        :
        super( worldComm ),
        M_precFeel( precFeel ),
        M_sub( sub ),
        M_prefix( prefix )
    {
        this->initVariableMap( _options );
    }
    ConfigurePCBase( PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm, std::string const& sub,
                     std::string const& prefix, std::vector<std::string> const& prefixOverwrite, po::options_description const& _options )
        :
        super( worldComm ),
        M_precFeel( precFeel ),
        M_sub( sub ),
        M_prefix( prefix ),
        M_prefixOverwrite( prefixOverwrite )
    {
        this->initVariableMap( _options );
    }

    void initVariableMap( po::options_description const& _options )
    {
        M_vm.clear();

        if ( Environment::vm().count( "show-preconditioner-options" ) && Environment::isMasterRank() )
            std::cout <<_options << "\n";

        auto mycmdparser = Environment::commandLineParser();
        po::parsed_options parsed = mycmdparser.options( _options ).
            style(po::command_line_style::allow_long | po::command_line_style::long_allow_adjacent | po::command_line_style::long_allow_next).
            allow_unregistered().run();
        po::store(parsed,M_vm);
        for ( auto & configFile : Environment::configFiles() )
        {
            std::istringstream & iss = std::get<1>( configFile );
            po::store(po::parse_config_file(iss, _options,true), M_vm);
        }
        po::notify(M_vm);
    }

    PreconditionerPetsc<double> * precFeel() const { return M_precFeel; }
    std::string const& prefix() const { return M_prefix; }
    std::string const& sub() const { return M_sub; }
    bool hasPrefixOverwrite() const { return M_prefixOverwrite.size() > 0; }
    std::vector<std::string> const& prefixOverwrite() const { return M_prefixOverwrite; }
    std::string const& prefixOverwrite( int k ) const { return M_prefixOverwrite[k]; }
    po::variables_map const& vm() const { return M_vm; }

    void check( int err ) const { CHKERRABORT( this->worldComm().globalComm(), err ); }

private :

    PreconditionerPetsc<double> * M_precFeel;
    std::string M_sub, M_prefix;
    std::vector<std::string> M_prefixOverwrite;
    po::variables_map M_vm;
};


/**
 * ConfigureKSP
 */
struct ConfigureKSP : public ConfigurePCBase
{
public :
    //ConfigureKSP( KSP& ksp,WorldComm const& worldComm, std::string const& sub, std::string const& prefix );
    ConfigureKSP( KSP& ksp,PreconditionerPetsc<double> * precFeel,worldcomm_ptr_t const& worldComm, std::string const& sub, std::string const& prefix,
                  std::vector<std::string> const& prefixOverwrite,
                  std::string const& kspType = "gmres", double rtol = 1e-8, size_type maxit=1000 );
    ConfigureKSP( PreconditionerPetsc<double> * precFeel,worldcomm_ptr_t const& worldComm, std::string const& sub, std::string const& prefix,
                  std::vector<std::string> const& prefixOverwrite,
                  std::string const& kspType = "gmres", double rtol = 1e-8, size_type maxit=1000 );
    bool kspView() const { return M_kspView; }

private :

    std::string M_type;
    bool M_useConfigDefaultPetsc;
    double M_rtol;
    size_type M_maxit;
    bool M_showMonitor,M_kspView;
    bool M_constantNullSpace;
    int M_nRestartGMRES;
    //private :
public :
    void run( KSP& ksp ) const;
};


/**
 * ConfigurePC
 */
class ConfigurePC : public ConfigurePCBase
{
public :
    // ConfigurePC( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
    //              std::string const& sub = "", std::string const& prefix = "" );
    // ConfigurePC( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
    //              std::string const& sub, std::string const& prefix,
    //              std::vector<std::string> const& prefixOverwrite );
    ConfigurePC( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                 std::string const& sub, std::string const& prefix,
                 std::vector<std::string> const& prefixOverwrite,
                 po::variables_map const& vm/* = Environment::vm()*/ );

    ConfigurePC( PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                 std::string const& sub, std::string const& prefix,
                 std::vector<std::string> const& prefixOverwrite,
                 po::variables_map const& vm/* = Environment::vm()*/ );

    ConfigurePC( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                 std::string const& sub, std::string const& prefix,
                 std::vector<std::string> const& prefixOverwrite,
                 po::options_description const& _options );


    bool view() const { return M_view; }

    void setFactorShiftType( std::string s )
    {
        CHECK( s == "none" || s == "nonzero" || s == "positive_definite" || s == "inblocks" ) << "invalid shift type : " << s;
        M_factorShiftType = s;
    }

    void run( PC& pc );
private :
    bool M_useConfigDefaultPetsc, M_view;
    std::string M_factorShiftType;
};

/**
 * ConfigurePCLU
 */
class ConfigurePCLU : public ConfigurePCBase
{
public :
    ConfigurePCLU( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                   std::string const& sub, std::string const& prefix, std::vector<std::string> const& prefixOverwrite );
private :
    std::string M_matSolverPackage;
    std::vector<std::pair<bool,int> > M_mumpsParameters;
private :
    void run( PC& pc );
};

/**
 * ConfigurePCILU
 */
class ConfigurePCILU : public ConfigurePCBase
{
public :
    ConfigurePCILU( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                    std::string const& sub, std::string const& prefix, std::vector<std::string> const& prefixOverwrite );
private :
    int M_levels;
    double M_fill;
private :
    void run( PC& pc );
};


/**
 * ConfigurePCSOR
 */
class ConfigurePCSOR : public ConfigurePCBase
{
public :
    ConfigurePCSOR( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                    std::string const& sub, std::string const& prefix,std::vector<std::string> const& prefixOverwrite );
private :
    std::string M_type;
    double M_omega;
    int M_nIteration, M_nLocalIteration;
private :
    void run( PC& pc );
};

/**
 * ConfigureGASM
 */
class ConfigurePCGASM : public ConfigurePCBase
{
public :
    ConfigurePCGASM( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                     std::string const& prefix, std::vector<std::string> const& prefixOverwrite );
private :
    std::string M_type;
    int M_overlap;
private :
    void run( PC& pc );
};

/**
 * ConfigureASM
 */
class ConfigurePCASM : public ConfigurePCBase
{
public :
    ConfigurePCASM( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                    std::string const& prefix, std::vector<std::string> const& prefixOverwrite );
private :
    std::string M_type;
    int M_overlap;
private :
    void run( PC& pc );
};

/**
 * ConfigureSubPC
 */
class ConfigureSubPC : public ConfigurePCBase
{
public :
    ConfigureSubPC( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                    std::string const& prefix, std::vector<std::string> const& prefixOverwrite );
private :
    std::string M_subPCtype, M_subMatSolverPackage;
    bool M_subPCview;
    std::string M_subPCfromPCtype;
    int M_nBlock;
    //private :
public :
    void run( KSP& ksp,ConfigureKSP const& kspConf );
};

/**
 * ConfigurePCML
 */
class ConfigurePCML : public ConfigurePCBase
{
public :
    ConfigurePCML( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                   std::string const& sub, std::string const& prefix );

private :
    void run( PC& pc );
    void configurePCMLCoarse( PC& pc );

private :
    std::string M_mgType;
    int M_nLevels;
    bool M_mlReuseInterp, M_mlKeepAggInfo, M_mlReusable, M_mlOldHierarchy;

    std::string M_prefixMGCoarse;
    std::string M_coarsePCtype, M_coarsePCMatSolverPackage;
    bool M_coarsePCview;
};

/**
 * ConfigurePCGAMG
 */
class ConfigurePCGAMG : public ConfigurePCBase
{
public :
    ConfigurePCGAMG( PC& pc, PreconditionerPetsc<double> * precFeel,worldcomm_ptr_t const& worldComm,
                     std::string const& sub, std::string const& prefix );

private :
    void run( PC& pc );
    void configurePCGAMGCoarse( PC& pc );

private :
    std::string M_mgType;
    std::string M_gamgType;
    int M_nLevels;
    int M_procEqLim, M_coarseEqLim;
    double M_threshold;
    bool M_setSymGraph, M_reuseInterpolation;
    int M_nSmooths;
    bool M_coarseGridUseConfigDefaultPetsc, M_gamgLevelsUseConfigDefaultPetsc;

    std::string M_prefixMGCoarse;
    std::string M_coarsePCtype, M_coarsePCMatSolverPackage;
    bool M_coarsePCview;
};

/**
 * ConfigurePCMGLevels
 */
class ConfigurePCMGLevels : public ConfigurePCBase
{
public :
    ConfigurePCMGLevels( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                         std::string const& sub, std::string const& prefix );

private :
    void run( PC& pc, int level );

private :
    std::vector<std::string> M_prefixMGLevels;

    // get number of levels (including coarse)
    int M_nLevels;
    // ksp/pc parameters on each levels (not including coarse)
    std::vector<bool> M_mgLevelsKSPview;
    std::vector<std::string> M_mgLevelsPCtype;
    std::vector<bool> M_mgLevelsPCview;
    std::vector<std::string> M_mgLevelsMatSolverPackage;
};

/**
 * ConfigurePCFieldSplit
 */
class ConfigurePCFieldSplit : public ConfigurePCBase
{
public :
    ConfigurePCFieldSplit( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                           std::string const& sub, std::string const& prefix );

private :
    void run( PC& pc );
    void runSchur( PC& pc );


    class ConfigureSubKSP : public ConfigurePCBase
    {
    public :
        ConfigureSubKSP( KSP ** subksps/*PC& pc*/, int nSplit, PreconditionerPetsc<double> * precFeel,
                         std::string const& typeFieldSplit, worldcomm_ptr_t const& worldComm, std::string const& sub, std::string const& prefix );
    private :
        void run(KSP& ksp, int splitId );
    private :
        int M_nSplit;
        std::vector<std::string> M_prefixSplit;
        std::string M_typeFieldSplit;
        std::vector<bool> M_subPCview;
        std::vector<std::string> M_subPCtype, M_subMatSolverPackage;
    };

private :
    std::string M_type;
    std::string M_schurFactType, M_schurPrecond;
};

/**
 * ConfigurePCLSC
 */
class ConfigurePCLSC : public ConfigurePCBase
{
public :
    ConfigurePCLSC( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                    std::string const& sub, std::string const& prefix, std::string lscVersion = "petsc" );

private :
    void run( PC& pc );

private :
    std::string M_version;
    std::string M_prefixLSC;
    bool M_scaleDiag;
    std::string M_subPCtype, M_subMatSolverPackage;
    bool M_subPCview;

};

/**
 * ConfigurePCPMM
 */
class ConfigurePCPMM : public ConfigurePCBase
{
public :
    ConfigurePCPMM( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                    std::string const& sub, std::string const& prefix );

private :
    void run( PC& pc );

private :
    std::string M_prefixPMM;
    std::string M_subPCtype, M_subMatSolverPackage;
    bool M_subPCview;

};

/**
 * ConfigurePCPCD
 */
class ConfigurePCPCD : public ConfigurePCBase
{
public :
    ConfigurePCPCD( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                    std::string const& sub, std::string const& prefix );

private :
    void run( PC& pc );

private :
    std::string M_prefixPCD_Ap, M_prefixPCD_Mp;
    std::string M_subPCtype_Ap, M_subPCtype_Mp, M_subMatSolverPackage_Ap, M_subMatSolverPackage_Mp;
    bool M_subPCview_Ap, M_subPCview_Mp;

};

/**
 * ConfigurePCHYPRE_BOOMERAMG
 */
class ConfigurePCHYPRE_BOOMERAMG : public ConfigurePCBase
{
public :
    ConfigurePCHYPRE_BOOMERAMG( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                             std::string const& sub, std::string const& prefix, std::vector<std::string> const& prefixOverwrite );

private :
    void run( PC& pc );

private :
    std::optional<int> M_max_iter;
    std::optional<double> M_tol;
    std::optional<std::string> M_cycle_type;
    std::optional<int> M_max_levels;
    std::optional<std::string> M_coarsen_type;
    std::optional<double> M_strong_threshold;
    std::optional<int> M_agg_nl;
    std::optional<std::string> M_relax_type_all;
    std::optional<std::string> M_interp_type;
};

/**
 * ConfigurePCHYPRE_AMS
 */
class ConfigurePCHYPRE_AMS : public ConfigurePCBase
{
public :
    ConfigurePCHYPRE_AMS( PC& pc, PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                             std::string const& sub, std::string const& prefix, std::vector<std::string> const& prefixOverwrite );

private :
    void run( PC& pc );

private :
    std::optional<int> M_print_level;
    std::optional<int> M_max_iter;
    std::optional<int> M_cycle_type;
    std::optional<double> M_tol;
    std::optional<int> M_relax_type;
    std::optional<int> M_relax_times;
    std::optional<double> M_relax_weight;
    std::optional<double> M_omega;
    //Mat M_G;
    //Vec M_vx;
    //Vec M_vy;
    //Vec M_vz;
};


/**
 * ConfigurePCRedundant
 */
class ConfigurePCRedundant : public ConfigurePCBase
{
public :
    ConfigurePCRedundant( PreconditionerPetsc<double> * precFeel, worldcomm_ptr_t const& worldComm,
                          std::string const& prefix, std::vector<std::string> const& prefixOverwrite );
    void run( PC& pc );

    void setFactorShiftType( std::string s )
    {
        CHECK( s == "none" || s == "nonzero" || s == "positive_definite" || s == "inblocks" ) << "invalid shift type : " << s;
        M_factorShiftType = s;
    }
private :
    std::string M_innerPCtype, M_innerPCMatSolverPackage;
    bool M_innerPCview;
    std::string M_factorShiftType;
};





} // Feel
#endif /* __PreconditionerPetsc_H */
