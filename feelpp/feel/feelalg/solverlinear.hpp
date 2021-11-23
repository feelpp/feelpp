// $Id: linear_solver.h,v 1.2 2005/02/22 22:17:34 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef FEELPP_ALG_SOLVERLINEAR_H
#define FEELPP_ALG_SOLVERLINEAR_H

#include <boost/mpi/communicator.hpp>

#include <feel/feelcore/parameter.hpp>
#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelalg/nullspace.hpp>
#include <feel/feelcore/traits.hpp>

namespace Feel
{
namespace mpi=boost::mpi;
template<typename T, typename SizeT> class Vector;
template<typename T> class MatrixSparse;

/**
 * This class provides a uniform interface for linear solvers.  This base
 * class is overloaded to provide linear solvers from different packages
 * like FEEL, GMM or PETSC
 *
 * @author Benjamin Kirk, 2003
 * @author Christophe Prud'homme, 2005
 */
template <typename T,typename SizeT = uint32_type>
class SolverLinear : public CommObject
{

public:

    using super = CommObject;
    typedef SolverLinear<T> self_type;
    typedef std::shared_ptr<SolverLinear<T> >  self_ptrtype;
    using size_type = SizeT;
    typedef T value_type;
    typedef typename type_traits<T>::real_type real_type;

    typedef std::shared_ptr<Preconditioner<T> > preconditioner_ptrtype;

    class SolveData : public boost::tuple<bool,size_type,value_type>
    {
        typedef boost::tuple<bool,size_type,value_type> super_type;
    public:
        SolveData() {}
        // rvalue: move constructor (fast)
        SolveData( super_type && i )
            :
            super_type( i )
        {}
        // copie
        SolveData( super_type const& i )
            :
            super_type( i )
        {}
        bool isConverged() const { return this->template get<0>(); }
        size_type nIterations() const { return this->template get<1>(); }
        value_type residual() const { return this->template get<2>(); }
    };

    // return type of solve()
    typedef SolveData solve_return_type;

    /**
     *  Constructor. Initializes Solver data structures
     */
    explicit SolverLinear ( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

    /**
     *  Constructor. Initializes Solver data structures
     */
    SolverLinear ( std::string const& prefix, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(), po::variables_map const& vm = Environment::vm() );

    /**
     * Destructor.
     */
    ~SolverLinear () override;

    /**
     * @returns true if the data structures are
     * initialized, false otherwise.
     */
    bool initialized () const
    {
        return M_is_initialized;
    }


    /**
     * Release all memory and clear data structures.
     */
    virtual void clear () {}

    /**
     * Initialize data structures if not done so already.
     */
    virtual void init () = 0;

    /**
     * return variables_map
     */
    // po::variables_map vm() const
    // {
    //     return M_vm;
    // }

    /**
     * \return the relative tolerance
     */
    value_type rTolerance() const
    {
        return M_rtolerance;
    }

    /**
     * \return the divergence tolerance
     */
    value_type dTolerance() const
    {
        return M_dtolerance;
    }

    /**
     * \return the absolute tolerance
     */
    value_type aTolerance() const
    {
        return M_atolerance;
    }

    /**
     * Returns the type of solver to use.
     */
    SolverType solverType () const
    {
        return M_solver_type;
    }

    /**
     * \return the maximum number of iterations
     */
    size_type maxIterations() const
    {
        return M_maxit;
    }

    /**
     * \return the prefix
     */
    std::string const& prefix() const{ return M_prefix; }

    /**
     * set the prefix of the solver (typically for command line options)
     */
    void setPrefix( std::string const& p ) { M_prefix = p; }

    /**
     * set tolerances: relative tolerance \p rtol, divergence tolerance \p dtol
     * and absolute tolerance \p atol
     */
    template <typename ... Ts>
    void setTolerances( Ts && ... v )
        {
            auto args = NA::make_arguments( std::forward<Ts>(v)... );
            double rtolerance = args.get(_rtolerance );
            size_type maxit = args.get_else(_maxit,1000);
            double atolerance = args.get_else(_atolerance,1e-50 );
            double dtolerance = args.get_else(_dtolerance,1e5 );
            M_rtolerance = rtolerance;
            M_dtolerance = dtolerance;
            M_atolerance = atolerance;
            M_maxit=maxit;
        }

    /**
     * Sets the type of solver to use.
     */
    void setSolverType ( const SolverType st )
    {
        M_solver_type = st;
    }

    /**
     * Returns the type of preconditioner to use.
     */
    PreconditionerType preconditionerType () const
    {
        if ( M_preconditioner )
            return M_preconditioner->type();

        return M_preconditioner_type;
    }

    /**
     * Sets the type of preconditioner to use.
     */
    void setPreconditionerType ( const PreconditionerType pct )
    {
        if ( M_preconditioner )
            M_preconditioner->setType( pct );

        else
            M_preconditioner_type = pct;
    }

    /**
     * Attaches a Preconditioner object to be used by the solver
     */
    void attachPreconditioner( preconditioner_ptrtype preconditioner )
    {
        if ( this->M_is_initialized )
        {
            std::cerr<<"Preconditioner must be attached before the solver is initialized!"<<std::endl;
        }

        M_preconditioner_type = SHELL_PRECOND;
        M_preconditioner = preconditioner;
    }

    void attachNullSpace( std::shared_ptr<NullSpace<value_type> > const& ns )
    {
        M_nullSpace = ns;
    }
    void attachNearNullSpace( std::shared_ptr<NullSpace<value_type> > const& ns )
    {
        M_nearNullSpace = ns;
    }

    void setFieldSplitType( const FieldSplitType fst )
    {
        M_fieldSplit_type = fst;
    }

    FieldSplitType fieldSplitType() const
    {
        return M_fieldSplit_type;
    }

    /**
     * Sets the type of preconditioner to use.
     */
    void setMatSolverPackageType ( const MatSolverPackageType mspackt )
    {
        M_matSolverPackage_type = mspackt;
    }
    /**
     * Returns the type of preconditioner to use.
     */
    MatSolverPackageType matSolverPackageType () const
    {
        return M_matSolverPackage_type;
    }

    /**
     * \return the preconditioner matrix structure
     * it may not be relevant to all non linear solvers
     */
    virtual MatrixStructure precMatrixStructure() const
    {
        return M_prec_matrix_structure;
    }

    /**
     * \return the preconditioner matrix structure
     * it may not be relevant to all non linear solvers
     */
    virtual void setPrecMatrixStructure( MatrixStructure mstruct  )
    {
        // warning : in boths cases!
        if ( M_preconditioner )
            M_preconditioner->setPrecMatrixStructure(mstruct);

        M_prec_matrix_structure = mstruct;
    }

    /**
     * show KSP monitor
     */
    bool showKSPMonitor() const { return M_showKSPMonitor; }
    void setShowKSPMonitor( bool b ) { M_showKSPMonitor=b; }

    /**
     * show KSP converged reason
     */
    bool showKSPConvergedReason() const { return M_showKSPConvergedReason; }
    void setShowKSPConvergedReason( bool b ) { M_showKSPConvergedReason=b; }

    /**
     * This function calls the solver "M_solver_type" preconditioned
     * with the "M_preconditioner_type" preconditioner.  Note that
     * this method will compute the preconditioner from the system
     * matrix.
     *
     * \param mat System Matrix
     * \param prec Preconditioning Matrix
     * \param x Solution vector
     * \param b RHS vector
     * \param tolerance Stopping tolerance
     * \param maxit maximum Number of Iterations
     * \param transpose true to solve the transpose system, false otherwise
     */
    virtual
    solve_return_type
    solve ( MatrixSparse<T> const& mat,
            Vector<T,size_type>& x,
            Vector<T,size_type> const& b,
            const double tolerance,
            const unsigned int maxit,
            bool transpose
          ) = 0;



    /**
     * This function calls the solver
     * "M_solver_type" preconditioned with the
     * "M_preconditioner_type" preconditioner.  Note that this method
     * will compute the preconditioner from the system matrix.
     *
     * \param mat System Matrix
     * \param prec Preconditioning Matrix
     * \param x Solution vector
     * \param b RHS vector
     * \param tolerance Stopping tolerance
     * \param maxit maximum Number of Iterations
     * \param transpose true to solve the transpose system, false otherwise
     */
    virtual
    solve_return_type
    solve ( MatrixSparse<T> const& mat,
            MatrixSparse<T> const& prec,
            Vector<T,size_type>& x,
            Vector<T,size_type> const& b,
            const double tolerance,
            const unsigned int maxit,
            bool transpose
          ) = 0;


protected:

    /**
     * set initialized only for subclasses
     */
    void setInitialized( bool init )
    {
        M_is_initialized = init;
    }

protected:

    ///
    // po::variables_map M_vm;

    std::string M_prefix;

    /// relative tolerance
    double M_rtolerance;

    /// divergence tolerance
    double M_dtolerance;

    /// absolute tolerance
    double M_atolerance;

    /// maximum number of iterations
    size_type M_maxit;

    /**
     * Enum stating which type of iterative solver to use.
     */
    SolverType M_solver_type;

    /**
     * Enum statitng with type of preconditioner to use.
     */
    PreconditionerType M_preconditioner_type;

    /**
     * Holds the Preconditioner object to be used for the linear solves.
     */
    preconditioner_ptrtype M_preconditioner;

    /**
     * Near Null Space
     */
    std::shared_ptr<NullSpace<value_type> > M_nullSpace, M_nearNullSpace;

    FieldSplitType M_fieldSplit_type;

    /**
     * Enum the software that is used to perform the factorization
     */
    MatSolverPackageType M_matSolverPackage_type;


    /**
     * Flag indicating if the data structures have been initialized.
     */
    bool M_is_initialized;

    MatrixStructure M_prec_matrix_structure;

    bool M_showKSPMonitor;
    bool M_showKSPConvergedReason;

    bool M_kspView;
    bool M_kspUseInitialGuessNonZero;
    std::string M_kspNormType;
    int M_nRestartGMRES, M_nRestartFGMRES, M_nRestartGCR;
};




/*----------------------- inline functions ----------------------------------*/
template <typename T, typename SizeT>
inline
SolverLinear<T,SizeT>::SolverLinear ( worldcomm_ptr_t const& worldComm ) :
    super( worldComm ),
    M_rtolerance( 1e-8 ),
    M_dtolerance( 1e5 ),
    M_atolerance( 1e-50 ),
    M_maxit( 1000 ),
    M_solver_type         ( GMRES ),
    M_preconditioner_type ( LU_PRECOND ),
    M_preconditioner(),
    M_is_initialized      ( false ),
    M_prec_matrix_structure( SAME_NONZERO_PATTERN ),
    M_showKSPMonitor( false ),
    M_showKSPConvergedReason( false ),
    M_kspView( false ), M_kspUseInitialGuessNonZero( false ),
    M_kspNormType( "default" ),
    M_nRestartGMRES( 30 ), M_nRestartFGMRES( 30 ), M_nRestartGCR( 30 )
{
}

template <typename T, typename SizeT>
inline
SolverLinear<T,SizeT>::SolverLinear ( std::string const& prefix, worldcomm_ptr_t const& worldComm, po::variables_map const& vm ) :
    super( worldComm ),
    //M_vm( vm ),
    M_prefix( prefix ),
    M_rtolerance( 1e-8 ),
    M_dtolerance( 1e5 ),
    M_atolerance( 1e-50 ),
    M_maxit( 1000 ),
    M_solver_type         ( GMRES ),
    M_preconditioner_type ( LU_PRECOND ),
    M_is_initialized      ( false ),
    M_prec_matrix_structure( SAME_NONZERO_PATTERN ),
    M_showKSPMonitor( false ),
    M_showKSPConvergedReason( false ),
    M_kspView( boption(_name="ksp-view",_prefix=prefix,_vm=vm) ),
    M_kspUseInitialGuessNonZero( boption(_name="ksp-use-initial-guess-nonzero", _prefix=prefix,_vm=vm) ),
    M_kspNormType( soption(_name="ksp-norm-type",_prefix=prefix,_vm=vm) ),
    M_nRestartGMRES( ioption(_name="gmres-restart", _prefix=prefix,_vm=vm ) ),
    M_nRestartFGMRES( ioption(_name="fgmres-restart", _prefix=prefix,_vm=vm ) ),
    M_nRestartGCR( ioption(_name="gcr-restart", _prefix=prefix,_vm=vm ) )
{
}



template <typename T, typename SizeT>
inline
SolverLinear<T,SizeT>::~SolverLinear ()
{
    this->clear ();
}


} // Feel

#endif // #ifdef __solver_h__
