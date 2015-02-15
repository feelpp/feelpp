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

#ifndef __linear_solver_h__
#define __linear_solver_h__

#include <boost/mpi/communicator.hpp>

#include <feel/feelcore/parameter.hpp>
#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelalg/nullspace.hpp>
#include <feel/feelcore/traits.hpp>

namespace Feel
{
namespace mpi=boost::mpi;
template<typename T> class Vector;
template<typename T> class MatrixSparse;

/**
 * This class provides a uniform interface for linear solvers.  This base
 * class is overloaded to provide linear solvers from different packages
 * like FEEL, GMM or PETSC
 *
 * @author Benjamin Kirk, 2003
 * @author Christophe Prud'homme, 2005
 */
template <typename T>
class SolverLinear
{

public:

    typedef SolverLinear<T> self_type;
    typedef boost::shared_ptr<SolverLinear<T> >  self_ptrtype;

    typedef T value_type;
    typedef typename type_traits<T>::real_type real_type;

    typedef boost::shared_ptr<Preconditioner<T> > preconditioner_ptrtype;

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
    SolverLinear ( WorldComm const& worldComm = Environment::worldComm() );

    /**
     *  Constructor. Initializes Solver data structures
     */
    SolverLinear ( po::variables_map const& vm, WorldComm const& worldComm = Environment::worldComm() );

    /**
     * Destructor.
     */
    virtual ~SolverLinear ();

    WorldComm const& worldComm() const
    {
        return M_worldComm;
    }
    void setWorldComm( WorldComm const& worldComm )
    {
        M_worldComm=worldComm;
    }

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
    po::variables_map vm() const
    {
        return M_vm;
    }

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
    BOOST_PARAMETER_MEMBER_FUNCTION( ( void ),
                                     setTolerances,
                                     tag,
                                     ( required
                                       ( rtolerance,( double ) )
                                     )
                                     ( optional
                                       ( maxit,( size_type ), 1000 )
                                       ( atolerance,( double ), 1e-50 )
                                       ( dtolerance,( double ), 1e5 )
                                     ) )
    {
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

    void attachNullSpace( boost::shared_ptr<NullSpace<value_type> > const& ns )
    {
        M_nullSpace = ns;
    }
    void attachNearNullSpace( boost::shared_ptr<NullSpace<value_type> > const& ns )
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
            Vector<T>& x,
            Vector<T> const& b,
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
            Vector<T>& x,
            Vector<T> const& b,
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

private:

    //mpi communicator
    WorldComm M_worldComm;

protected:

    ///
    po::variables_map M_vm;

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
    boost::shared_ptr<NullSpace<value_type> > M_nullSpace, M_nearNullSpace;

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
};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
SolverLinear<T>::SolverLinear ( WorldComm const& worldComm ) :
    M_worldComm( worldComm ),
    M_solver_type         ( GMRES ),
    M_preconditioner_type ( LU_PRECOND ),
    M_preconditioner(),
    M_is_initialized      ( false ),
    M_prec_matrix_structure( SAME_NONZERO_PATTERN ),
    M_showKSPMonitor( false ),
    M_showKSPConvergedReason( false )
{
}

template <typename T>
inline
SolverLinear<T>::SolverLinear ( po::variables_map const& vm, WorldComm const& worldComm ) :
    M_worldComm( worldComm ),
    M_vm( vm ),
    M_solver_type         ( GMRES ),
    M_preconditioner_type ( LU_PRECOND ),
    M_is_initialized      ( false ),
    M_prec_matrix_structure( SAME_NONZERO_PATTERN ),
    M_showKSPMonitor( false ),
    M_showKSPConvergedReason( false )
{
}



template <typename T>
inline
SolverLinear<T>::~SolverLinear ()
{
    this->clear ();
}


} // Feel

#endif // #ifdef __solver_h__
