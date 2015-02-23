/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-02

  Copyright (C) 2007-2011 Universite Joseph Fourier (Grenoble I)

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
/**
   \file solvernonlinear.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-02
 */
#ifndef __SolverNonLinear_H
#define __SolverNonLinear_H 1

#include <functional>

//#include <boost/function.hpp>
//#include <boost/bind.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/preconditioner.hpp>


namespace Feel
{

template<typename T> class Vector;
template<typename T> class MatrixSparse;


/**
 * \class SolverNonLinear
 * \brief Non linear solver base interface
 *
 * This class provides a uniform interface for nonlinear solvers.
 * This base class is overloaded to provide nonlinear solvers from
 * different packages like PETSC.
 *
 *  @author Christophe Prud'homme
 */
template <typename T>
class SolverNonLinear
{
public:


    /** @name Typedefs
     */
    //@{

    typedef SolverNonLinear<T> self_type;
    typedef boost::shared_ptr<SolverNonLinear<T> > self_ptrtype;
    typedef self_type solvernonlinear_type;
    typedef boost::shared_ptr<self_type> solvernonlinear_ptrtype;

    typedef T value_type;
    typedef typename type_traits<T>::real_type real_type;

    typedef boost::shared_ptr<Preconditioner<T> > preconditioner_ptrtype;

    typedef boost::shared_ptr<Vector<value_type> > vector_ptrtype;
    typedef boost::shared_ptr<MatrixSparse<value_type> > sparse_matrix_ptrtype;

    typedef ublas::matrix<value_type> dense_matrix_type;
    typedef ublas::vector<value_type> dense_vector_type;

    typedef Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > map_dense_matrix_type;
    typedef Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > map_dense_vector_type;

    typedef std::function<void ( const vector_ptrtype& X,
                                   vector_ptrtype& R )> residual_function_type;
    typedef std::function<void ( const vector_ptrtype& X,
                                   sparse_matrix_ptrtype& J )> jacobian_function_type;
    typedef std::function<void ( const vector_ptrtype& X,
                                   vector_ptrtype& R,
                                   sparse_matrix_ptrtype& J )> matvec_function_type;

    typedef std::function<void ( dense_vector_type const& X,
                                   dense_vector_type & R )> dense_residual_function_type;
    typedef std::function<void ( dense_vector_type const& X,
                                   dense_matrix_type& J )> dense_jacobian_function_type;
    typedef std::function<void ( dense_vector_type const& X,
                                   dense_vector_type& R,
                                   dense_matrix_type& J )> dense_matvec_function_type;

    //eigen
    typedef std::function<void ( map_dense_vector_type const& X,
                                   map_dense_vector_type & R )> map_dense_residual_function_type;
    typedef std::function<void ( map_dense_vector_type const& X,
                                   map_dense_matrix_type& J )> map_dense_jacobian_function_type;
    typedef std::function<void ( map_dense_vector_type const& X,
                                   map_dense_vector_type& R,
                                   map_dense_matrix_type& J )> map_dense_matvec_function_type;

    class NLSolveData : public boost::tuple<bool,size_type,value_type>
    {
        typedef boost::tuple<bool,size_type,value_type> super_type;
    public:
        NLSolveData() {}
        // rvalue: move constructor (fast)
        NLSolveData( super_type && i )
            :
            super_type( i )
        {}
        // copie
        NLSolveData( super_type const& i )
            :
            super_type( i )
        {}
        bool isConverged() const { return this->template get<0>(); }
        size_type nIterations() const { return this->template get<1>(); }
        value_type residual() const { return this->template get<2>(); }
    };

    // return type of solve()
    typedef NLSolveData solve_return_type;
    using SolveData = NLSolveData;
    
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     *  Constructor. Initializes Solver data structures
     */
    SolverNonLinear( std::string const& prefix = "", WorldComm const& worldComm = Environment::worldComm() );

    /**
     * copy constructor
     */
    SolverNonLinear( SolverNonLinear const & );

    /**
     * Destructor.
     */
    virtual ~SolverNonLinear();


    /**
     * Builds a \p NonlinearSolver using the nonlinear solver package specified by
     * the \p variables_map \p vm and  \p prefix
     */
    static FEELPP_DEPRECATED solvernonlinear_ptrtype build( po::variables_map const& vm, std::string const& prefix = "",
                                                            WorldComm const& worldComm = Environment::worldComm() );
    static solvernonlinear_ptrtype build( std::string const& prefix = "",
                                          WorldComm const& worldComm = Environment::worldComm() );
    static solvernonlinear_ptrtype build( std::string const& kind, std::string const& prefix = "",
                                          WorldComm const& worldComm = Environment::worldComm() );

    /**
     * Builds a \p NonlinearSolver using the nonlinear solver package specified by
     * \p solver_package
     */
    static FEELPP_DEPRECATED solvernonlinear_ptrtype build( SolverPackage solver_package, WorldComm const& worldComm = Environment::worldComm() );

    /**
     * Initialize data structures if not done so already.
     */
    virtual void init () = 0;


    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the communicator
     */
    WorldComm const& comm() const { return M_worldComm; }
    WorldComm const& worldComm() const { return M_worldComm; }

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
     * \return the preconditioner matrix structure
     * it may not be relevant to all non linear solvers
     */
    virtual MatrixStructure precMatrixStructure() const
    {
        return M_prec_matrix_structure;
    }

    SolverNonLinearType getType() const
    {
        return M_snl_type;
    }

    double getAbsoluteResidualTol() const
    {
        return M_absoluteResidualTol;
    }
    double getRelativeResidualTol() const
    {
        return M_relativeResidualTol;
    }
    double getAbsoluteSolutionTol() const
    {
        return M_absoluteSolutionTol;
    }

    uint getNbItMax() const
    {
        return M_nbItMax;
    }

    int reuseJacobian() const
    {
        return M_reuse_jac;
    }
    int reusePreconditioner() const
    {
        return M_reuse_prec;
    }

    /**
     * \return the prefix
     */
    std::string const& prefix() const{ return M_prefix; }


    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the prefix of the solver (typically for command line options)
     */
    void setPrefix( std::string const& p ) { M_prefix = p; }

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
     * Select type of non linear solver : LINEAR_SEARCH, TRUST_REGION, ...
     */
    void setType( SolverNonLinearType snl_type )
    {
        M_snl_type=snl_type;
    }
    void setType( std::string const& snl_type )
    {
        M_snl_type=snesTypeConvertStrToEnum(snl_type);
    }
    /**
     * Returns the type of solver to use.
     */
    SolverNonLinearType nlSolverType () const
    {
        return M_snl_type;
    }

    /**
     * Sets the type of solver to use.
     */
    void setKspSolverType ( const SolverType st )
    {
        M_kspSolver_type = st;
    }

    /**
     * Returns the type of solver to use.
     */
    SolverType kspSolverType () const
    {
        return M_kspSolver_type;
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
     * set reuse jacobian and/or preconditioner
     *  - jac=-1: means never rebuilt (preconditioner is not rebuit either)
     *  - jac=-2: build at each new nonlinear iterations (preconditioner is not rebuit either)
     *  - jac= 1: build at each new nonlinear iterations (rebuild jacobian every single nonlinear iteration)
     *  - jac= <n> (n>1): build every n iterations
     * when jac >= n, prec=-1 (rebuilt once at each new nonlinear iteration)
     * when jac >= n, prec=-1 (rebuilt once at each new nonlinear iteration)
     */
    virtual void setReuse( int jac=1, int prec=1 )
    {
        M_reuse_jac=jac;
        M_reuse_prec=prec;
    }

    /**
     * Define values of tolerance for the non linear solver
     */
    void setRelativeResidualTol( double tol )
    {
        M_relativeResidualTol = tol;
    }
    void setAbsoluteResidualTol( double tol )
    {
        M_absoluteResidualTol = tol;
    }
    void setAbsoluteSolutionTol( double tol )
    {
        M_absoluteSolutionTol = tol;
    }

    /**
     * Define the number max of iteration
     */
    void setNbItMax( uint n )
    {
        M_nbItMax=n;
    }
    //@}


    /**
     * show SNES monitor
     */
    bool showSNESMonitor() const { return M_showSNESMonitor; }
    void setShowSNESMonitor(bool b) { M_showSNESMonitor=b; }

    /**
     * show KSP monitor
     */
    bool showKSPMonitor() const { return M_showKSPMonitor; }
    void setShowKSPMonitor(bool b) { M_showKSPMonitor=b; }

    /**
     * show SNES converged reason
     */
    bool showSNESConvergedReason() const { return M_showSNESConvergedReason; }
    void setShowSNESConvergedReason( bool b ) { M_showSNESConvergedReason=b; }

    /**
     * show KSP converged reason
     */
    bool showKSPConvergedReason() const { return M_showKSPConvergedReason; }
    void setShowKSPConvergedReason( bool b ) { M_showKSPConvergedReason=b; }

    bool viewSNESInfo() const { return M_viewSNESInfo; }
    void setViewSNESInfo( bool b ) { M_viewSNESInfo=b; }
    /**
     * KSP relative tolerance
     */
    double rtoleranceKSP() const { return M_rtoleranceKSP; }
    void setRtoleranceKSP( double tol ) { M_rtoleranceKSP=tol; }
    /**
     * KSP divergence tolerance
     */
    double dtoleranceKSP() const { return M_dtoleranceKSP; }
    void setDtoleranceKSP( double tol ) { M_dtoleranceKSP=tol; }
    /**
     * KSP absolute tolerance
     */
    double atoleranceKSP() const { return M_atoleranceKSP; }
    void setAtoleranceKSP( double tol ) { M_atoleranceKSP=tol; }
    /**
     * KSP maximum number of iterations
     */
    size_type maxitKSP() const { return M_maxitKSP; }
    void setMaxitKSP( size_type n) { M_maxitKSP=n; }


    /** @name  Methods
     */
    //@{

    /**
     * Solves a sparse nonlinear system.
     */
    virtual solve_return_type solve ( sparse_matrix_ptrtype&,  // System Jacobian Matrix
            vector_ptrtype&, // Solution vector
            vector_ptrtype&, // Residual vector
            const double,      // Stopping tolerance
            const unsigned int ) = 0; // N. Iterations

    /**
     * Solves a dense nonlinear system.
     */
    virtual std::pair<unsigned int, real_type> solve ( dense_matrix_type&,  // System Jacobian Matrix
            dense_vector_type&, // Solution vector
            dense_vector_type&, // Residual vector
            const double,      // Stopping tolerance
            const unsigned int ) = 0; // N. Iterations

    /**
     * Solves a dense nonlinear system ( using eigen).
     */
    virtual std::pair<unsigned int, real_type> solve ( map_dense_matrix_type&,  // System Jacobian Matrix
            map_dense_vector_type&, // Solution vector
            map_dense_vector_type&, // Residual vector
            const double,      // Stopping tolerance
            const unsigned int ) = 0; // N. Iterations

    /**
     * Function that computes the residual \p R(X) of the nonlinear system
     * at the input iterate \p X.
     */
    residual_function_type residual;

    /**
     * Function that computes the Jacobian \p J(X) of the nonlinear system
     * at the input iterate \p X.
     */
    jacobian_function_type jacobian;

    /**
     * Function that computes either the residual \f$ R(X) \f$ or the
     * Jacobian \f$ J(X) \f$ of the nonlinear system at the input
     * iterate \f$ X \f$.  Note that either \p R or \p J could be
     * \p XSNULL.
     */
    matvec_function_type matvec;

    /**
     * Function that computes the residual \p R(X) of the nonlinear system
     * at the input iterate \p X.
     */
    dense_residual_function_type dense_residual;

    /**
     * Function that computes the residual \p R(X) of the nonlinear system
     * at the input iterate \p X using eigen.
     */
    map_dense_residual_function_type map_dense_residual;

    /**
     * Function that computes the Jacobian \p J(X) of the nonlinear system
     * at the input iterate \p X.
     */
    dense_jacobian_function_type dense_jacobian;

    /**
     * Function that computes the Jacobian \p J(X) of the nonlinear system
     * at the input iterate \p X using eigen.
     */
    map_dense_jacobian_function_type map_dense_jacobian;

    /**
     * Function that computes either the residual \f$ R(X) \f$ or the
     * Jacobian \f$ J(X) \f$ of the nonlinear system at the input
     * iterate \f$ X \f$.  Note that either \p R or \p J could be
     * \p XSNULL.
     */
    dense_matvec_function_type dense_matvec;

    /**
     * Function that computes either the residual \f$ R(X) \f$ or the
     * Jacobian \f$ J(X) \f$ of the nonlinear system at the input
     * iterate \f$ X \f$.  Note that either \p R or \p J could be
     * \p XSNULL.
     */
    map_dense_matvec_function_type map_dense_matvec;

    //@}



protected:

    std::string M_prefix;

    WorldComm M_worldComm;

    /**
     * Flag indicating if the data structures have been initialized.
     */
    bool M_is_initialized;

    MatrixStructure M_prec_matrix_structure;


    /**
     * Define the type of non linear solver
     */
    SolverNonLinearType M_snl_type;

    /**
     * Enum stating which type of iterative linear solver to use.
     */
    SolverType M_kspSolver_type;

    /**
     * Enum statitng with type of preconditioner to use.
     */
    PreconditionerType M_preconditioner_type;

    /**
     * Holds the Preconditioner object to be used for the linear solves.
     */
    preconditioner_ptrtype M_preconditioner;

    /**
     * Enum the software that is used to perform the factorization
     */
    MatSolverPackageType M_matSolverPackage_type;

    /**
     * Two differents tolerances on the residual for the resolution of non linear system
     */
    double M_relativeResidualTol;
    double M_absoluteResidualTol;

    /**
     * Absolute tolerances between successive iteration
     */
    double M_absoluteSolutionTol;

    /**
     * The maximum numbers of allowable nonlinear iterations
     */
    uint M_nbItMax;

    /**
     * reuse jac level
     */
    int M_reuse_jac;

    /**
     * reuse preconditioner level
     */
    int M_reuse_prec;

    bool M_showKSPMonitor, M_showSNESMonitor;
    bool M_showKSPConvergedReason, M_showSNESConvergedReason;
    bool M_viewSNESInfo;

    /**
     * KSP relative tolerance
     */
    double M_rtoleranceKSP;
    /**
     * KSP divergence tolerance
     */
    double M_dtoleranceKSP;
    /**
     * KSP absolute tolerance
     */
    double M_atoleranceKSP;
    /**
     * KSP maximum number of iterations
     */
    size_type M_maxitKSP;


};

#if !defined( FEELPP_INSTANTIATE_SOLVERNONLINEAR )
extern template class SolverNonLinear<double>;
extern template class SolverNonLinear<std::complex<double>>;
#endif

} // Feel
#endif /* __SolverNonLinear_H */
