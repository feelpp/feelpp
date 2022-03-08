/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-04

  Copyright (C) 2007-2011 Université Joseph Fourier (Grenoble I)

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
   \file solvereigen.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-04
 */
#ifndef FEELPP_ALG_SOLVEREIGEN_H
#define FEELPP_ALG_SOLVEREIGEN_H 1

#include <boost/tuple/tuple.hpp>
#include <feel/feelcore/parameter.hpp>

#include <feel/feelalg/enums.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/matrixsparse.hpp>

namespace Feel
{
/**
 * \class SolverEigen
 * \brief base class for eigen solvers
 *
 *  @author Christophe Prud'homme
 */
template<typename T>
class SolverEigen
{
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef typename type_traits<T>::real_type real_type;

    typedef SolverEigen<value_type> solvereigen_type;
    typedef std::shared_ptr<solvereigen_type> solvereigen_ptrtype;
    typedef DataMap<> datamap_type;
    typedef std::shared_ptr<datamap_type> datamap_ptrtype;
    using size_type = typename datamap_type::size_type;
    typedef boost::tuple<size_type, size_type, std::vector<double> > solve_return_type;

    typedef Vector<value_type> vector_type;
    typedef std::shared_ptr<vector_type> vector_ptrtype;
    typedef MatrixSparse<value_type> sparse_matrix_type;
    typedef std::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;

    typedef boost::tuple<real_type, real_type, vector_ptrtype> eigenpair_type;
    typedef std::map<real_type, eigenpair_type> eigenmodes_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     *  Constructor. Initializes Solver data structures
     */
    SolverEigen();

    /**
     *  Constructor. Initializes Solver data structures
     *
     * The \p prefix parameter allows to set different eigensolver options for
     * different eigensolver. It allows to distinguish between these options
     * \code
     * // register two slepc eigensolver options
     * add_options( solvereigen_options( "eigen1" ) ).add_options( solvereigen_options( "eigen2" ));
     * // build an eigen solver associated with option set eigen1
     * SolverEigen<double>::build( vm, "eigen1" );
     * \endcode
     *
     * \param vm variables map
     * \param prefix string that allows for various options of the same type
     */
    SolverEigen( po::variables_map const& vm, std::string const& prefix = "" );

    /**
     * copy constructor
     */
    SolverEigen( SolverEigen const & );

    /**
     * Destructor.
     */
    virtual ~SolverEigen();

    /**
     * Builds a \p SolverEigen using the linear solver package
     * specified by \p solver_package
     */
    static std::shared_ptr<SolverEigen<value_type> > build( const SolverPackage solver_package = SOLVERS_SLEPC );

    /**
     * Builds a \p SolverEigen using the linear solver package
     * specified by \p vm
     * \param vm variables_map that contains the command line options and their defaults
     * \param prefix string that allows for various options of the same type
     */
    static std::shared_ptr<SolverEigen<value_type> > build( po::variables_map const& vm,
            std::string const& prefix = std::string() );

    /**
     * Returns the \p ith eigenvalue (real and imaginary part), and
     * copies the \ ith eigen vector to the solution vector.
     */
    virtual eigenpair_type eigenPair ( unsigned int i ) = 0;

    /*
     * Returns the eigen modes in a map
     */
    virtual eigenmodes_type eigenModes () = 0;

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the prefix that differentiates the eigen solvers
     */
    std::string const& prefix() const
    {
        return M_prefix;
    }

    /**
     * Returns the type of eigensolver to use.
     */
    EigenSolverType eigenSolverType () const
    {
        return M_eigen_solver_type;
    }

    /**
     * Returns the type of the eigen problem.
     */
    EigenProblemType eigenProblemType () const
    {
        return M_eigen_problem_type;
    }

    /**
     * Returns the position of the spectrum to compute.
     */
    PositionOfSpectrum postitionOfSpectrum () const
    {
        return M_position_of_spectrum;
    }

    /**
     * Returns the spectral transforms
     */
    SpectralTransformType spectralTransform() const
    {
        return M_spectral_transform;
    }

    /**
     * Returns the number of eigenvalues to compute
     */
    size_type numberOfEigenvalues() const
    {
        return M_nev;
    }

    /**
     * Returns the dimension of the subspace
     */
    size_type numberOfEigenvaluesConverged() const
    {
        return M_ncv;
    }

    /**
     * Returns the maximum projected dimension
     */
    size_type maximumProjectedDimension() const
    {
        return M_mpd;
    }

    /**
     * Returns the start of the interval
     */
    double intervalA() const
    {
        return M_interval_a;
    }

    /**
     * Returns the end of the interval
     */
    double intervalB() const
        {
            return M_interval_b;
        }

    /**
     * Returns the tolerance to be reached by eigenvalue solver
     */
    value_type tolerance() const
    {
        return M_tolerance;
    }

    /**
     * Returns the maximum number of iterations
     */
    size_type maxIterations() const
    {
        return M_maxit;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    /**
     * Sets the type of eigensolver to use.
     */
    void setEigenSolverType ( const EigenSolverType est )
    {
        M_eigen_solver_type = est;
    }

    /**
     * Sets the type of the eigenproblem.
     */
    void setEigenProblemType ( EigenProblemType ept )
    {
        M_eigen_problem_type = ept;
    }

    /**
     * Sets the position of the spectrum.
     */
    void setPositionOfSpectrum ( PositionOfSpectrum pos )
    {
        M_position_of_spectrum= pos;
    }

    /**
     * set the spectral transform
     */
    void setSpectralTransform( SpectralTransformType st )
    {
        M_spectral_transform = st;
    }

    /**
     * set the tolerance
     * \param tol tolerance to reach
     */
    void setTolerance( value_type tol )
    {
        M_tolerance = tol;
    }

    /**
     * set the max number of iterations
     * \param maxiter maximum number of iterations
     */
    void setMaxIterations( size_type maxiter )
    {
        M_maxit = maxiter;
    }

    /**
     * set the number of eigenvalues to compute
     */
    void setNumberOfEigenValues( size_type nev )
    {
        M_nev = nev;
    }

    /**
     * set the dimension of the subspace
     */
    void setNumberOfEigenValuesConverged( size_type ncv )
    {
        M_ncv = ncv;
    }

    /**
     * set the maximum projected dimension
     */
    void setMaximumProjectedDimension( size_type mpd )
    {
        M_mpd = mpd;
    }

    /**
     * set the interval
     */
    void setInterval( double interval_a, double interval_b )
    {
        M_interval_a = interval_a;
        M_interval_b = interval_b;
    }

    datamap_type const& mapRow() const
    {
        return *M_mapRow;
    }
    datamap_type const& mapCol() const
    {
        return *M_mapCol;
    }
    datamap_ptrtype const& mapRowPtr() const
    {
        return M_mapRow;
    }
    datamap_ptrtype const& mapColPtr() const
    {
        return M_mapCol;
    }

    void setMapRow( datamap_ptrtype const& d )
    {
        M_mapRow=d;
    }
    void setMapCol( datamap_ptrtype const& d )
    {
        M_mapCol=d;
    }

    //@}

    /** @name  Methods
     */
    //@{

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


    template <typename ... Ts,typename  = typename std::enable_if_t< sizeof...(Ts) != 0 && ( NA::is_named_argument_v<Ts> && ...) > >
    solve_return_type solve( Ts && ... v )
        {
            auto args = NA::make_arguments( std::forward<Ts>(v)... );
            sparse_matrix_ptrtype matrixA = args.get(_matrixA);
            sparse_matrix_ptrtype matrixB = args.get(_matrixB);
            size_type maxit = args.get_else(_maxit,1000);
            double tolerance = args.get_else(_tolerance,1e-11 );
            this->setMaxIterations( maxit );
            this->setTolerance( tolerance );
            solve_return_type ret =  solve( matrixA, matrixB, M_nev, M_ncv, this->tolerance(), this->maxIterations() );
            return ret;
        }


    /**
     * Solves the standard eigen problem and returns the
     * number of converged eigenpairs and the number
     * of iterations.
     */
    virtual solve_return_type solve ( MatrixSparse<value_type> &matrix_A,
                                      int nev,
                                      int ncv,
                                      const double tol,
                                      const unsigned int m_its ) = 0;

    /**
     * Solves the standard eigen problem and returns the
     * number of converged eigenpairs and the number
     * of iterations.
     */
    solve_return_type solve ( std::shared_ptr<MatrixSparse<value_type> > &matrix_A,
                              int nev,
                              int ncv,
                              const double tol,
                              const unsigned int m_its )
    {
        return this->solve( *matrix_A, nev, ncv, tol, m_its );
    }


    /**
     * Solves the generalized eigen value problem \f$A x = \lambda
     * Bx\f$ and returns the number of converged eigenpairs and the
     * number of iterations.
     */
    virtual solve_return_type solve ( MatrixSparse<value_type> &matrix_A,
                                      MatrixSparse<value_type> &matrix_B,
                                      int nev,
                                      int ncv,
                                      const double tol,
                                      const unsigned int m_its ) = 0;

    /**
     * Solves the generalized eigen value problem \f$A x = \lambda
     * Bx\f$ and returns the number of converged eigenpairs and the
     * number of iterations.
     */
    solve_return_type solve ( std::shared_ptr<MatrixSparse<value_type> > &matrix_A,
                              std::shared_ptr<MatrixSparse<value_type> > &matrix_B,
                              int nev,
                              int ncv,
                              const double tol,
                              const unsigned int m_its )
    {
        return this->solve( *matrix_A, *matrix_B, nev, ncv, tol, m_its );
    }

    //@}





protected:

    /** prefix to differentiate eigen solvers */
    std::string M_prefix;

    /**
     * Enum stating which type of eigensolver to use.
     */
    EigenSolverType M_eigen_solver_type;

    /**
     * Enum stating which type of eigen problem we deal with.
     */
    EigenProblemType M_eigen_problem_type;

    /**
     * Enum stating where to evaluate the spectrum.
     */
    PositionOfSpectrum M_position_of_spectrum;

    /**
     * spectral transformation type
     */
    SpectralTransformType M_spectral_transform;

    /**
     * Flag indicating if the data structures have been initialized.
     */
    bool M_is_initialized;

    //! number of eigenvalues
    size_type M_nev;

    //! dimension of the subspace
    size_type M_ncv;

    //! maximum projected dimension
    size_type M_mpd;

    //! start of the interval
    double M_interval_a;
    //! end of the interval
    double M_interval_b;

    //! max number of iterations
    size_type M_maxit;

    //! tolerance
    value_type M_tolerance;

    // datamap
    datamap_ptrtype M_mapRow,M_mapCol;

};

/**
 * \fn eigs
 * \brief solve for \f$A x =\lambda B x\f$
 *
 * \param A sparse matrix
 * \param B sparse matrix
 * \param ncv number of converged eigen values (dimension of Krylov space)
 * \param nev number of eigen values
 * \param backend type of backend (default BACKEND_PETSC)
 * \param solver type of solver (default KRYLOVSCHUR)
 * \param problem type of problem (default GHEP)
 * \param transform type of spectral transform (default SHIFT)
 * \param spectrum type of eigenvalue (default LARGEST_MAGNITUDE)
 * \param maxit number of maximum iterations
 * \param tolerance tolerance for the eigenvalue solver (default 1e-11)
 *
 * \return eigen modes
 */
template <typename ... Ts>
typename SolverEigen<double>::eigenmodes_type eigs( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    d_sparse_matrix_ptrtype matrixA = args.get(_matrixA);
    d_sparse_matrix_ptrtype matrixB = args.get(_matrixB);
    int nev = args.get_else_invocable(_nev,[](){ return ioption(_name="solvereigen.nev"); } );
    int ncv = args.get_else_invocable(_ncv,[](){ return ioption(_name="solvereigen.ncv"); } );
    int mpd = args.get_else_invocable(_mpd,[](){ return ioption(_name="solvereigen.mpd"); } );
    double interval_a = args.get_else_invocable(_interval_a,[](){ return doption(_name="solvereigen.interval-a"); } );
    double interval_b = args.get_else_invocable(_interval_b,[](){ return doption(_name="solvereigen.interval-b"); } );
    BackendType backend = args.get_else(_backend,BACKEND_DEFAULT);
    EigenSolverType solver = args.get_else(_solver,KRYLOVSCHUR);
    EigenProblemType problem = args.get_else(_problem,GHEP);
    SpectralTransformType transform = args.get_else(_transform,SHIFT);
    PositionOfSpectrum spectrum = args.get_else(_spectrum,LARGEST_MAGNITUDE);
    int maxit = args.get_else_invocable(_maxit,[](){ return ioption(_name="solvereigen.maxiter"); } );
    double tolerance = args.get_else_invocable(_tolerance,[](){ return doption(_name="solvereigen.tolerance"); } );
    bool verbose = args.get_else_invocable(_verbose,[](){ return boption(_name="solvereigen.verbose"); } );

    typedef std::shared_ptr<Vector<double> > vector_ptrtype;
    //std::shared_ptr<SolverEigen<double> > eigen = SolverEigen<double>::build(  backend );
    std::shared_ptr<SolverEigen<double> > eigen = SolverEigen<double>::build();
    eigen->setEigenSolverType( solver );
    eigen->setEigenProblemType( problem );
    eigen->setPositionOfSpectrum( spectrum );
    eigen->setNumberOfEigenValues( nev );
    eigen->setNumberOfEigenValuesConverged( ncv );
    eigen->setMaximumProjectedDimension( (mpd>0)?mpd:invalid_v<uint32_type> );
    eigen->setInterval( interval_a, interval_b );
    eigen->setMaxIterations( maxit );
    eigen->setSpectralTransform( transform );
    eigen->setTolerance( tolerance );
    eigen->setMapRow( matrixA->mapRowPtr() );
    eigen->setMapCol( matrixA->mapColPtr() );

    LOG(INFO) << "number of eigen values = " << nev << "\n";
    LOG(INFO) << "number of eigen values converged = " << ncv << "\n";
    LOG(INFO) << "number of eigen value solver iterations = " << maxit << "\n";
    LOG(INFO) << "eigenvalue tolerance = " << tolerance << "\n";

    unsigned int nconv, nits;
    std::vector<double> err( ncv );
    boost::tie( nconv, nits, err )  = eigen->solve( _matrixA=matrixA,
                                                    _matrixB=matrixB,
                                                    _maxit=maxit,
                                                    _tolerance=tolerance );

    if ( verbose )
    {
        int i = 0;
        for( auto e : err )
        {
            if ( Environment::worldComm().isMasterRank() )
            {
                std::cout << i++ << "-th mode error ||A x - lambda B x ||/||x|| = " << e << "\n";
            }
        }
    }

    return eigen->eigenModes();
}

template <typename ... Ts>
auto veigs( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && formA = args.get(_formA);
    auto && formB = args.get(_formB);
    int nev = args.get_else_invocable(_nev,[](){ return ioption(_name="solvereigen.nev"); } );
    int ncv = args.get_else_invocable(_ncv,[](){ return ioption(_name="solvereigen.ncv"); } );
    int mpd = args.get_else_invocable(_mpd,[](){ return ioption(_name="solvereigen.mpd"); } );
    double interval_a = args.get_else_invocable(_interval_a,[](){ return doption(_name="solvereigen.interval-a"); } );
    double interval_b = args.get_else_invocable(_interval_b,[](){ return doption(_name="solvereigen.interval-b"); } );

    std::string const& solver = args.get_else_invocable(_solver,[](){ return soption(_name="solvereigen.solver"); } );
    std::string const& problem = args.get_else_invocable(_problem,[](){ return soption(_name="solvereigen.problem"); } );
    std::string const& transform = args.get_else_invocable(_transform,[](){ return soption(_name="solvereigen.transform"); } );
    std::string const& spectrum = args.get_else_invocable(_spectrum,[](){ return soption(_name="solvereigen.spectrum"); } );
    int maxit = args.get_else_invocable(_maxit,[](){ return ioption(_name="solvereigen.maxiter"); } );
    double tolerance = args.get_else_invocable(_tolerance,[](){ return doption(_name="solvereigen.tolerance"); } );
    bool verbose = args.get_else_invocable(_verbose,[](){ return boption(_name="solvereigen.verbose"); } );

    typedef std::shared_ptr<Vector<double> > vector_ptrtype;
    //std::shared_ptr<SolverEigen<double> > eigen = SolverEigen<double>::build(  backend );
    using size_type = typename SolverEigen<double>::size_type;
    std::shared_ptr<SolverEigen<double> > eigen = SolverEigen<double>::build();
    eigen->setEigenSolverType( (EigenSolverType)EigenMap[solver] );
    eigen->setEigenProblemType( (EigenProblemType)EigenMap[problem] );
    eigen->setPositionOfSpectrum( (PositionOfSpectrum)EigenMap[spectrum] );
    eigen->setNumberOfEigenValues( nev );
    eigen->setNumberOfEigenValuesConverged( ncv );
    eigen->setMaximumProjectedDimension( (mpd>0)?mpd:invalid_v<size_type> );
    eigen->setInterval( interval_a, interval_b );
    eigen->setMaxIterations( maxit );
    eigen->setSpectralTransform( (SpectralTransformType)EigenMap[transform] );
    eigen->setTolerance( tolerance );
    eigen->setMapRow( formA.matrixPtr()->mapRowPtr() );
    eigen->setMapCol( formA.matrixPtr()->mapColPtr() );

    LOG(INFO) << "number of eigen values = " << nev << "\n";
    LOG(INFO) << "number of eigen values converged = " << ncv << "\n";
    LOG(INFO) << "number of eigen value solver iterations = " << maxit << "\n";
    LOG(INFO) << "eigenvalue tolerance = " << tolerance << "\n";

    unsigned int nconv, nits;
    std::vector<double> err( ncv );
    boost::tie( nconv, nits, err )  = eigen->solve( _matrixA=formA.matrixPtr(),
                                                    _matrixB=formB.matrixPtr(),
                                                    _maxit=maxit,
                                                    _tolerance=tolerance );

    if ( verbose )
    {
        int i = 0;
        for( auto e : err )
        {
            if ( Environment::worldComm().isMasterRank() )
            {
                std::cout << i++ << "-th mode error ||A x - lambda B x ||/||x|| = " << e << "\n";
            }
        }
    }

    auto modes = eigen->eigenModes();
    auto Xh = formA.trialSpace();
    auto femodes = std::vector<std::pair<double, decltype( Xh->element() )> >( modes.size(), std::make_pair(0.,Xh->element() ) );

    if ( !modes.empty() )
    {
        LOG(INFO) << "eigenvalue " << 0 << " = (" << modes.begin()->second.get<0>() << "," <<  modes.begin()->second.get<1>() << ")\n";

        int i = 0;
        for( auto const& mode : modes )
        {
            if ( Environment::worldComm().isMasterRank() )
                std::cout << " -- eigenvalue " << i << " = (" << mode.second.get<0>() << "," <<  mode.second.get<1>() << ")\n";
            femodes[i].first = mode.second.get<0>();
            femodes[i].second = *mode.second.get<2>();
            femodes[i].second.close();
            ++i;
        }
    }
    return femodes;

}


} // Feel
#endif /* __SolverEigen_H */
