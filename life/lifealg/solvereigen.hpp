/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-04

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-04
 */
#ifndef __SolverEigen_H
#define __SolverEigen_H 1

#include <life/lifealg/enums.hpp>
#include <life/lifecore/traits.hpp>
#include <life/lifealg/vector.hpp>
#include <life/lifealg/matrixsparse.hpp>

namespace Life
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
    typedef boost::shared_ptr<solvereigen_type> solvereigen_ptrtype;
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
    static boost::shared_ptr<SolverEigen<value_type> > build(const SolverPackage solver_package = SOLVERS_SLEPC);

    /**
     * Builds a \p SolverEigen using the linear solver package
     * specified by \p vm
     * \param vm variables_map that contains the command line options and their defaults
     * \param prefix string that allows for various options of the same type 
     */
    static boost::shared_ptr<SolverEigen<value_type> > build( po::variables_map const& vm, 
                                                              std::string const& prefix = std::string() );

    /**
     * Returns the \p ith eigenvalue (real and imaginary part), and
     * copies the \ ith eigen vector to the solution vector.
     */
    virtual std::pair<real_type, real_type> eigenPair (unsigned int i, Vector<value_type> &solution) = 0;

    /**
     * Returns the \p ith eigenvalue (real and imaginary part), and
     * copies the \ ith eigen vector to the solution vector.
     */
    virtual std::pair<real_type, real_type> eigenPair (unsigned int i,
                                                       boost::shared_ptr<Vector<value_type> > &solution)
        {
            return this->eigenPair( i, *solution );
        }

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * Returns the type of eigensolver to use.
     */
    EigenSolverType eigenSolverType () const { return M_eigen_solver_type; }

    /**
     * Returns the type of the eigen problem.
     */
    EigenProblemType eigenProblemType () const { return M_eigen_problem_type;}

    /**
     * Returns the position of the spectrum to compute.
     */
    PositionOfSpectrum postitionOfSpectrum () const { return M_position_of_spectrum;}

    //@}

    /** @name  Mutators
     */
    //@{


    /**
     * Sets the type of eigensolver to use.
     */
    void setEigenSolverType (const EigenSolverType est) { M_eigen_solver_type = est; }

    /**
     * Sets the type of the eigenproblem.
     */
    void setEigenProblemType ( EigenProblemType ept) { M_eigen_problem_type = ept;}

    /**
     * Sets the position of the spectrum.
     */
    void setPositionOfSpectrum (PositionOfSpectrum pos) { M_position_of_spectrum= pos; }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * @returns true if the data structures are
     * initialized, false otherwise.
     */
    bool initialized () const { return M_is_initialized; }

    /**
     * Release all memory and clear data structures.
     */
    virtual void clear () {}

    /**
     * Initialize data structures if not done so already.
     */
    virtual void init () = 0;



    /**
     * Solves the standard eigen problem and returns the
     * number of converged eigenpairs and the number
     * of iterations.
     */
    virtual std::pair<unsigned int, unsigned int> solve (MatrixSparse<value_type> &matrix_A,
                                                         int nev,
                                                         int ncv,
                                                         const double tol,
                                                         const unsigned int m_its) = 0;

    /**
     * Solves the standard eigen problem and returns the
     * number of converged eigenpairs and the number
     * of iterations.
     */
    std::pair<unsigned int, unsigned int> solve ( boost::shared_ptr<MatrixSparse<value_type> > &matrix_A,
                                                  int nev,
                                                  int ncv,
                                                  const double tol,
                                                  const unsigned int m_its) 
        {
            return this->solve( *matrix_A, nev, ncv, tol, m_its );
        }

    /**
     * Solves the generalized eigen value problem \f$A x = \lambda
     * Bx\f$ and returns the number of converged eigenpairs and the
     * number of iterations.
     */
    virtual std::pair<unsigned int, unsigned int> solve (MatrixSparse<value_type> &matrix_A,
                                                         MatrixSparse<value_type> &matrix_B,
                                                         int nev,
                                                         int ncv,
                                                         const double tol,
                                                         const unsigned int m_its) = 0;

    /**
     * Solves the generalized eigen value problem \f$A x = \lambda
     * Bx\f$ and returns the number of converged eigenpairs and the
     * number of iterations.
     */
    std::pair<unsigned int, unsigned int> solve ( boost::shared_ptr<MatrixSparse<value_type> > &matrix_A,
                                                  boost::shared_ptr<MatrixSparse<value_type> > &matrix_B,
                                                  int nev,
                                                  int ncv,
                                                  const double tol,
                                                  const unsigned int m_its)
        {
            return this->solve( *matrix_A, *matrix_B, nev, ncv, tol, m_its );
        }
    



    //@}





protected:


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
     * Flag indicating if the data structures have been initialized.
     */
    bool M_is_initialized;

};
/**
 * defines solver eigen options
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
 * \param prefix prefix allows to prefix options 
 */
po::options_description solvereigen_options( std::string const& prefix = "" );
} // Life
#endif /* __SolverEigen_H */
