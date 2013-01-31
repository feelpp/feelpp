/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-10-04

  Copyright (C) 2004 EPFL
  Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)

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
   \file solverumfpack.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-10-04
 */
#ifndef __SolverUMFPACK_H
#define __SolverUMFPACK_H 1

#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/matrixtriplet.hpp>

#if defined(FEELPP_HAS_UMFPACK)

extern "C"
{
#if defined (FEELPP_HAS_SUITESPARSE_UMFPACK_H)
#include <suitesparse/umfpack.h>
#elif defined (FEELPP_HAS_UFSPARSE_UMFPACK_H)
#include <ufsparse/umfpack.h>
#elif defined (FEELPP_HAS_UMFPACK_UMFPACK_H)
#include <umfpack/umfpack.h>
#else
#include <umfpack.h>
#endif
};
#endif


namespace Feel
{
/*!
  \class SolverUMFPACK
  \brief Interface for the UMFPACK Solver

  UMFPACK is a direct Solver for (un)symmetric problem \f$ A x = b \f$.

  @author Christophe Prud'homme
*/
class SolverUMFPACK
{
public:


    /** @name Typedefs
     */
    //@{

    typedef double value_type;

    //typedef boost::numeric::ublas::compressed_matrix<value_type, boost::numeric::ublas::column_major> matrix_type;
    typedef MatrixTriplet<double> matrix_type;

    typedef ublas::vector<value_type> array_type;



    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
       default constructor

       it sets the umfpack print level to the maximum (ie 6)
     */
    SolverUMFPACK();

    SolverUMFPACK( SolverUMFPACK const & umfpackSolver );

    ~SolverUMFPACK();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{

    void setMatrix( matrix_type const& m );

    /**
     * set the umfpack strategy
     * possible values are :
     * - UMFPACK_STRATEGY_AUTO
     * - UMFPACK_STRATEGY_UNSYMMETRIC
     * - UMFPACK_STRATEGY_SYMMETRIC
     * - UMFPACK_STRATEGY_2BY2
     *
     * @see UMFPACK user guide for further details
     */
    void setStrategy( int strategy );

    //@}

    /** @name  Methods
     */
    //@{


    //! solve A X = B
    /*!

    \param __X  the solution
    \param __B    the right hand side
    \return the number of iterations
    */
    void solve( array_type& __X, array_type const& __B );

    //! report some info about umfpack
    void reportInfo();


    /**
     *  report status of umfpack
     *
     *  \param status status integer returned by umfpack routines
     */
    void reportStatus( int status );

    //@}

private:

    void prepareSolve();


private:

    class Pimpl;
    boost::shared_ptr<Pimpl> M_p;

    bool M_matrix_reset;
    bool M_matrix_values_reset;

    void *M_symbolic;
    void *M_numeric;

    double* M_Control;
    double* M_Info;

};
}
#endif /* __SolverUMFPACK_H */

