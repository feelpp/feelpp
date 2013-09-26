/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-04-07

  Copyright (C) 2005,2006 EPFL

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
   \file iteration.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-04-07
 */
#ifndef __Iteration_H
#define __Iteration_H 1

#include <boost/shared_ptr.hpp>

namespace Feel
{
/*!
  \class Iteration
  \brief brief description

  The Iteration object calculates whether the solution has reached
  the desired accuracy, or whether the maximum number of iterations
  has been reached. The method \c isFinished() checks both convergence and
  number of iterations. The method \c isConverged() only checks convergence.
  The \c isFirst() method is used to determine the first iteration of the loop.


  The following notation will be used
  @li \f$r\f$ the residual
  @li \f$\epsilon\f$ the relative precision
  @li \f$I\f$ the number of already performed iterations
  @li \f$M\f$ the maximum number of iterations allowed
  @author Christophe Prud'homme
  @see
  @version $Id: Iteration.hpp,v 1.6 2002/08/22 13:09:56 prudhomm Exp $
*/
template<typename Real>
class Iteration
{
public:


    /** @name Typedefs
     */
    //@{

    /**
       \brief Numerical Type
    */
    typedef Real NumericalType;
    typedef  typename ublas::type_traits<Real>::value_type value_type;
    typedef typename ublas::type_traits<Real>::real_type real_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     *  \brief create a new instance
     */
    static Iteration<NumericalType>* New()
    {
        return new Iteration<NumericalType>;
    }

    Iteration( Iteration const& iter )
        :
        __iterations( iter.__iterations ),
        __max_iter( iter.__max_iter ),
        __residual( iter.__residual ),
        __precision( iter.__precision ),
        __norm_init( iter.__norm_init )
    {
        // do nothing here
    }


    //! destructor
    virtual ~Iteration()
    {
        // do nothing here
    }

    //@}

    /** @name Operators
     */
    //@{

    //! copy operator
    Iteration& operator=( Iteration<NumericalType> const& iter )
    {
        if ( this == &iter )
        {
            return *this;
        }

        __precision = iter.__precision;
        __max_iter = iter.__max_iter;
        __residual = iter.__residual;
        __norm_init = iter.__norm_init;
        return *this;
    }

    //! prefix ++ operator
    void operator++() throw()
    {
        ++__iterations;
    }

    //@}

    /** @name Accessors
     */
    //@{

    int numberOfIterations() const
    {
        return __iterations;
    }

    /**
       \brief get the Residual
    */
    real_type residual() const throw()
    {
        return __residual;
    }

    real_type relativePrecision() const
    {
        return __precision;
    }

    int maximumNumberOfIterations() const
    {
        return __max_iter;
    }

    NumericalType initialResidual() const
    {
        return __norm_init;
    }

    real_type relaxation () const
    {
        return M_relaxation;
    }
    int iteration() const
    {
        return __iterations;
    }

    /** @name  Mutators
     */
    //@{

    //! set the Max number of iterations
    /**
       @param m max number of iterations to perform
    */
    void setMaximumNumberOfIterations( int m ) throw()
    {
        __max_iter = m;
    }

    //! set the relative precision to reach
    /**
       @param p precision
    */
    void setRelativePrecision( NumericalType p ) throw()
    {
        __precision = p;
    }

    //! initial norm for the residual
    /*!
      \param ninit initial norm for  the residual
    */
    void setInitialResidual( NumericalType ninit ) throw()
    {
        __norm_init = ninit;
    }

    void setRelaxation ( real_type __w )
    {
        M_relaxation = __w;
    }
    //@}

    /** @name  Methods
     */
    //@{

    /**
     * \brief tells if  the iteration finished
     *
     * Three cases can occur:
     * @li if \f$ r < \epsilon \f$ then the iteration is over
     * @li if \f$ I > M \f$ then the iteration is over
     * @li else the iteration must continue
     *
     * @param r residual to test the convergence
     * @param verbose true for verbose output, false otherwise
     * @return false if not finished and true otherwise
     */
    bool isFinished( NumericalType r, bool verbose = false ) //throw(SExceptionSolverHasNotConverged)
    {
        bool ret = false;

        if ( this->isConverged( r ) )
        {
            ret = true;
        }

        else if ( __iterations >= __max_iter )
        {
            ret = true;
        }

        else
        {
#if 0
            handleEvents( true );
            std::string why = "Solver has not converged";
            SExceptionSolverHasNotConverged __e( why.c_str(), LOCATION );
            __e.setNumberOfIterations( __iterations );
            __e.setResidual( __residual );
            throw __e;
#endif
        }

        handleEvents( ret, verbose );

        return ret;
    }


    template<typename VectorX>
    bool isFinished( const VectorX& r, bool verbose = false ) //throw(SExceptionSolverHasNotConverged)
    {
        bool ret = false;

        if ( this->isConverged( r ) )
        {
            ret = true;
        }

        else if ( __iterations >= __max_iter )
        {
            ret = true;
        }

        else
        {
#if 0
            handleEvents( true );
            std::string why = "Solver has not converged";
            SExceptionSolverHasNotConverged __e( why.c_str(), LOCATION );
            __e.setNumberOfIterations( __iterations );
            __e.setResidual( __residual );
            throw __e;
#endif
        }

        handleEvents( ret, verbose );

        return ret;
    }

    bool isConverged( NumericalType r ) throw()
    {
        __residual = r / __norm_init;
        return ( __residual <= __precision );
    }

    template<typename VectorX> bool isConverged( VectorX const& x ) throw()
    {
        __residual = ublas::norm_2( x ) / __norm_init;
        return ( __residual <= __precision );
    }

    bool isFirst() const
    {
        return ( __iterations == 0 );
    }

    void reset()
    {
        __iterations = 0;
    }

    //@}

protected:

    /**
       Default constructor.

    */
    Iteration()
        :
        __iterations( 0 ),
        __max_iter( 0 ),
        __residual( 0 ),
        __precision( 0 ),
        __norm_init( 1.0 ),
        M_relaxation( 1.0 )
    {
        // do nothing here
    }

    virtual void handleEvents( bool __is_finished, bool verbose )
    {
        if ( __iterations == 0 )
        {
            //SEvent  startEv( IterationStarted, new IterationStarted<value_type>( __residual, __iterations) );
            //SSubject::notifyObservers( &startEv );
            if ( verbose )
                std::cout << "iteration " << __iterations << " : " << residual() << "\n";
        }

        if ( __is_finished == true )
        {
            //SEvent  finishEv( IterationFinished, new IterationFinished<value_type>( __residual, __iterations) );
            //SSubject::notifyObservers( &finishEv );
            if ( verbose )
                std::cout << "iteration " << __iterations << " : " << residual() << "\n";
        }

        else
        {
            if ( verbose )
                std::cout << "iteration " << __iterations << " : " << residual() << "\n";

            //SEvent  aEvent( IterationUpdated, new IterationUpdated<value_type>( __residual, __iterations) );
            //SSubject::notifyObservers( &aEvent );
        }
    }
private:

    int __iterations;
    int __max_iter;

    real_type __residual;
    real_type __precision;
    real_type __norm_init;

    real_type M_relaxation;
};

typedef Iteration<double> iteration_type;
typedef boost::shared_ptr<iteration_type> iteration_ptrtype;
}


#endif /* __Iteration_H */
