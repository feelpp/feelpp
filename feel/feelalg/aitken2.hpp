/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes
       Date: 24-03-2011

  Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)

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


#ifndef __AitkenExtrapolation2
#define __AitkenExtrapolation2 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelalg/vector.hpp>


namespace Feel
{
/**
 * \class Aitken
 * \brief Aitken relaxation method for fixed point iterations
 *
 * \code
 * space_ptrtype Xh;
 * space_type::element_type residual( Xh );
 * space_type::element_type u_old( Xh );
 * space_type::element_type u_new( Xh );
 * Aitken<space_type> aitken( Xh );
 * // initialize aitken
 * aitken.initialize( residual, u_new );
 * // reset aitken parameter before entering the fixed point loop
 * aitken.resetPreviousParameter();
 * // fixed point loop
 * for( int i = 0; i < niter; ++i )
 * {
 *   // do some computation
 *   aitken.SetElement( residual, u_new );
 *   theta = aitken.calculateParameter();
 *   // exploit aitken relaxation parameter
 *   u = theta * u_new + (1-theta)* u_old;
 *   aitken.shiftRight();
 *   u_old = u_new;
 * }
 * \endcode
 *
 * \author Vincent Chabannes
 */
template< typename fs_type >
class Aitken
{

public:

    typedef fs_type functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;

    typedef typename functionspace_type::template Element<typename functionspace_type::value_type,
            typename VectorUblas<typename functionspace_type::value_type>::range::type > element_range_type;


    /**
     * Constructor
     *
     * the \p _failsafeParameter is set to 1 by default. The _failsafeParameter
     * parameter provides an upper value for the relaxation parameter that will
     * not be exceeded.
     *
     * \param _Xh the function space from which the element will be used
     * \param _failsafeParameter fail safe parameter value
     */
    Aitken( functionspace_ptrtype _Xh, double _failsafeParameter = 0.1  )
        :
        Xh( _Xh ),
        failsafeParameter( _failsafeParameter ),
        previousParameter( _failsafeParameter ),
        previousResidual( Xh, "previous residual" ),
        previousElement( Xh, "previous element" ),
        currentResidual( Xh, "current residual" ),
        currentElement( Xh, "current element" )
    {
    }

    /**
     * copy constructor
     */
    Aitken( Aitken const& tc )
        :
        Xh( tc.Xh ),
        failsafeParameter( tc.failsafeParameter ),
        previousParameter( tc.previousParameter ),
        previousResidual( tc.previousResidual ),
        previousElement( tc.previousElement ),
        currentResidual( tc.currentResidual ),
        currentElement( tc.currentElement )
    {
    }

    /**
     * destructor
     */
    ~Aitken() {}

    /**
     * initiliaze the aitken algorithm
     * \param residual  previous residual
     * \param elem previous element
     */
    void initialize( element_type const& residual, element_type const& elem )
    {
        previousResidual = residual;
        previousElement = elem;
    }

    void initialize( element_type const& residual, element_range_type const& elem )
    {
        previousResidual = residual;
        previousElement.zero();
        previousElement.add( 1.,elem );
        /*previousElement = vf::project(previousElement.functionSpace(),
                                      elements(previousElement.mesh()),
                                      vf::idv(elem) );*/
    }

    /**
     * Set the current element
     * \param residual current residual
     * \param elem current element
     */
    void setElement( element_type const& residual, element_type const& elem )
    {
        currentResidual = residual;
        currentElement = elem;
    }

    void setElement( element_type const& residual, element_range_type const& elem )
    {
        currentResidual = residual;
        currentElement.zero();
        currentElement.add( 1.,elem );
        /*currentElement = vf::project(currentElement.functionSpace(),
                                     elements(currentElement.mesh()),
                                     vf::idv(elem) );*/
    }


    /**
     * \return the Aitken parameter
     */
    double calculateParameter();

    /**
     * Do a relaxation step
     * \param new_elem new element to compute the relaxation step
     */
    void relaxationStep( element_type& new_elem )
    {
        new_elem = currentResidual;
        new_elem.scale( -previousParameter );

        new_elem += currentElement;
    }

    /**
     * shift current step to previous step. After the call, we are ready for the
     * next step.
     */
    void shiftRight()
    {
        previousResidual = currentResidual;
        previousElement = currentElement;
    }

    /**
     * reset the previous parameter
     */
    void resetPreviousParameter()
    {
        previousParameter = failsafeParameter;
    }

private:

    /**
     * function space
     */
    functionspace_ptrtype Xh;

    double failsafeParameter, previousParameter;

    element_type previousResidual, previousElement, currentResidual, currentElement;

};


template< typename fs_type >
double
Aitken<fs_type>::calculateParameter()
{
    element_type aux( Xh, "aux" );

    aux = currentResidual;
    aux -= previousResidual;

    double scalar = inner_product( aux, aux );

    aux.scale( 1.0/scalar );

    element_type aux2( Xh, "aux2" );

    aux2 = currentElement;
    aux2 -= previousElement;

    scalar = inner_product( previousResidual , aux );
    scalar = -previousParameter*scalar;
#if 1

    if ( scalar > 1 )
        scalar = 1;//-failsafeParameter;

    if ( scalar < 0 )
        scalar = failsafeParameter;

#endif
    previousParameter = scalar;

    return scalar;
}


} // End namespace Feel

#endif
