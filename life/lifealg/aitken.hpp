/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Goncalo Pena <goncalo.pena@epfl.ch>
       Date: 15-07-2008

  Copyright (C) 2007-2008 Universite Joseph Fourier (Grenoble I)

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


#ifndef __AitkenExtrapolation
#define __AitkenExtrapolation 1

#include <life/lifediscr/functionspace.hpp>
#include <life/lifealg/vector.hpp>


namespace Life
{
template< typename fs_type >

class Aitken
{

public:

    typedef fs_type functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;

    Aitken( functionspace_ptrtype& _Xh, double _failsafeParameter = 1  )
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

    ~Aitken() {}

    void initialize( element_type& residual, element_type& elem )
    {
        previousResidual = residual;
        previousElement = elem;
    }

    void setElement( element_type& residual, element_type& elem )
    {
        currentResidual = residual;
        currentElement = elem;
    }

    double calculateParameter()
    {
        element_type aux( Xh, "aux");

        aux = currentResidual;
        aux -= previousResidual;

        double scalar = inner_product( aux, aux );

        aux.scale( 1.0/scalar );

        element_type aux2( Xh, "aux2");

        aux2 = currentElement;
        aux2 -= previousElement;

        scalar = inner_product( aux2, aux );

        if ( scalar > failsafeParameter )
            scalar = previousParameter;

        if ( scalar < 1-failsafeParameter )
            scalar = previousParameter;

        previousParameter = scalar;

        return scalar;
    }

    void relaxationStep( element_type& new_elem )
    {
        new_elem = currentResidual;
        new_elem.scale( -previousParameter );

        new_elem += currentElement;
    }

    void shiftRight()
    {
        previousResidual = currentResidual;
        previousElement = currentElement;
    }

    void resetPreviousParameter()
    {
        previousParameter = failsafeParameter;
    }

private:

    functionspace_ptrtype Xh;

    double failsafeParameter, previousParameter;

    element_type previousResidual, previousElement, currentResidual, currentElement;

};


} // End namespace Life

#endif
