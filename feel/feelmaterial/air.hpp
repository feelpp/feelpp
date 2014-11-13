/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-05-25

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file air.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-25
 */
#include <cmath>

#include <boost/any.hpp>
#include <map>
#include <utility>

#include <feel/feelmaterial/material.hpp>

namespace Feel
{

/**
 * \class Air
 *  \brief Air material
 *
 *  @author Christophe Prud'homme
 * @see
 */
class Air : public Material
{
public:


    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    Air() : Material( "Air" ) {}
    Air( Air const & m ): Material( m ) {}
    ~Air() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //! thermal conductivity in \f$ W/(m*K) \f$
    virtual double k() const
    {
        double T = 273; // K (default)
        return pow( 10,( 0.8616*log10( std::abs( T ) )-3.7142 ) );
    }

    //! density in \f$ kg/m^3 \f$
    virtual double rho() const
    {

        return 1.2;
    }

    //! thermal capacity in \f$ J/(kg*K) \f$
    virtual double C() const
    {
        double T=273;
        return 0.0769*T+1076.9;
    }

    //! Poisson coefficient
    virtual double nu() const
    {
        return 1.7*1e-5;
    }

    //! Young modulus in \f$ Pa \f$
    virtual double E() const
    {
        return -1;
    }


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

private:

};

} // Feel

