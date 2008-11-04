/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-05-25

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file castiron.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-05-25
 */
#include <boost/any.hpp>
#include <map>
#include <utility>

#include <boost/plugin/export_plugin.hpp>

#include <life/lifematerial/material.hpp>


namespace Life
{

/**
 * \class CastIron
 *  \brief Cast iron material
 *
 *  @author Christophe Prud'homme
 * @see
 */
class CastIron : public Material
{
public:


    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    CastIron() : Material( "Cast Iron" ) {}
    CastIron( CastIron const & m ): Material( m ) {}
    ~CastIron() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //! thermal conductivity in \f$ W/(m*K) \f$
    virtual double k() const { return -1;}

    //! density in \f$ kg/m^3 \f$
    virtual double rho() const { return 7000; }

    //! thermal capacity in \f$ J/(kg*K) \f$
    virtual double C() const { return -1; }

    //! Poisson coefficient
    virtual double nu() const { return 0.25; }

    //! Young modulus in \f$ Pa \f$
    virtual double E() const { return 140e9; }


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
BOOST_PLUGIN_EXPORT(Material, CastIron, "Cast Iron");


} // Life
