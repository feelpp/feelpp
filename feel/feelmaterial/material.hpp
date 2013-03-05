/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-05-25

  Copyright (C) 2008 Universite Joseph Fourier (Grenoble I)

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
   \file material.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-25
 */
#ifndef __Material_H
#define __Material_H 1

#include <string>
#include <boost/shared_ptr.hpp>

namespace Feel
{
/**
 * \class Material
 * \brief Material base class
 *
 * @author Christophe Prud'homme
 * @see
 */
class Material
{
public:


    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    Material( std::string const& name );
    Material( Material const & );
    virtual ~Material();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    //! \return the name of the material
    std::string const& name() const
    {
        return M_name;
    }

    //! thermal conductivity in \f$ W/(m*K) \f$
    virtual double k() const = 0;

    //! density in \f$ kg/m^3 \f$
    virtual double rho() const = 0;

    //! thermal capacity in \f$ J/(kg*K) \f$
    virtual double C() const = 0;

    //! Poisson coefficient
    virtual double nu() const = 0;

    //! Young modulus in \f$ Pa \f$
    virtual double E() const = 0;

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

    //! name of the material
    std::string M_name;

private:

};

typedef Material material_type;
typedef boost::shared_ptr<Material> material_ptrtype;

} // Feel


#endif /* __Material_H */
