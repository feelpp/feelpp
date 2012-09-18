/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-01-04

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file systemexplicit.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#ifndef __SystemExplicit_H
#define __SystemExplicit_H 1

#include <feel/feelsystem/system.hpp>

namespace Feel
{
/**
 * \class SystemExplicit
 * \brief describes an explicit system
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename SpaceType>
class SystemExplicit : public System<SpaceType>
{
    typedef System<SpaceType> super;
public:


    /** @name Constants
     */
    //@{

    static const uint16_type Dim = super::Dim;

    //@}

    /** @name Typedefs
     */
    //@{

    typedef SystemExplicit<SpaceType> system_type;

    typedef typename super::value_type value_type;
    typedef typename super::functionspace_type functionspace_type;
    typedef typename super::functionspace_type functionspace_ptrtype;
    typedef typename super::element_type element_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    SystemExplicit( functionspace_ptrtype Xh, po::variables_map const& vm ) : super( Xh, vm ) {}
    SystemExplicit( SystemExplicit const & se ) : super( se ) {}
    ~SystemExplicit() {}

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


    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

private:

};
} // Feel
#endif /* __SystemExplicit_H */
