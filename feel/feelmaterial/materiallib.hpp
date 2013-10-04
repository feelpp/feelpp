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
   \file materiallib.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-25
 */
#ifndef __MaterialLib_H
#define __MaterialLib_H 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>
#include <feel/feelmaterial/material.hpp>



namespace Feel
{
namespace detail
{
template<typename T>
Material* createMaterial()
{
    return new T;
}
}
/**
 * \class MaterialLib
 * \brief Material library
 *
 *  @author Christophe Prud'homme
 *  @see
 */
class MaterialLib
{
public:


    /** @name Typedefs
     */
    //@{

    typedef Singleton< Feel::Factory< Material, std::string > > factory_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    MaterialLib();
    MaterialLib( po::variables_map const& vm );
    MaterialLib( MaterialLib const & );
    ~MaterialLib();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    static material_ptrtype material( std::string const& name );


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

typedef Singleton< Feel::Factory< Material, std::string > > MaterialFactory;

} // Feel
#endif /* __MaterialLib_H */
