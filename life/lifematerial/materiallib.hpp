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
   \file materiallib.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-05-25
 */
#ifndef __MaterialLib_H
#define __MaterialLib_H 1

#include <life/lifecore/life.hpp>
#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>
#include <life/lifematerial/material.hpp>



namespace Life
{
namespace detail
{
template<typename T>
Material* createMaterial() { return new T; }
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

    typedef Singleton< Life::Factory< Material, std::string > > factory_type;

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

typedef Singleton< Life::Factory< Material, std::string > > MaterialFactory;
po::options_description material_options( std::string const& prefix = "" );

} // Life
#endif /* __MaterialLib_H */
