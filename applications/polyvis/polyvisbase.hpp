/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-04-21

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file polyvisbase.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-04-21
 */
#ifndef __PolyvisBase_H
#define __PolyvisBase_H 1

#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>

namespace Life
{
/**
 * \class PolyvisBase
 * \brief Polynomial visualisation base class
 *
 * @author Christophe Prud'homme
 * @see
 */
class PolyvisBase
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{

    //! polyvis factory
    typedef Life::Singleton< Life::Factory< PolyvisBase, std::string > > factory_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    PolyvisBase();
    //! copy constructor
    PolyvisBase( PolyvisBase const & );
    //! destructor
    virtual ~PolyvisBase() {}

    static PolyvisBase* New( po::variables_map const& vm )
    {
        //return factory_type::instance().createObject( vm[] );
        //return p;
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    PolyvisBase& operator=( PolyvisBase const & o)
    {
        if (this != &o )
            {
            }
        return *this;
    }
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

    virtual void run() = 0;

    //@}



protected:

private:

};
} // Life
#endif /* __PolyvisBase_H */
