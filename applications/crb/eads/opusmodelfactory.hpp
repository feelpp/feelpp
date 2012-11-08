/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-07-20

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
   \file opusmodelfactory.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-07-20
 */
#ifndef __OpusModelFactory_H
#define __OpusModelFactory_H 1

#include <opusmodelbase.hpp>

namespace Feel
{
/**
 * \addtogroup Models
 * \\@{
 */
/**
 * \class OpusModelFactory
 * \brief Opus model factory class
 *
 * This class is in charge of creating new OpusModel classes which
 * have been registered through the Factory pattern.
 *
 * @author Christophe Prud'homme
 */
class OpusModelFactory
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{

    typedef OpusModelBase opusmodel_type;
    typedef boost::shared_ptr<opusmodel_type> opusmodel_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    OpusModelFactory();
    //! copy constructor
    OpusModelFactory( OpusModelFactory const & );
    //! destructor
    ~OpusModelFactory();

    //! instantiate new opus model from vm
    static opusmodel_ptrtype New( Feel::po::variables_map const& vm );

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    OpusModelFactory& operator=( OpusModelFactory const & o )
    {
        if ( this != &o )
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


    //@}



protected:

private:

};
/** \\@} */
} // Feel
#endif /* __OpusModelFactory_H */
