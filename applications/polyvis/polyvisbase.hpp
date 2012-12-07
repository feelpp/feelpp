/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-21

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
   \file polyvisbase.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-21
 */
#ifndef __PolyvisBase_H
#define __PolyvisBase_H 1

#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

namespace Feel
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
    typedef Feel::Singleton< Feel::Factory< PolyvisBase, std::string > > factory_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    PolyvisBase();

    //! copy constructor
    PolyvisBase( PolyvisBase const & );

    //! destructor
    virtual ~PolyvisBase();

    /**
     * initialize the polynomial visualisation class
     * \param vm variables map from the command line
     */
    virtual void init( po::variables_map const& vm );

    /**
     * create a new polynomial visualisation class
     * \param vm variables map from the command line
     */
    static PolyvisBase* New( po::variables_map const& vm );

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    PolyvisBase& operator=( PolyvisBase const & o );

    //@}

    /** @name Accessors
     */
    //@{

    //! return the variables map
    po::variables_map vm() const
    {
        return M_vm;
    }

    //! return the polynomial name
    std::string name() const
    {
        return M_pname;
    }

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

    //! variables map
    po::variables_map M_vm;

    //! polynomial family name
    std::string M_pname;
private:

};
} // Feel
#endif /* __PolyvisBase_H */
