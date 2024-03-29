/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-12-31

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file functionspacebase.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-12-31
 */
#ifndef __FunctionSpaceBase_H
#define __FunctionSpaceBase_H 1

#include <feel/feelcore/journalwatcher.hpp>
#include <feel/feelcore/commobject.hpp>
#include <feel/feelalg/datamap.hpp>

namespace Feel
{

/**
 * @class FunctionSpaceBase
 * @brief base class for FunctionSpace
 *
 * @author Christophe Prud'homme
 * @see
 */
class FEELPP_EXPORT FunctionSpaceBase
    : public CommObject,
      public JournalWatcher
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{
    using super = CommObject;
    using super2 = JournalWatcher;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    FunctionSpaceBase()
        :
        super( Environment::worldCommPtr() ), super2( "FunctionSpace" )
    {}
    explicit FunctionSpaceBase( std::string const& name, worldcomm_ptr_t const& w )
        :
        super( w ), super2( "FunctionSpace", name )
    {}
    FunctionSpaceBase( FunctionSpaceBase const& ) = default;
    FunctionSpaceBase( FunctionSpaceBase && ) = default;
    //! destructor
    ~FunctionSpaceBase() override = default;

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    FunctionSpaceBase& operator=( FunctionSpaceBase const & o )
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

    virtual datamap_ptr_t<> mapPtr() const = 0;

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    class ElementBase {};

private:
    //! Private Methods
    //! @{

    //! @}
};

}
#endif /* __FunctionSpaceBase_H */
