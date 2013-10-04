/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-10-13

  Copyright (C) 2009 Université de Grenoble 1
  Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano

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
   \file typeInfo.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-10-13
 */
#ifndef __TypeInfo_H
#define __TypeInfo_H 1

#include <typeinfo>

namespace Feel
{
/*!
  \class TypeInfo
 *\ingroup Core
 *\brief wrapper for std::type_info

  \sa SFactory, type_info

  @author Christophe Prud'homme
  @version $Id: typeInfo.hpp,v 1.1 2004/10/13 10:18:52 prudhomm Exp $
*/
class TypeInfo
{
public:


    /** @name Typedefs
     */
    //@{

    //@}

    /** @name Constructors, destructor
     */
    //@{

    TypeInfo();
    TypeInfo( const std::type_info& ); // non-explicit
    TypeInfo( TypeInfo const & );
    ~TypeInfo();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    //! Access for the wrapped \c std::type_info
    const std::type_info& typeInfo() const;

    //! get the \c name()
    const char* name() const;

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //! Compatibility functions
    bool before( const TypeInfo& rhs ) const;


    //@}

private:

    const std::type_info* M_info;
};

inline bool operator==( const TypeInfo& lhs, const TypeInfo& rhs )
{
    return lhs.typeInfo() == rhs.typeInfo();
}

inline bool operator<( const TypeInfo& lhs, const TypeInfo& rhs )
{
    return lhs.before( rhs );
}

inline bool operator!=( const TypeInfo& lhs, const TypeInfo& rhs )
{
    return !( lhs == rhs );
}

inline bool operator>( const TypeInfo& lhs, const TypeInfo& rhs )
{
    return rhs < lhs;
}

inline bool operator<=( const TypeInfo& lhs, const TypeInfo& rhs )
{
    return !( lhs > rhs );
}

inline bool operator>=( const TypeInfo& lhs, const TypeInfo& rhs )
{
    return !( lhs < rhs );
}

}// end namespace Feel

#endif /* __TypeInfo_H */
