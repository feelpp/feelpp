/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-10-13

  Copyright (C) 2007,2009 Université de Grenoble 1
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
   \file typeInfo.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-10-13
 */
#include <cassert>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/typeinfo.hpp>

namespace Feel
{
TypeInfo::TypeInfo()
{
    class Nil {};
    M_info = &typeid( Nil );
    assert( M_info != 0 );

}

TypeInfo::TypeInfo( const std::type_info& ti )
    :
    M_info( &ti )
{
    assert( M_info != 0 );
}

TypeInfo::TypeInfo( TypeInfo const& ti )
    :
    M_info( ti.M_info )
{
    assert( M_info != 0 );
}

TypeInfo::~TypeInfo()
{
    ; // do nothing here
}

bool
TypeInfo::before( const TypeInfo& rhs ) const
{
    assert( M_info != 0 );
    return M_info->before( *rhs.M_info );
}

const std::type_info&
TypeInfo::typeInfo() const
{
    assert( M_info != 0 );
    return *M_info;
}

const char*
TypeInfo::name() const
{
    assert( M_info != 0 );
    return M_info->name();
}

// Comparison operators


}
