/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2004-10-13

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2004-10-13
 */
#include <cassert>

#include <life/lifecore/life.hpp>
#include <life/lifecore/typeinfo.hpp>

namespace Life
{
TypeInfo::TypeInfo()
{
    class Nil {};
    _M_info = &typeid(Nil);
    assert( _M_info != 0 );

}

TypeInfo::TypeInfo(const std::type_info& ti)
    :
    _M_info(&ti)
{
    assert( _M_info != 0 );
}

TypeInfo::TypeInfo( TypeInfo const& ti )
    :
    _M_info( ti._M_info )
{
    assert( _M_info != 0 );
}

TypeInfo::~TypeInfo()
{
    ; // do nothing here
}

bool
TypeInfo::before(const TypeInfo& rhs) const
{
    assert( _M_info != 0 );
    return _M_info->before(*rhs._M_info);
}

const std::type_info&
TypeInfo::typeInfo() const
{
    assert( _M_info != 0 );
    return *_M_info;
}

const char*
TypeInfo::name() const
{
    assert( _M_info != 0 );
    return _M_info->name();
}

// Comparison operators


}
