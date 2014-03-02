/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-03-01

  Copyright (C) 2014 Feel++ Consortium

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
   \file exprbase.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-04-21
 */
#include <typeinfo>
#include <feel/feelvf/exprbase.hpp>

namespace Feel
{
ExprBase::ExprBase()
{}

ExprBase::~ExprBase()
{}

std::string
ExprBase::typeName() const
{
    return typeid( *this ).name();
}

std::string
ExprBase::toString() const
{
    std::ostringstream oss;
    this->toText( oss, false );
    return oss.str();
}

} // Feel
