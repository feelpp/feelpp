/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-23

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
   \file polyvisbase.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-23
 */
#include "polyvisbase.hpp"

namespace Feel
{

PolyvisBase::PolyvisBase()
{}

PolyvisBase::PolyvisBase( PolyvisBase const & )
{}

PolyvisBase::~PolyvisBase()
{}

void
PolyvisBase::init( po::variables_map const& vm )
{
    M_vm = vm;
    std::ostringstream pname;
    pname << vm["poly"].as<std::string>()
          << "("
          << vm["dim"].as<uint16_type>() << ","
          << vm["order"].as<uint16_type>() << ","
          << vm["convex"].as<std::string>()
          << ")";
    M_pname = pname.str();
}

PolyvisBase*
PolyvisBase::New( po::variables_map const& vm )
{
    std::ostringstream pname;
    pname << vm["poly"].as<std::string>()
          << "("
          << vm["dim"].as<uint16_type>() << ","
          << vm["order"].as<uint16_type>()<< ","
          << vm["convex"].as<std::string>()
          << ")";
    PolyvisBase* pbase = factory_type::instance().createObject( pname.str() );
    pbase->init( vm );
    return pbase;
}

PolyvisBase&
PolyvisBase::operator=( PolyvisBase const & o )
{
    if ( this != &o )
    {
        M_vm = o.M_vm;
        M_pname = o.M_pname;
    }

    return *this;
}

}

