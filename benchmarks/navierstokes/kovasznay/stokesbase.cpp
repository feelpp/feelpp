/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2010-01-21

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file stokesbase.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-01-21
 */
#include <stokesbase.hpp>

namespace Feel
{

StokesBase::StokesBase()
{}

StokesBase::StokesBase( StokesBase const & )
{}

StokesBase::~StokesBase()
{}

void
StokesBase::init( po::variables_map const& vm )
{
    M_vm = vm;
    std::ostringstream name;
    name << "stokes("
         << vm["orderU"].as<uint16_type>() << ","
         << vm["orderP"].as<uint16_type>()
         << ")";
    M_name = name.str();
}

StokesBase*
StokesBase::New( po::variables_map const& vm )
{
    std::ostringstream name;
    name << "stokes("
         << vm["stokes-order-u"].as<uint16_type>() << ","
         << vm["stokes-order-p"].as<uint16_type>()
         << ")";
    StokesBase* pbase = factory_type::instance().createObject( name.str() );
    pbase->init( vm );
    return pbase;
}

StokesBase&
StokesBase::operator=( StokesBase const & o)
{
    if (this != &o )
        {
            M_vm = o.M_vm;
            M_name = o.M_name;
        }
    return *this;
}

}


