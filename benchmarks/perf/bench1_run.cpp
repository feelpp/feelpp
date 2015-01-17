/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-06-14

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
   \file bench1_run.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-06-14
 */
#include <bench1.hpp>

namespace Feel
{
Bench1::Bench1( int argc,
                char** argv,
                AboutData const& ad,
                po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( soption("backend") ) ),
    meshSize( vm()["hsize"].as<double>() )
{
}

Bench1::~Bench1()
{}

void
Bench1::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    if ( this->vm().count( "nochdir" ) )
    {
        this->changeRepository( boost::format( "/benchmarks/perf/%1%/%2$dD/%3$.3f" )
                                % this->about().appName()
                                % this->vm()["dim"].as<int>()
                                % this->vm()["hsize"].as<double>() );
    }

    switch (  vm()["dim"].as<int>() )
    {
    case 1:
        run1d();
        break;

    case 2:
        run2d();
        break;

    case 3:
        run3d();
        break;

    default:
        std::cout << this->optionsDescription() << "\n";
        return;
    }
}

}
