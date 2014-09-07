/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4 expandtab

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
	     Guillaume Dollé <guillaume.dolle@math.unistra.fr>
  Date: 2013-02-07

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2013-2014 Feel++ Consortium

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
// [marker]
#include <feel/feelcore/environment.hpp>

int main( int argc, char* argv[] )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="myapp",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org") );
    std::cout << "proc " << Environment::rank()
              <<" of "<< Environment::numberOfProcessors()
              << std::endl;

}
// [marker]
