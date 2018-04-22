/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-05-02

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file magnetostatic_main2.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-05-02
 */
#include <iostream>

#include <feel/feelmodels/hifimagnet/MagnetostaticModel/magnetostatic.hpp>


int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=MagnetostaticOptions()
                              .add(MagnetostaticOptionsLib())
                              .add(HifiMagnetOptions()),
                     _about=makeAboutMagnetostatic<FEELPP_DIM,FEELPP_ORDER,FEELPP_G_ORDER>() );

    typedef magnetostatic<FEELPP_DIM,FEELPP_ORDER,FEELPP_G_ORDER> magnetostatic_type;

    magnetostatic_type model = magnetostatic_type();
    model.init();
    model.solve();
    model.exportResults();
}
