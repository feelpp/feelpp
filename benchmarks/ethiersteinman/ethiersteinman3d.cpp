/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2007-04-03

  Copyright (C) 2007 EPFL

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
   \file ethiersteinman3d.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-04-03
 */
#include "ethiersteinman.hpp"

int
main( int argc, char** argv )
{
    Life::EthierSteinman ethierSteinman( argc, argv,
                                         makeAbout(), makeOptions());
    ethierSteinman.run();
}
