/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-06-20

  Copyright (C) 2006 EPFL

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
   \file ale.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2006-06-20
 */
#include "ale.hpp"

int
main( int argc, char** argv )
{
    //MPI_Init(&argc, &argv);

    Feel::TestALE<1> mesh1( argc, argv, makeAbout(), makeOptions() );
    mesh1.run();
    /*
    Feel::TestALE<2> mesh2( argc, argv, makeAbout(), makeOptions());
    mesh2.run();
    */
    /*
    Feel::TestALE<3> mesh3( argc, argv, makeAbout(), makeOptions());
    mesh3.run();

    Feel::TestALE<4> mesh4( argc, argv, makeAbout(), makeOptions());
    mesh4.run();
    */
}




