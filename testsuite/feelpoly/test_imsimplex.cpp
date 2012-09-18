/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-06-19

  Copyright (C) 2008-2012 University Joseph Fourier (Grenoble I)

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
   \file test_imsimplex.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-06-19
 */
#include <feel/feelpoly/im.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    IMSimplex<2,1,double> im21;
    IMSimplex<2,2,double> im22;
    IMSimplex<2,4,double> im24;
    IMSimplex<2,6,double> im26;

    IMSimplex<3,1,double> im31;
    IMSimplex<3,2,double> im32;
    IMSimplex<3,4,double> im34;

    IMSimplex<2,1,double> im212;
    IMSimplex<3,1,double> im312;

}
