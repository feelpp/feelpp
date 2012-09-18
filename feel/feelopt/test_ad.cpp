/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-05-25

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file test_ad.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-05-25
 */
#include <feel/feelopt/adtype.hpp>
int main()
{
    using namespace Feel;

    ADType<double,3,2,0> x( 1. );
    ADType<double,3,2,1> y( 1. );
    ADType<double,3,2,2> z( 1. );
    ADType<double,3,2> __g = sin( M_PI*x )*cos( M_PI*y )*cos( M_PI*z );
    std::cout << "g=" << __g << "\n";
}
