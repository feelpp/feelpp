/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4 expandtab

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
	     Guillaume Dollé <guillaume.dolle@math.unistra.fr>
  Date: 2013-02-07

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2013 Feel++ Consortium

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
   \file myblas.cpp
   \author Vincent HUBER <vincent.huber@cemosis.fr>
   \date 2013-12-06
 */

#include <feel/feel.hpp>
using namespace Feel;

int main( int argc, char* argv[] )
{
    // initialize feel++ environment
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about( _name="myblas",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org") );

    auto mesh = unitSquare();

    auto m = mat<2,2>(cst(1.),cst(0.),cst(0.),cst(1.));
    auto v1 = vec(cst(1),cst(1));
    auto v2 = mat<2,1>(cst(1),cst(1));
    auto v3 = trans(v2); 
    auto v4 = trans(v1); 
    auto mv1 = m*v1; //OK
    auto mv2 = m*v2; //OK
    auto mv3 = m*v3; // NOT OK
    auto mv4 = m*v4; // NOT OK
    
    double int_1 = integrate( _range = elements( mesh ), _expr= trans(v1)*mv1 ).evaluate()(0,0);
    double int_2 = integrate( _range = elements( mesh ), _expr= trans(v2)*mv2 ).evaluate()(0,0);
    double int_3 = integrate( _range = elements( mesh ), _expr= trans(v3)*mv3 ).evaluate()(0,0);
    double int_4 = integrate( _range = elements( mesh ), _expr= trans(v4)*mv4 ).evaluate()(0,0);

    std::cout 
      << int_1 << "\t"
      << int_2 << "\t"
      << int_3 << "\t"
      << int_4 
      << std::endl;
} 

