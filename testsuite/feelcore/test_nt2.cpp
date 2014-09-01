/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
   
   This file is part of the Feel library
   
   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 31 Aug 2014
   
   Copyright (C) 2014 Feel++ Consortium
   
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
#include <nt2/table.hpp>
#include <nt2/include/functions/plus.hpp>
#include <nt2/include/functions/multiplies.hpp>
#include <nt2/include/functions/mtimes.hpp>
#include <nt2/linalg/functions/lapack/mldivide.hpp>
//nt2/linalg/functions/lapack/mldivide.hpp>
#include <nt2/include/functions/ones.hpp>
#include <nt2/include/functions/eye.hpp>


int main(int argc, char *argv[])
{
    using namespace nt2;
nt2::table<double> A,x,f;
f=ones(10,1);
A=eye(10,10);
x=mtimes(A,f);
x=mldivide(A,f);
NT2_DISPLAY(A);
NT2_DISPLAY(f);
NT2_DISPLAY(x);
    return 0;
}


