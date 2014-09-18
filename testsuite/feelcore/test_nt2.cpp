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
#include <nt2/include/functions/lu.hpp>
#include <nt2/include/functions/randn.hpp>
#include <nt2/include/functions/sum.hpp>
#include <nt2/include/functions/numel.hpp>
#include <nt2/include/functions/sqr.hpp>
//#include <nt2/linalg/functions/lapack/mldivide.hpp>
//nt2/linalg/functions/lapack/mldivide.hpp>
#include <nt2/include/functions/ones.hpp>
#include <nt2/include/functions/eye.hpp>


int main(int argc, char *argv[])
{
    using namespace nt2;
    table < double > A1 = _ (1. ,1000.) ;
    table < double > A2 = A1 + randn ( size ( A1 ));
    table < double > X = lu ( mtimes ( A1 , trans ( A1 ) ) );
//double rms = sqrt ( sum ( sqr ( A1 ( _) - A2 (_)) ) / numel ( A1 ) );
    NT2_DISPLAY(X);
    return 0;
}
