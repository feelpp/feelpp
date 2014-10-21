/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-06-20

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
/**
   \file test_prodspace.cpp
   \author Cecile Daversin <daversin@math.unistra.fr>
   \date 2014-06-20
 */
#include <feel/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcrb/parameterspace.hpp>

int main( int argc, char **argv)
{
    using namespace Feel;
    Environment env(  _argc=argc, _argv=argv );

    /*mesh*/
    typedef Simplex<2,1> entity_type;
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef Lagrange<1, Scalar> basis_type;
    typedef bases<basis_type, basis_type> prod_basis_type;

    /*spaces*/
    typedef FunctionSpace<mesh_type, bases<basis_type>, double> space_type;
    typedef FunctionSpace<mesh_type, prod_basis_type, double > prod_space_type;

    auto mesh=unitCircle();
    auto pXh = prod_space_type::New( mesh );
    auto Xh = space_type::New( mesh );

    auto U1 = pXh->element("M_U1");
    auto U2 = Xh->element("M_U2");

    std::cout << "Xh nDof = " << Xh->nDof() << std::endl;
    std::cout << "pXh nDof = " << pXh->nDof() << std::endl;

    auto myproj1 = vf::project(_space=Xh, _expr= idv(U2));
    std::cout << "classical proj ok" << std::endl;

    auto myexpr_noview = idv(U2);
    auto myproj2 = vf::project(_space=Xh, _expr= myexpr_noview);
    std::cout << "classical proj (with copy) ok" << std::endl;

    auto myproj3 = vf::project(_space=Xh, _expr= idv(U1.element<0>()));
    std::cout << "view proj (directly in expr) ok" << std::endl;

    auto myexpr_view = idv(U1.element<0>());
    auto myproj4 = vf::project(_space=Xh, _expr=myexpr_view );
    std::cout << "view proj (with copy) ok" << std::endl;

    return 0;
}
