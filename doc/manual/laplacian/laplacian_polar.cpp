/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
             Adnane Hamiaz <hamiaz@math.unistra.fr>
       Date: 2013-04-04

  Copyright (C) 2013 Université de Strasbourg

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
   \file laplacian_polar.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \author Adnane Hamiaz <hamiaz@math.unistra.fr>
   \date 2013-04-04
 */
#include <feel/feel.hpp>

namespace Feel
{
/**
   \page LaplacianCoordinateSystem Laplacian in polar and cartesian coordinate systems
   \author Feel++ Consortium

   <br>
   <br>

   Two examples are available to solve the Laplacian in two coordinate systems (cartesian and polar) having the same solution.

   - `feelpp_doc_laplacian_polar` polar coodinate system
   - `feelpp_doc_laplacian_cartesian` cartesian coodinate system

   The default arguments(and configuration) should work seamlessly.

   \section LaplacianCoordinateSystem_Implementation

   the implementation is available in \ref doc/manual/laplacian/laplacian_polar.cpp
   \snippet laplacian_polar.cpp marker1

 */
}
int
main(int argc, char**argv )
{
    /// [marker1]
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name=
#if defined(FEELPP_POLAR)
                                  "laplacian_polar",
#else
                                  "laplacian_cartesian",
#endif
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
#if FEELPP_POLAR
    typedef Mesh<Hypercube<2> > mesh_type;
#else
    typedef Mesh<Simplex<2> > mesh_type;
#endif

    auto mesh =  loadMesh( _mesh=new mesh_type );

    auto Vh = Pch<1>( mesh );
    LOG(INFO) << "dim Vh : " << Vh->nDof();

    auto u = Vh->element();
    auto v = Vh->element();

#if FEELPP_POLAR
    auto eta1=Px();
    auto eta2=Py();

    auto Jacmat=mat<2,2> (cos(eta2),-eta1*sin(eta2),sin(eta2),eta1*cos(eta2));
    auto Jac=det(Jacmat);// it is eta1 but compute the det() anyway for the sake of demonstration;
    //auto Jac=eta1;
#else
    auto eta1=sqrt(Px()*Px()+Py()*Py());
    auto eta2=acos(Px()/eta1);
    auto Jacmat=mat<2,2> (cst(1.),cst(0.),cst(0.), cst(1.));
    auto Jac=det(Jacmat);
#endif
    auto Jacmatinv=inv(Jacmat);
    auto Jacmattrans=trans(Jacmat);
    auto Jacmattransinv=inv(Jacmattrans);

    auto uexact = pow(eta1,4./3.)*sin(4*eta2/3);

    auto l = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=trans(Jacmattransinv*trans(gradt(u)))*(Jacmattransinv*trans(grad(v)))*Jac );
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=uexact );

    a.solve(_rhs=l,_solution=u);

    LOG(INFO) << "Error L2 norm: " << normL2( _range=elements(mesh), _expr=uexact-idv(u));

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
    /// [marker1]
}
