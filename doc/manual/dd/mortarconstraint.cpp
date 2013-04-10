/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-04-10

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
   \file mortarconstraint.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-04-10
 */
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="qs_laplacian",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    typedef Mesh<Simplex<2>,1> mesh1_type;
    typedef Mesh<Simplex<2>,2> mesh2_type;
    auto master_mesh = loadMesh(_mesh=new mesh1_type, _filename="master.geo");
    auto slave_mesh = loadMesh(_mesh=new mesh2_type, _filename="slave.geo");
    auto mesh = fusion::make_vector(master_mesh, slave_mesh);

    typedef meshes<mesh1_type,mesh2_type> mesh_type;
    typedef bases<Lagrange<1>,Lagrange<1>> basis_type;
    typedef FunctionSpace< mesh_type, basis_type > space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    space_ptrtype Vh = space_type::New( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    auto u1 = u.element<0>();
    auto u2 = u.element<1>();

    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh1),
                  _expr=id(v1));
    l += integrate(_range=elements(mesh2),
                   _expr=id(v2));

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh1),
                  _expr=gradt(u1)*trans(grad(u1)) );
    a += integrate(_range=elements(mesh2),
                   _expr=gradt(u2)*trans(grad(u2)) );

    a+=on(_range=boundaryfaces(mesh1,"Dirichlet"), _rhs=l, _element=u1,
          _expr=constant(0.) );
    a+=on(_range=boundaryfaces(mesh2,"Dirichlet"), _rhs=l, _element=u2,
          _expr=constant(0.) );

    auto Vh1 = Vh->functionSpace<0>();
    auto Vh2 = Vh->functionSpace<1>();
    auto Vhm = Vh1->trace( markedfaces(mesh1, "Master"));
    auto Vhs = Vh2->trace( markedfaces(mesh2, "Slave"));
    auto c_m = form2( _trial=Vhm, _test=Vhs );
    c_m = integrate(_range=elements(mesh1), _expr=idt(u1)*trans(id(u2)) );
    auto c_s = form2( _trial=Vhs, _test=Vhs );
    c_s = integrate(_range=elements(mesh2), _expr=idt(u2)*trans(id(u2)) );

    // add product and inverse (shell matrix)
    auto Q = transfer( product( inverse( c_s.matrix() ), c_m.matrix() ), "Master", "Slave" );

    // solve Q^T A Q = Q^T f, add trans()
    auto Atilde = product( transpose( Q ), product( A, Q ) );
    auto ftilde = transpose( Q ).apply( f );
    auto utilde = transpose( Q ).apply( u );

    cg( Atilde, ftilde, utilde );

    u = Q.apply( utilde );

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
    return 0;
}
