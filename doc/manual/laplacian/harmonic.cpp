/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-03-13

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
   \file harmonic.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-03-13
 */
#include <feel/feel.hpp>
#include <feel/feelmesh/meshmover.hpp>

/**
   \page LaplacianHarmonic Harmonic extension of a boundary displacement
   \author Christophe Prud'homme
   \date 2013-03-13
   \tableofcontents

   <br>
   <br>

   \section LaplacianHarmonicProblem Problem

   This  program `feelpp_doc_harmonic`  shows how to solve for the harmonic
   extension of a displacement given at a boundary of a domain.
   Given \f$\eta\f$ on \f$\partial \Omega\f$, find \f$\eta^H\f$ such that

   \f[
   \begin{split}
   -\Delta \eta^H &= 0 \mbox{ in } \Omega\\
   \eta^H &=\eta \mbox{ on } \partial \Omega.
   \end{split}
   \f]

   \f$\eta^H\f$ is the \b harmonic \b extension of \f$\eta\f$ in \f$\Omega\f$,
   then we move the mesh vertices using \f$\eta^H\f$.

   \remark Note that the degrees of freedom of \f$\eta^H\f$ must coincide with
   the mesh vertices which implies that the mesh order (the order of the
   geometric transformation) should be the same as \f$\eta^H\f$. If \f$\eta^H\f$
   is piecewise polynomial of degree \f$N\f$ then the geometric transformation
   associated to the mesh should be of degree \f$N\f$ too, i.e. \c Mesh<Simplex<2,N>>.

   \section LaplacianHarmonicResults Results

   We consider the following displacement

   \f[
   \eta = ( 0, 0.08 (x+0.5) (x-1) (x^2-1) )^T
   \f]

   on the bottom boundary and \f$\eta\f$ being 0 on the remaining borders and we
   look for \f$\eta^H\f$ as a \f$P_2\f$ piecewise polynomial function in
   \f$\Omega\f$ and the mesh associated is of the same order i.e. \c Mesh<Simplex<2,2>>.

   \section LaplacianHarmonicImplementation Implementation

   The implementation is as follows
   \snippet harmonic.cpp marker1
 */
int main(int argc, char**argv )
{
    /// [marker1]
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="harmonic",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    auto mesh = createGMSHMesh( _mesh=new Mesh<Simplex<2,2> >,
                                _desc=domain( _name="harmonic",
                                              _usenames=false,
                                              _shape="hypercube",
                                              _h=option(_name="gmsh.hsize").as<double>(),
                                              _xmin=-0.5, _xmax=1,
                                              _ymin=-0.5, _ymax=1.5 ) );
    auto Vh = Pchv<2>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    auto l = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=trace(gradt(u)*trans(grad(v))) );

    // boundary conditions
    a+=on(_range=markedfaces(mesh,1), _rhs=l, _element=u,
          _expr=zero<2,1>() );
    a+=on(_range=markedfaces(mesh,3), _rhs=l, _element=u,
          _expr=zero<2,1>() );
    a+=on(_range=markedfaces(mesh,4), _rhs=l, _element=u,
          _expr=zero<2,1>() );
    a+=on(_range=markedfaces(mesh,2), _rhs=l, _element=u,
          _expr=vec(cst(0.),0.08*(Px()+0.5)*(Px()-1)*(Px()*Px()-1)));

    a.solve(_rhs=l,_solution=u);

    auto m1 = lagrangeP1(_space=Vh)->mesh();
    auto XhVisu = Pchv<1>(m1);

    auto opIVisu = opInterpolation(_domainSpace=Vh,
                                   _imageSpace=XhVisu,
                                   _type=InterpolationNonConforme(false,true,false) );
    auto uVisu = opIVisu->operator()(u);

    // exporter mesh and harmonic extension
    auto e = exporter( _mesh=m1, _name="initial" );
    e->step(0)->setMesh( m1 );
    e->step(0)->add( "uinit", uVisu );
    e->save();

    // move the mesh vertices
    meshMove( m1, uVisu );

    // export mesh after moving the vertices
    auto e1 = exporter( _mesh=m1, _name="moved" );
    e1->step(0)->setMesh( m1  );
    e1->step(0)->add( "umoved", uVisu );
    e1->save();
    /// [marker1]

#if 0
	auto e = exporter( _mesh=m1, _name="initial" );
    e->step(0)->setMesh( m1 );
    e->step(0)->add( "u", u );
    e->save();

    meshMove( mesh, u );

    auto Vh2 = Pchv<2>( mesh );
    auto m2 = lagrangeP1(_space=Vh2)->mesh();
    auto e1 = exporter( _mesh=m2, _name="moved" );
    e1->step(0)->setMesh( m2  );
    e1->step(0)->add( "u", u );
    e1->save();
#endif


    return 0;
}
