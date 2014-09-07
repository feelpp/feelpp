/*  -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-08-26

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
   \file laplacian_dirac.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-08-26
 */

#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options().add(backend_options( "diffusion" )),
                     _about=about(_name="laplacian_dirac",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Vh = Pch<3>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    auto ts = bdf( _space=Vh );
    ts->start();
    auto l = form1( _test=Vh );
    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=idt(u)*id(u)*ts->polyDerivCoefficient( 0 ) + gradt(u)*trans(grad(v)) );
    a+=integrate(_range=boundaryfaces(mesh),
                 _expr=-gradt(u)*N()*id(v)-gradt(v)*N()*id(u)
                 +50*idt(u)*id(v)/hFace() );
    auto backend_diffusion = backend(_prefix="diffusion");
    auto prec_diffusion = preconditioner( _prefix="diffusion",_matrix=a.matrixPtr(),
                                          _pc=backend_diffusion->pcEnumType()/*LU_PRECOND*/,
                                          _pcfactormatsolverpackage=backend_diffusion->matSolverPackageEnumType(),
                                          _backend=backend_diffusion );
    auto e = exporter( _mesh=mesh );
    ts->initialize( u );
    for( ts->start(); !ts->isFinished(); ts->next() )
    {
        auto time = ts->time();
        LOG(INFO) << "Time = " << time << "s";
        if ( ts->iteration() == 1 )
        {
            BOOST_FOREACH( auto i, Vh->dof()->markerToDof( std::string("center") ) )
            {
                l.add( i.second, 1 );
            }
        }
        l += integrate( _range=elements(mesh),
                        _expr= idv(ts->polyDeriv())*id(v) );
        a.solveb(_rhs=l,_solution=u, _backend=backend_diffusion, _prec=prec_diffusion);
        ts->shiftRight( u );

        e->step(time)->add( "u", u );
        e->save();
        l.zero();
    }
    return 0;

}
