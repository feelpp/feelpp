/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2011-06-21

   Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 \file testload.cpp
 \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 \date 16-06-2011
 */
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/integrate.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="loadmesh",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );


    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);

    LOG(INFO) << "mesh " << option(_name="gmsh.filename").as<std::string>() << " loaded";

    LOG(INFO) << "volume =" << integrate( elements( mesh ), cst( 1. ) ).evaluate();
    LOG(INFO) << "surface = " << integrate( boundaryfaces( mesh ), cst( 1. ) ).evaluate();
}
