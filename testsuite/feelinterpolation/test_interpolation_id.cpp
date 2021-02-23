/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 21 août 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;


int main(int argc, char** argv )
{
    Environment env( argc, argv );
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pchv<1>( mesh );
    auto Yh = Pchv<1>( mesh );
    auto Ih = I( _domainSpace = Xh, _imageSpace = Yh );
    Ih.matPtr()->printMatlab( "Ih.m");
}
