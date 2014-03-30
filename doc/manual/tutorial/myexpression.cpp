/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

   This file is part of the Feel library

   Author(s): Vincent HUBER <vincent.huber@cemosis.fr>

   Date 2013-02-18

   Copyright (C) 2013 Université de Strasbourg

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
//! [all]
#include <feel/feel.hpp>
using namespace Feel;

int main(int argc, char**argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="myexpression",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    //! [mesh]
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    //! [mesh]

    //! [expr]
    auto g = soption(_name="functions.g");
    auto f = soption(_name="functions.f");
    //! [expr]

    //! [export]
    auto e = exporter(_mesh=mesh);
    e->add( "g", g );
    e->add( "f", f );
    e->save();
    //! [export]
}
//! [all]
