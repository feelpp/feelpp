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

#include <vector>
#include <Eigen/StdVector>

//! [all]
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/function.hpp>
using namespace Feel;

namespace Feel {
struct MyFunctor
{
    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef Feel::uint16_type uint16_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 1;
    static const bool imIsPoly = true;
    MyFunctor( int i ): coord(i){}
    double operator()( uint16_type a, uint16_type b, ublas::vector<double> const& x, ublas::vector<double> const& n ) const
        {
            return x[coord];
        }
    double operator()( uint16_type a, uint16_type b, ublas::vector<double> const& x, eigen_vector_type<3,value_type> const& n ) const
        {
            return x[coord];
        }
    int coord = 0;
};
}
int main(int argc, char**argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="myfunctor",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);

    auto Xh= Pch<3>(mesh); 

    MyFunctor My_X(0);

    auto u = Xh->element();//idf(My_X)); // Will contain x

    u.on( _range=elements(mesh), _expr=idf(My_X));

    auto ex1 = exporter(_mesh=mesh,_name="ex1");
    ex1->add("u",u);
    ex1->save();
}
//! [all]
