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
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/ginac.hpp>
using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description EXPRoptions( "DAR options" );
    EXPRoptions.add_options()
    ( "a", po::value<double>()->default_value( 1 ), "a parameter" )
    ( "b", po::value<double>()->default_value( 2 ), "a parameter" )
    ;
    return EXPRoptions;
}


int main(int argc, char**argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="myexpression",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    //! [mesh]
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    //! [mesh]

    //! [expr]
    auto g = expr(soption(_name="functions.g"));
    std::cout << "g=" << g << std::endl;

    auto f = expr<2,1>(soption(_name="functions.f"));
    std::cout << "f=" << f << std::endl;

    double aVal = doption("a")+doption("b");
    std::map<std::string,double> myMap; myMap["aVal"]=aVal;
    auto i = expr(soption("functions.i"),myMap);
    std::cout << "i=" << i << std::endl;
    //! [expr]

    //! [grad]
    auto grad_g=grad<2>(g);
    auto grad_f=grad(f);
    std::cout << "grad(g)=" << grad_g << std::endl;
    std::cout << "grad(f)=" << grad_f << std::endl;
    //! [grad]

    //! [laplacian]
    auto laplacian_g=laplacian(g);
    std::cout << "laplacian(g)=" << laplacian_g << std::endl;

    auto laplacian_f=laplacian(f);
    std::cout << "laplacian(f)=" << laplacian_f << std::endl;

    //! [laplacian]

    //! [div]
    auto div_f=div(f);
    std::cout << "div(f)=" << div_f << std::endl;
    //! [div]

    //! [curl]
    auto curl_f=curl(f);
    std::cout << "curl(f)=" << curl_f << std::endl;
    //! [curl]

    //! [eval]
    std::cout << "Evaluation  at  (" << doption("x") << "," << doption("y") << "):" << std::endl;
    std::cout << "           g(x,y)=" << g.evaluate() << std::endl;
    std::cout << "           f(x,y)=" << f.evaluate() << std::endl;
    std::cout << "           i(x,y)=" << i.evaluate() << std::endl;
    std::cout << "Gradient:\n";
    std::cout << "     grad(g)(x,y)=" << grad_g.evaluate() << std::endl;
    std::cout << "     grad(f)(x,y)=" << grad_f.evaluate() << std::endl;
    std::cout << "Divergence:\n";
    std::cout << "      div(f)(x,y)=" << div_f.evaluate() << std::endl;
    std::cout << "Curl:\n";
    std::cout << "     curl(f)(x,y)=" << curl_f.evaluate() << std::endl;
    std::cout << "Laplacian:\n";
    std::cout << "laplacian(g)(x,y)=" << laplacian_g.evaluate() << std::endl;
    std::cout << "laplacian(f)(x,y)=" << laplacian_f.evaluate() << std::endl;
    //! [eval]
}
//! [all]
