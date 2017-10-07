//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 29 Aug 2017
//! @copyright 2017 Feel++ Consortium
//!
#if 1
#include <iostream>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelio.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelpython/pyexpr.hpp>

namespace py = pybind11;

int main(int argc, char** argv)
{

    // tag::env[]
    using namespace Feel;
    using Feel::cout;

	Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="qs_ginac",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    // end::env[]
    using namespace py::literals;
    using Feel::cout;
    auto dict = Feel::pyexprFromFile( "elasticity", {"displ", "strain", "stress", "stressn", "f", "c1", "c2"}  );
    cout << "displ : " << dict.at("displ") << std::endl;
    cout << "strain : " << dict.at("strain") << std::endl;
    cout << "stress : " << dict.at("stress") << std::endl;
    cout << "f : " << dict.at("f") << std::endl;
    cout << "c1 : " << dict.at("c1") << std::endl;
    cout << "c2 : " << dict.at("c2") << std::endl;
    
}
#else


// tag::global[]
#include <feel/feelcore/environment.hpp>
#include <feel/feelvf/vf.hpp>

int main(int argc, char**argv )
{
    // tag::env[]
    using namespace Feel;
    using Feel::cout;

	Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="qs_ginac",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    // end::env[]
    constexpr int Dim = FEELPP_DIM;
    auto p = soption("functions.p");
    auto k = soption("functions.k");
    cout << "k : " << k << std::endl;

    GinacExpr<>  pex = expr(p);
    cout << "p : " << pex << std::endl;
    
    GinacMatrixExpr<1,Dim> gradp = grad<Dim>(pex);
    cout << "gradp : " << gradp << std::endl;
    
    GinacMatrixExpr<1,Dim> fluxp = -expr(k)*grad<Dim>(pex);
    // fluxp = -expr(k)*grad<Dim>(pex);
    //cout << "fluxp : " << fluxp << std::endl;
    
    //GinacMatrixExpr<1,1> div_fluxp = div( fluxp );
    //cout << "div_fluxp : " << div_fluxp << std::endl;
    
    return 0;

}
// end::global[]

#endif
