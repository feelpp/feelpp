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

#include <pybind11/embed.h>
namespace py = pybind11;

int main() {
    
    using namespace py::literals;
    
    py::scoped_interpreter guard{};
    py::module sys = py::module::import("sys");
    py::print(sys.attr("path"));
    //py::module sympy = py::module::import("sympy");
    auto locals = py::dict("name"_a="World", "number"_a=42);
    py::exec(R"(
        from sympy import *;
        x, y, z, t = symbols('x y z t');
        int_cos=integrate(cos(x), x);
        print('print:',int_cos);
        E=ccode(int_cos);
    )", py::globals(), locals);
    auto message = locals["E"].cast<std::string>();
    std::cout << message;
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
