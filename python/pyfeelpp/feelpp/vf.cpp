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
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <feel/feelpython/pybind11/pybind11.h>

#include <mpi4py/mpi4py.h>
#include<feel/feelvf/expr.hpp>
#include<feel/feelvf/ginac.hpp>

namespace py = pybind11;

using namespace Feel;

void defExpr(py::module &m)
{
    using namespace Feel;
 
}
    

PYBIND11_MODULE(_vf, m )
{
    using namespace Feel;
    using namespace Feel::vf;

    if (import_mpi4py()<0) return ;

    std::string pyclass_name = "ExprGinacEx2";
    py::class_<Expr<GinacEx<2>> >(m,pyclass_name.c_str())
        .def(py::init<>());
    pyclass_name = "ExprGinacMatrix212";
    py::class_<Expr<GinacMatrix<2,1,2>>>( m, pyclass_name.c_str() )
        .def( py::init<>() );
    pyclass_name = "ExprGinacMatrix312";
    py::class_<Expr<GinacMatrix<3,1,2>>>( m, pyclass_name.c_str() )
        .def( py::init<>() );

    m.def( "expr_", static_cast<Expr<GinacEx<2>> (*)( std::string const&, std::string const&, WorldComm const&, std::string const&)>(&expr),
           py::arg("expr"),
           py::arg("filename")="",
           py::arg("worldComm"),
           py::arg("dir")="",
           "create an expression out of a string" );
    m.def( "exprv2_", static_cast<Expr<GinacMatrix<2,1,2>> ( * )( std::string const&, std::string , WorldComm const&, std::string const& )>( &expr<2,1,2> ),
           py::arg( "expr" ),
           py::arg( "filename" ) = "",
           py::arg( "worldComm" ),
           py::arg( "dir" ) = "",
           "create an 2D vectorial expression out of a string" );
    m.def( "exprv3_", static_cast<Expr<GinacMatrix<3, 1, 2>> ( * )( std::string const&, std::string , WorldComm const&, std::string const& )>( &expr<3,1,2> ),
           py::arg( "expr" ),
           py::arg( "filename" ) = "",
           py::arg( "worldComm" ),
           py::arg( "dir" ) = "",
           "create an 3D vectorial expression out of a string" );
}
