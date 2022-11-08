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

#include <feel/feelpython/pybind11/eigen.h>
#include <feel/feelpython/pybind11/json.h>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/ginac.hpp>
#include <mpi4py/mpi4py.h>
#include <pybind11/stl.h>

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
    py::class_<Expr<GinacEx<2>>>( m, pyclass_name.c_str() )
        .def( py::init<>() )
        .def( "setParameterValues", []( Expr<GinacEx<2>>& e, std::pair<std::string, double> const& p ){
                return e.setParameterValues( p );
            }, "set parameter value for the expression", py::arg( "mp" ) )
        .def( "setParameterValues", []( Expr<GinacEx<2>>& e, std::map<std::string,double/*value_type*/> const& m ){
                return e.setParameterValues( m );
            }, "set parameter value for the expression", py::arg( "mp" ) )
        .def( "evaluate", [](  Expr<GinacEx<2>>& e, bool parallel, worldcomm_ptr_t wc ) {
                return e.evaluate( parallel, wc );
            }, "evaluate the expression", py::arg( "parallel" ) = true, py::arg( "worldcomm" ) = Environment::worldCommPtr() )
        .def( "evaluate", [](  Expr<GinacEx<2>>& e, std::map<std::string,double/*value_type*/> const& m ) {
                return e.evaluate( m );
            }, "evaluate the expression", py::arg( "mp" ) )
        .def( "evaluate", [](  Expr<GinacEx<2>>& e, std::string const& s, Eigen::VectorXd const& x, bool parallel, worldcomm_ptr_t wc ) {
                Eigen::VectorXd y(x.size());
                std::transform( x.begin(), x.end(),  y.begin(), 
                                  [&e,s,parallel,&wc](auto x) { 
                                      e.setParameterValues( { { s, x } } );
                                      return e.evaluate(parallel, wc)( 0, 0 ); 
                                  } );
                return y;
            }, "evaluate the expression", py::arg( "parameter" ), py::arg( "values" ), py::arg( "parallel" ) = true, py::arg( "worldcomm" ) = Environment::worldCommPtr() )
        .def("diff", []( Expr<GinacEx<2>>& e, std::string const& s ) { return e.diff<1>( s ); }, "differentiate the expression with respect to symbol s", py::arg( "symbol" ) )
        .def("diff2", []( Expr<GinacEx<2>>& e, std::string const& s ) { return e.diff<2>( s ); }, "differentiate twice the expression with respect to symbol s", py::arg( "symbol" ) )
        .def("__str__", []( Expr<GinacEx<2>> const& e ) { return str(e.expression()); }, "get the string representation" )
        ;
    pyclass_name = "ExprGinacMatrix212";
    py::class_<Expr<GinacMatrix<2, 1, 2>>>( m, pyclass_name.c_str() )
        .def( py::init<>() )
        .def(
            "setParameterValues", []( Expr<GinacMatrix<2,1,2>>& e, std::pair<std::string, double> const& p )
            { return e.setParameterValues( p ); },
            "set parameter value for the expression", py::arg( "mp" ) )
        .def(
            "setParameterValues", []( Expr<GinacMatrix<2,1,2>>& e, std::map<std::string, double /*value_type*/> const& m )
            { return e.setParameterValues( m ); },
            "set parameter value for the expression", py::arg( "mp" ) )
        .def(
            "evaluate", []( Expr<GinacMatrix<2,1,2>>& e, bool parallel, worldcomm_ptr_t wc )
            { return e.evaluate( parallel, wc ); },
            "evaluate the expression", py::arg( "parallel" ) = true, py::arg( "worldcomm" ) = Environment::worldCommPtr() )
        .def(
            "evaluate", []( Expr<GinacMatrix<2,1,2>>& e, std::map<std::string, double /*value_type*/> const& m )
            { return e.evaluate( m ); },
            "evaluate the expression", py::arg( "mp" ) )
        ;
    pyclass_name = "ExprGinacMatrix312";
    py::class_<Expr<GinacMatrix<3, 1, 2>>>( m, pyclass_name.c_str() )
        .def( py::init<>() )
        .def(
            "setParameterValues", []( Expr<GinacMatrix<3, 1, 2>>& e, std::pair<std::string, double> const& p )
            { return e.setParameterValues( p ); },
            "set parameter value for the expression", py::arg( "mp" ) )
        .def(
            "setParameterValues", []( Expr<GinacMatrix<3, 1, 2>>& e, std::map<std::string, double /*value_type*/> const& m )
            { return e.setParameterValues( m ); },
            "set parameter value for the expression", py::arg( "mp" ) )
        .def(
            "evaluate", []( Expr<GinacMatrix<3, 1, 2>>& e, bool parallel, worldcomm_ptr_t wc )
            { return e.evaluate( parallel, wc ); },
            "evaluate the expression", py::arg( "parallel" ) = true, py::arg( "worldcomm" ) = Environment::worldCommPtr() )
        .def(
            "evaluate", []( Expr<GinacMatrix<3, 1, 2>>& e, std::map<std::string, double /*value_type*/> const& m )
            { return e.evaluate( m ); },
            "evaluate the expression", py::arg( "mp" ) )
        ;

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
