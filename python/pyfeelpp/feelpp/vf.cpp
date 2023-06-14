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
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/vonmises.hpp>
#include <mpi4py/mpi4py.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Feel;

void defExpr(py::module &m)
{
    using namespace Feel;
 
}

template<int M, int N, int Order>
void addGinacMatrix( py::module& m )
{
    using namespace Feel;
    using namespace Feel::vf;
    std::string pyclass_name = fmt::format( "ExprGinacMatrix{}{}{}", M, N, Order );
    py::class_<Expr<GinacMatrix<M, N, Order>>>( m, pyclass_name.c_str() )
        .def( py::init<>() )
        .def(
            "setParameterValues", []( Expr<GinacMatrix<M, N, Order>>& e, std::pair<std::string, double> const& p )
            { return e.setParameterValues( p ); },
            "set parameter value for the expression", py::arg( "mp" ) )
        .def(
            "setParameterValues", []( Expr<GinacMatrix<M, N, Order>>& e, std::map<std::string, double /*value_type*/> const& m )
            { return e.setParameterValues( m ); },
            "set parameter value for the expression", py::arg( "mp" ) )
        .def(
            "evaluate", []( Expr<GinacMatrix<M, N, Order>>& e, bool parallel, worldcomm_ptr_t wc )
            { return e.evaluate( parallel, wc ); },
            "evaluate the expression", py::arg( "parallel" ) = true, py::arg( "worldcomm" ) = Environment::worldCommPtr() )
        .def(
            "evaluate", []( Expr<GinacMatrix<M, N, Order>>& e, std::map<std::string, double /*value_type*/> const& m )
            { return e.evaluate( m ); },
            "evaluate the expression", py::arg( "mp" ) )
        .def(
            "evaluate", []( Expr<GinacMatrix<M, N, Order>>& e, std::string const& s, Eigen::VectorXd const& x, bool parallel, worldcomm_ptr_t wc )
            {
                Eigen::VectorXd y(x.size());
                std::transform( x.begin(), x.end(),  y.begin(), 
                                  [&e,s,parallel,&wc](auto x) { 
                                      e.setParameterValues( { { s, x } } );
                                      return e.evaluate(parallel, wc)( 0, 0 ); 
                                  } );
                return y; },
            "evaluate the expression", py::arg( "parameter" ), py::arg( "values" ), py::arg( "parallel" ) = true, py::arg( "worldcomm" ) = Environment::worldCommPtr() )
        .def(
            "diff", []( Expr<GinacMatrix<M, N, Order>>& e, std::string const& s )
            { return e.template diff<1>( s ); },
            "differentiate the expression with respect to symbol s", py::arg( "symbol" ) )
        .def(
            "diff2", []( Expr<GinacMatrix<M, N, Order>>& e, std::string const& s )
            { return e.template diff<2>( s ); },
            "differentiate twice the expression with respect to symbol s", py::arg( "symbol" ) )
        .def(
            "__str__", []( Expr<GinacMatrix<M, N, Order>> const& e )
            { return str( e.expression() ); },
            "get the string representation" );
    std::string e_fn = fmt::format( "expr{}{}{}_", M, N, Order );
    m.def( e_fn.c_str(), static_cast<Expr<GinacMatrix<M, N, Order>> ( * )( std::string const&, std::string, WorldComm const&, std::string const& )>( &expr<M, N, Order> ),
           py::arg( "expr" ),
           py::arg( "filename" ) = "",
           py::arg( "worldComm" ),
           py::arg( "dir" ) = "",
           fmt::format( "create an {}x{}D expression out of a string", M, N ).c_str() );
}

PYBIND11_MODULE(_vf, m )
{
    if (import_mpi4py()<0) return ;

    addGinacMatrix<1, 1, 2>( m );
    addGinacMatrix<2, 1, 2>( m );
    addGinacMatrix<3, 1, 2>( m );
    addGinacMatrix<2, 2, 2>( m );
    addGinacMatrix<3, 3, 2>( m );

    using namespace Feel;
    using namespace Feel::vf;

    auto dimt = hana::make_tuple(2_c,3_c);
    auto ordert = hana::make_tuple( 0_c, 1_c, 2_c, 3_c );

    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,ordert)), [&m]( auto const& d )
                       {
                            constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                            constexpr int _order = std::decay_t<decltype(hana::at_c<1>(d))>::value;
                            using mesh_t = Mesh<Simplex<_dim, 1>>;
                            using mesh_ptr_t = std::shared_ptr<mesh_t>;
                            
                            m.def( "vonmises", []( Pchv_element_t<mesh_t, 1> const& d, nl::json const& model )
                                   {
                                        auto Xh = Pch<_order>(d.functionSpace()->mesh());
                                        auto r = Xh->element();
                                        if ( model["/model/type"_json_pointer] == "linear-elasticity" )
                                        {
                                            double mu = model["/model/parameters/mu"_json_pointer];
                                            double lambda = model["/model/parameters/lambda"_json_pointer];

                                            auto def = (gradv(d)+trans(gradv(d)))/2;
                                            r.on( _range=elements(d.functionSpace()->mesh()), _expr=vonmises( 2*mu*def+lambda*divv(d)*eye<_dim,_dim>() ) );
                                        }

                                        return r; 
                                    },
                                    "compute von mises stress", py::arg( "displacement" ), py::arg( "model" ) );
                        } );
}
