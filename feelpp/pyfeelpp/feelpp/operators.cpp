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
//! @date 23 Jul 2021
//! @copyright 2017 Feel++ Consortium
//!
#include <pybind11/pybind11.h>
#include <fmt/core.h>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelvf/vf.hpp>
#include <mpi4py/mpi4py.h>
#include <pybind11/eigen.h>
#include <boost/algorithm/string.hpp>
namespace py = pybind11;

template <typename SpaceT>
void defOperator( py::module& m )
{
    using namespace Feel;
    using space_t = SpaceT;
    using space_ptr_t = std::shared_ptr<space_t>;
    using mesh_t = typename SpaceT::mesh_type;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    constexpr int Dim = mesh_t::nDim;
    constexpr int Order = space_t::basis_type::nOrder;
    constexpr int RealDim = mesh_t::nRealDim;
    using size_type = uint32_type;
    

    std::string suffix;
    if ( space_t::is_continuous && space_t::is_scalar )
        suffix = std::string("Pch");
    if ( space_t::is_continuous && space_t::is_vectorial )
        suffix = std::string("Pchv");
    if ( !space_t::is_continuous && space_t::is_scalar )
        suffix = std::string("Pdh");
    if ( !space_t::is_continuous && space_t::is_vectorial )
        suffix = std::string("Pdhv");
    std::string pyclass_name = fmt::format( "Mass_{}_{}D_P{}", suffix, Dim, Order );
    VLOG(2) << fmt::format("[pyfeelpp] class name: {}", pyclass_name ) << std::endl;
    using iterator_range_t = elements_reference_wrapper_t<mesh_t>;
    using faces_iterator_range_t = faces_reference_wrapper_t<mesh_t>;

    m.def(
        "mass", []( space_ptr_t const& domain, space_ptr_t const& image, iterator_range_t const& r, std::string const& c = "1" )
        { 
            auto a=form2(_test=image,_trial=domain);
            auto u=domain->element();
            auto v=image->element();

            a=integrate(_range=r,_expr=expr(c)*trace(trans(idt(u))*id(v)));
            return a.matrixPtr(); },
        py::arg( "trial" ), py::arg( "test" ), py::arg( "range" ), py::arg( "coeff" )=std::string{"1"}, "create a mass matrix" );
    m.def(
        "mass", []( space_ptr_t const& domain, space_ptr_t const& image, faces_iterator_range_t const& r, std::string const& c = "1" )
        { 
            auto a=form2(_test=image,_trial=domain);
            auto u=domain->element();
            auto v=image->element();

            a=integrate(_range=r,_expr=expr(c)*trace(trans(idt(u))*id(v)));
            return a.matrixPtr(); },
        py::arg( "trial" ), py::arg( "test" ), py::arg( "range" ), py::arg( "coeff" )=std::string{"1"}, "create a mass matrix" );
    m.def(
        "stiffness", []( space_ptr_t const& domain, space_ptr_t const& image, iterator_range_t const& r, std::string const& c = "1" )
        { 
            auto a=form2(_test=image,_trial=domain);
            auto u=domain->element();
            auto v=image->element();
            if ( boost::algorithm::contains(c, "{") && boost::algorithm::contains(c, "}") )
                a = integrate( _range = r, _expr = trace( expr<RealDim, RealDim>( c ) * trans( gradt( u ) ) * grad( v ) ) );
            else
                a=integrate(_range=r,_expr=trace(expr(c)*trans(gradt(u))*grad(v)));
            return a.matrixPtr(); },
        py::arg( "trial" ), py::arg( "test" ), py::arg( "range" ), py::arg( "coeff" )=std::string{"1"},"create a stiffness matrix" );
}

PYBIND11_MODULE( _operators, m )
{
    using namespace Feel;
    using namespace hana::literals;
    if ( import_mpi4py() < 0 ) return;

    auto ordert = hana::make_tuple( 1_c, 2_c );
    hana::for_each( ordert, [&m](auto const& o ){
        constexpr int _order = std::decay_t<decltype(o)>::value;
        // 1D
        //std::cout << fmt::format("-- Pch 1D P{}", _order ) << std::endl;
        defOperator<Pch_type<Mesh<Simplex<1>>, _order>>( m );
        // 2D
        //std::cout << fmt::format("-- Pch 2D P{}", _order ) << std::endl;
        defOperator<Pch_type<Mesh<Simplex<2>>, _order>>( m );
        //std::cout << fmt::format("-- Pchv 2D P{}", _order ) << std::endl;
        defOperator<Pchv_type<Mesh<Simplex<2>>, _order>>( m );
        // 3D
        //std::cout << fmt::format("-- Pch 3D P{}", _order ) << std::endl;
        defOperator<Pch_type<Mesh<Simplex<3>>, _order>>( m );
    });
    
}
