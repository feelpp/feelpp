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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fmt/core.h>
#include <feel/feelcore/environment.hpp>
#include <feel/feelts/tsbase.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <mpi4py/mpi4py.h>

#include <boost/shared_ptr.hpp>

namespace py = pybind11;

using namespace Feel;

template<typename SpaceT>
void defBDF( py::module& m )
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
        suffix = std::string( "Pch" );
    if ( space_t::is_continuous && space_t::is_vectorial )
        suffix = std::string( "Pchv" );
    if ( !space_t::is_continuous && space_t::is_scalar )
        suffix = std::string( "Pdh" );
    if ( !space_t::is_continuous && space_t::is_vectorial )
        suffix = std::string( "Pdhv" );
    std::string pyclass_name = fmt::format( "BDF_{}_{}D_P{}", suffix, Dim, Order );
    std::cout << fmt::format( "class name: {}", pyclass_name ) << std::endl;
    using bdf_t = Bdf<space_t>;
    using bdf_ptr_t = std::shared_ptr<bdf_t>;

    //Bdf( space_ptrtype const& space, std::string const& name, std::string const& prefix="", po::variables_map const& vm =  Environment::vm() );

    py::class_<bdf_t, TSBase, std::shared_ptr<bdf_t>>( m, pyclass_name.c_str() )
        .def( py::init<space_ptr_t const&, std::string const&, std::string const&, po::variables_map const&>(),
              py::arg( "space" ), py::arg( "name" ), py::arg( "prefix" )="", py::arg( "vm" ), "Initialize a Bdf" )
        .def ("priorTimes", &bdf_t::priorTimes, "get the times prior to time initial")
        .def ("initialize", static_cast<void (bdf_t::*)( typename bdf_t::element_type const& )>(&bdf_t::initialize), py::arg("field"), "Initialize all the entries of the unknown vector to be derived with the vector u0")
        .def ("initialize", static_cast<void (bdf_t::*)( typename bdf_t::unknowns_type const& )>(&bdf_t::initialize), py::arg("fields"), "Initialize all the entries of the unknown vector to be derived with the vector set")
        .def ("isFinished", &bdf_t::isFinished, "return if time stepping is finished, false otherwise")
        .def ("iteration", &bdf_t::iteration, "return the iteration");

    m.def(
        "bdf", []( space_ptr_t const& space, std::string const& name, std::string const& prefix )
        { return Bdf<space_t>{ space, name, prefix }; },
        py::arg( "space" ), py::arg( "name" ), py::arg( "prefix" )="", "create a BDF time stepper" );
}

class PyTSBase : public TSBase
{
public:
    using TSBase::TSBase;
    
};
PYBIND11_MODULE(_ts, m )
{
    using namespace Feel;
    using namespace hana::literals;
    if (import_mpi4py()<0) return ;

    py::class_<TSBase, PyTSBase, std::shared_ptr<TSBase>>(m,"TSBase")
        .def(py::init<>())
        .def("isFinished", &TSBase::isFinished, "return if time stepping is finished, false otherwise")
        ;
    auto ordert = hana::make_tuple( 1_c, 2_c );
    hana::for_each( ordert, [&m](auto const& o ){
        constexpr int _order = std::decay_t<decltype(o)>::value;
        // 1D
        //std::cout << fmt::format("-- BDF Pch 1D P{}", _order ) << std::endl;
        defBDF<Pch_type<Mesh<Simplex<1>>, _order>>( m );
        // 2D
        //std::cout << fmt::format("-- BDF Pch 2D P{}", _order ) << std::endl;
        defBDF<Pch_type<Mesh<Simplex<2>>, _order>>( m );
        //std::cout << fmt::format("-- BDF Pchv 2D P{}", _order ) << std::endl;
        defBDF<Pchv_type<Mesh<Simplex<2>>, _order>>( m );
        // 3D
        //std::cout << fmt::format("-- BDF Pch 3D P{}", _order ) << std::endl;
        defBDF<Pch_type<Mesh<Simplex<3>>, _order>>( m );
        defBDF<Pchv_type<Mesh<Simplex<3>>, _order>>( m );
    });

}
