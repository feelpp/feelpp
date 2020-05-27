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
//! @date 25 Jul 2018
//! @copyright 2018 Feel++ Consortium
//!
#include <pybind11/pybind11.h>

#include <feel/feelmodels/hdg/mixedpoisson.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>


namespace py = pybind11;
using namespace Feel;

template <int nDim, int Order>
void defHDGPoisson( py::module& m )
{
    typedef FeelModels::MixedPoisson<nDim, Order> mp_type;

#if 0
    auto MP = mp_type::New("mixedpoisson");
    auto mesh = loadMesh( _mesh=new typename mp_type::mesh_type );
    decltype( IPtr( _domainSpace=Pdh<OrderT>(mesh), _imageSpace=Pdh<OrderT>(mesh) ) ) Idh ;
    decltype( IPtr( _domainSpace=Pdhv<OrderT>(mesh), _imageSpace=Pdhv<OrderT>(mesh) ) ) Idhv;
    if ( soption( "gmsh.submesh" ).empty() )
        MP -> init(mesh);
    else
    {
        Feel::cout << "Using submesh: " << soption("gmsh.submesh") << std::endl;
        auto cmesh = createSubmesh( mesh, markedelements(mesh,soption("gmsh.submesh")), Environment::worldComm() );
        Idh = IPtr( _domainSpace=Pdh<OrderT>(cmesh), _imageSpace=Pdh<OrderT>(mesh) );
        Idhv = IPtr( _domainSpace=Pdhv<OrderT>(cmesh), _imageSpace=Pdhv<OrderT>(mesh) );
        MP -> init( cmesh, mesh );
    }
#endif
    using namespace Feel;
    using namespace Feel::FeelModels;
    using toolbox_t = MixedPoisson< nDim, Order>;
    using element_flux_t = typename toolbox_t::Vh_element_t;
    using element_potential_t = typename toolbox_t::Wh_element_t;
    using element_trace_t = typename toolbox_t::Mh_element_t;
    using element_constant_t = typename toolbox_t::Ch_element_t;
    using mesh_ptr_t = typename toolbox_t::mesh_ptrtype;
    using op_interp_ptr_t = typename toolbox_t::op_interp_ptrtype;
    using opv_interp_ptr_t = typename toolbox_t::opv_interp_ptrtype;

    py::enum_<MixedPoissonPhysics>(m,"MixedPoissonPhysics")
        .value("none", MixedPoissonPhysics::None)
        .value("electric", MixedPoissonPhysics::Electric)
        .value("heat", MixedPoissonPhysics::Heat);

    std::string pyclass_name = std::string("HDGPoisson_") + std::to_string(nDim) + std::string("DP") + std::to_string(Order);
    py::class_<toolbox_t,std::shared_ptr<toolbox_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,
             MixedPoissonPhysics const&,
             //bool,
             worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("physic")=MixedPoissonPhysics::None,
             //py::arg("buildmesh")=true,
             py::arg("worldComm"),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize the heat mechanics toolbox"
             )
        //.def("init",static_cast<void (toolbox_t::*)()>(&toolbox_t::init), "initialize the HDG toolbox")
        .def("init",static_cast<void (toolbox_t::*)(mesh_ptr_t, mesh_ptr_t)>(&toolbox_t::init), "initialize the HDG toolbox",py::arg("mesh")=(mesh_ptr_t)nullptr,py::arg("mesh_visu")=(mesh_ptr_t)nullptr)
        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        // elements
        .def( "fluxSpace", &toolbox_t::fluxSpace, "get the flux function space")
        .def( "fluxField", static_cast<element_flux_t const& (toolbox_t::*)() const>(&toolbox_t::fluxField), "returns the flux field" )
        .def( "potentialSpace", &toolbox_t::potentialSpace, "get the potential function space")
        .def( "potentialField", static_cast<element_potential_t const& (toolbox_t::*)() const>(&toolbox_t::potentialField), "returns the potential field" )
        .def( "traceSpace", &toolbox_t::traceSpace, "get the trace function space")
        //.def( "fieldTrace", static_cast<element_trace_t const& (toolbox_t::*)() const>(&toolbox_t::traceField), "returns the trace field" )
        .def( "constantSpace", &toolbox_t::constantSpace, "get the constant function space")
        // assembly
        .def("assembleAll",&toolbox_t::assembleAll, "assemble the HDG Poisson model")
        // solve
        .def("solve",&toolbox_t::solve, "solve the HDG poisson problem")

        //.def("exportResults",&toolbox_t::exportResults, "export the results of the heat mechanics problem",py::arg("mesh")=(mesh_ptr_t)nullptr)
        .def("exportResults",static_cast<void (toolbox_t::*)(mesh_ptr_t, op_interp_ptr_t, opv_interp_ptr_t)>(&toolbox_t::exportResults), "export the results of the heat mechanics problem",py::arg("mesh")=(mesh_ptr_t)nullptr,py::arg("Idh")=(op_interp_ptr_t)nullptr,py::arg("Idhv")=(opv_interp_ptr_t)nullptr)
        //.def("exportResults",static_cast<void (toolbox_t::*)( double, mesh_ptr_t, op_interp_ptr_t, opv_interp_ptr_t )>(&toolbox_t::exportResults), "export the results of the heat mechanics problem", py::arg("time"),py::arg("mesh")=(mesh_ptr_t)nullptr,py::arg("Idh")=(op_interp_ptr_t)nullptr,py::arg("Idhv")=(opv_interp_ptr_t)nullptr)
        ;

}
