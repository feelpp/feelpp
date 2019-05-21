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

#include <feel/feelmodels/hdg/mixedelasticity.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <mpi4py/mpi4py.h>

namespace py = pybind11;
using namespace Feel;

template <int nDim, int Order>
void defHDGElasticity( py::module& m )
{
    typedef FeelModels::MixedElasticity<nDim, Order> mp_type;

#if 0
    auto MP = mp_type::New("mixedelasticity");
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
    using toolbox_t = MixedElasticity< nDim, Order>;
    using element_flux_t = typename toolbox_t::Vh_element_t;
    using element_potential_t = typename toolbox_t::Wh_element_t;
    using element_trace_t = typename toolbox_t::Mh_element_t;
    using element_constant_t = typename toolbox_t::Ch_element_t;
    using mesh_ptr_t = typename toolbox_t::mesh_ptrtype;
    using op_interp_ptr_t = typename toolbox_t::op_interp_ptrtype;
    using opv_interp_ptr_t = typename toolbox_t::opv_interp_ptrtype;

    std::string pyclass_name = std::string("HDGElasticity_") + std::to_string(nDim) + std::string("DP") + std::to_string(Order);
    py::class_<toolbox_t,std::shared_ptr<toolbox_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,
             //bool,
             worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             //py::arg("buildmesh")=true,
             py::arg("worldComm")=Environment::worldCommPtr(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize the heat mechanics toolbox"
             )
        //.def("init",static_cast<void (toolbox_t::*)()>(&toolbox_t::init), "initialize the HDG toolbox")
        .def("init",static_cast<void (toolbox_t::*)(mesh_ptr_t, mesh_ptr_t)>(&toolbox_t::init), "initialize the HDG toolbox",py::arg("mesh")=(mesh_ptr_t)nullptr,py::arg("mesh_visu")=(mesh_ptr_t)nullptr)
        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        // elements
        .def( "stressSpace", &toolbox_t::stressSpace, "get the stress function space")
        .def( "stressField", static_cast<element_stress_t const& (toolbox_t::*)() const>(&toolbox_t::stressField), "returns the stress field" )
        .def( "displacementSpace", &toolbox_t::displacementSpace, "get the displacement function space")
        .def( "displacementField", static_cast<element_displacement_t const& (toolbox_t::*)() const>(&toolbox_t::displacementField), "returns the displacement field" )
        // assembly
        .def("assembleAll",&toolbox_t::assembleAll, "assemble the HDG Elasticity model")
        // solve
        .def("solve",&toolbox_t::solve, "solve the HDG elasticity problem")

        //.def("exportResults",&toolbox_t::exportResults, "export the results of the heat mechanics problem",py::arg("mesh")=(mesh_ptr_t)nullptr)
        //.def("exportResults",static_cast<void (toolbox_t::*)(mesh_ptr_t, op_interp_ptr_t, opv_interp_ptr_t)>(&toolbox_t::exportResults), "export the results of the heat mechanics problem",py::arg("mesh")=(mesh_ptr_t)nullptr,py::arg("Idh")=(op_interp_ptr_t)nullptr,py::arg("Idhv")=(opv_interp_ptr_t)nullptr)
        //.def("exportResults",static_cast<void (toolbox_t::*)( double, mesh_ptr_t, op_interp_ptr_t, opv_interp_ptr_t )>(&toolbox_t::exportResults), "export the results of the heat mechanics problem", py::arg("time"),py::arg("mesh")=(mesh_ptr_t)nullptr,py::arg("Idh")=(op_interp_ptr_t)nullptr,py::arg("Idhv")=(opv_interp_ptr_t)nullptr)
        ;

}
