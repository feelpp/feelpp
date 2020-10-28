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

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/elasticity/elasticity.hpp>

namespace py = pybind11;
using namespace Feel;

template<int nDim, int Order>
void defSM(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using toolbox_t = Elasticity< Simplex<nDim,1>,
                           Lagrange<Order, Scalar,Continuous,PointSetFekete> >;
    using element_temperature_t = typename toolbox_t::element_temperature_type;
    using element_temperature_ptr_t = typename toolbox_t::element_temperature_ptrtype;
    using element_velocityconvection_t = typename toolbox_t::element_velocityconvection_type;
    using element_velocityconvection_ptr_t = typename toolbox_t::element_velocityconvection_ptrtype;

    std::string pyclass_name = std::string("Elasticity_") + std::to_string(nDim) + std::string("DP") + std::to_string(Order);
    py::class_<toolbox_t,std::shared_ptr<toolbox_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,bool,worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("buildmesh")=true,
             py::arg("worldComm")=Environment::worldCommPtr(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize the elasticity mechanics toolbox"
             )
        .def("init",&toolbox_t::init, "initialize the elasticity  toolbox",py::arg("buildModelAlgebraicFactory")= true)

        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        .def( "rangeMeshElements", &toolbox_t::rangeMeshElements, "get the range of mesh elements" )
        
        // elements
        .def( "spaceTemperature", &toolbox_t::spaceTemperature, "get the temperature function space")
        .def( "fieldTemperature", static_cast<element_temperature_t const& (toolbox_t::*)() const>(&toolbox_t::fieldTemperature), "returns the temperature field" )
        .def( "fieldTemperaturePtr", static_cast<element_temperature_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldTemperaturePtr), "returns the temperature field shared_ptr" )
        //.def( "spaceVelocityConvection", &toolbox_t::spaceVelocityConvection, "get the field function space")
        .def( "fieldVelocityConvection", static_cast<element_velocityconvection_t const& (toolbox_t::*)() const>(&toolbox_t::fieldVelocityConvection), "returns the convection velocity field" )
        .def( "fieldVelocityConvectionPtr", static_cast<element_velocityconvection_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldVelocityConvectionPtr), "returns the convection velocity field shared_ptr" )
        
        // solve
        .def("solve",&toolbox_t::solve, "solve the elasticity mechanics problem, set boolean to true to update velocity and acceleration")
        .def("exportResults",static_cast<void (toolbox_t::*)()>(&toolbox_t::exportResults), "export the results of the elasticity mechanics problem")
        .def("exportResults",static_cast<void (toolbox_t::*)( double )>(&toolbox_t::exportResults), "export the results of the elasticity mechanics problem", py::arg("time"))
        ;
        
}
    
