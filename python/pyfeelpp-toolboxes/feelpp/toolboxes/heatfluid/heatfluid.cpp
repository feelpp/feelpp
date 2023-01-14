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
#include <feel/feelpython/pybind11/pybind11.h>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/heatfluid/heatfluid.hpp>

namespace py = pybind11;
using namespace Feel;

template<int nDim, int OrderT, int OrderV, int OrderP>
void defToolbox(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;

    using model_heat_type = FeelModels::Heat< Simplex<nDim,1>,
                                              Lagrange<OrderT, Scalar,Continuous,PointSetFekete> >;
    using model_fluid_type = FeelModels::FluidMechanics< Simplex<nDim,1>,
                                                        Lagrange<OrderV, Vectorial,Continuous,PointSetFekete>,
                                                        Lagrange<OrderP, Scalar,Continuous,PointSetFekete> >;

    using toolbox_t = FeelModels::HeatFluid< model_heat_type,model_fluid_type>;
    using toolbox_ptr_t = std::shared_ptr<toolbox_t>;

    using space_temperature_t = typename toolbox_t::heat_model_type::space_temperature_type;
    using element_temperature_t = typename toolbox_t::heat_model_type::element_temperature_type;
    using element_temperature_ptr_t = typename toolbox_t::heat_model_type::element_temperature_ptrtype;
//    using element_heatfluidpotential_ptr_t = typename toolbox_t::electric_model_type::element_electricpotential_ptrtype;

    std::string pyclass_name = fmt::format("HeatFluid_{}D_P{}_P{}P{}",nDim,OrderT,OrderV,OrderP);

    py::class_<toolbox_t,std::shared_ptr<toolbox_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,std::string const&,worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("keyword")=std::string("thermo-electric"),
             py::arg("worldComm")=Environment::worldCommPtr(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize the heatfluid mechanics toolbox"
             )
        .def("init",&toolbox_t::init, "initialize the heatfluid mechanics toolbox",py::arg("buildModelAlgebraicFactory")= true)

        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        .def( "setMesh", &toolbox_t::setMesh, "set the mesh", py::arg( "mesh" ) )
        .def( "updateParameterValues", &toolbox_t::updateParameterValues, "update parameter values" )
        //.def( "rangeMeshElements", &toolbox_t::rangeMeshElements, "get the range of mesh elements" )

        // temperature space and field
        .def( "modelHeat", []( toolbox_ptr_t& t ) { return t->heatModel(); } , "get the fluid model" )
        .def( "spaceTemperature", []( toolbox_ptr_t& t ) { return t->heatModel()->spaceTemperature(); } , "get the temperature function space")
        .def( "fieldTemperature", []( toolbox_ptr_t& t ) { return t->heatModel()->fieldTemperature(); } , "get the temperature function space")
        .def( "fieldTemperaturePtr", []( toolbox_ptr_t& t ) { return t->heatModel()->fieldTemperaturePtr(); }, "returns the temperature field shared_ptr" )

        // fluid space and fields
        .def( "modelFluid", []( toolbox_ptr_t& t ) { return t->fluidModel(); } , "get the fluid model" )
        .def( "spaceVelocity", []( toolbox_ptr_t& t ) { return t->fluidModel()->functionSpaceVelocity(); } , "get the velocity function space" )
        .def( "spacePressure", []( toolbox_ptr_t& t ) { return t->fluidModel()->functionSpacePressure(); } , "get the pressure function space" )
        .def( "fieldVelocity", []( toolbox_ptr_t& t ) { return t->fluidModel()->fieldVelocity(); } , "get the velocity field" )
        .def( "fieldPressure", []( toolbox_ptr_t& t ) { return t->fluidModel()->fieldPressure(); } , "get the pressure field" )

        // solve
        .def("solve",&toolbox_t::solve, "solve the heatfluid mechanics problem, set boolean to true to update velocity and acceleration")
        .def("exportResults",static_cast<void (toolbox_t::*)()>(&toolbox_t::exportResults), "export the results of the heatfluid mechanics problem")
        .def("exportResults",static_cast<void (toolbox_t::*)( double )>(&toolbox_t::exportResults), "export the results of the heatfluid mechanics problem", py::arg("time"))
        ;

}
    

PYBIND11_MODULE(_heatfluid, m )
{
    using namespace Feel;

    defToolbox<2,1,1,1>(m);
    defToolbox<2,1,2,1>(m);
    defToolbox<3,1,1,1>(m);
    defToolbox<3,1,2,1>(m);

}

