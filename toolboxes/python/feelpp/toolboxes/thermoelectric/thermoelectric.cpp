/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2018-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include <feel/feelpython/pybind11/pybind11.h>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/thermoelectric/thermoelectric.hpp>

namespace py = pybind11;
using namespace Feel;

template<int nDim>
void defToolbox(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using model_heat_type = FeelModels::Heat< Simplex<nDim,Dynamic>,
                                              Lagrange<Dynamic, Scalar,Continuous,PointSetFekete> >;
    using model_electric_type = FeelModels::Electric< Simplex<nDim,Dynamic>,
                                                      Lagrange<Dynamic, Scalar,Continuous,PointSetFekete> >;

    using toolbox_t = ThermoElectric< model_heat_type,model_electric_type >;
    using element_thermoelectricpotential_t = typename toolbox_t::electric_model_type::element_electricpotential_type;
    using element_thermoelectricpotential_ptr_t = typename toolbox_t::electric_model_type::element_electricpotential_ptrtype;

    std::string pyclass_name = std::string("Thermoelectric_") + std::to_string(nDim) + "DDynamic";
    py::class_<toolbox_t,std::shared_ptr<toolbox_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init([](std::string const& prefix, std::string const& keyword, py::object worldComm, std::string const& subprefix, ModelBaseRepository const& modelRep, int orderTemperature, int orderPotential, int orderGeometry) {
                 worldcomm_ptr_t wc = worldComm.is_none() ? Environment::worldCommPtr() : worldComm.cast<worldcomm_ptr_t>();
                 return new toolbox_t(prefix, keyword, wc, subprefix, modelRep, RuntimeOrder{orderTemperature}, RuntimeOrder{orderPotential}, RuntimeOrder{orderGeometry});
             }),
             py::arg("prefix"),
             py::arg("keyword")=std::string("thermo-electric"),
             py::arg("worldComm")=py::none(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             py::arg("orderTemperature") = 1,
             py::arg("orderPotential") = 1,
             py::arg("orderGeometry") = 1,
             "Initialize the thermoelectric mechanics toolbox"
             )
        .def("init",&toolbox_t::init, "initialize the thermoelectric mechanics toolbox",py::arg("buildModelAlgebraicFactory")= true)

        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        .def( "heatModel", static_cast<typename toolbox_t::heat_model_ptrtype (toolbox_t::*)()>(&toolbox_t::heatModel) )
        .def( "electricModel", static_cast<typename toolbox_t::electric_model_ptrtype (toolbox_t::*)()>(&toolbox_t::electricModel) )
        //.def( "rangeMeshElements", &toolbox_t::rangeMeshElements, "get the range of mesh elements" )
#if 0
        // elements
        .def( "spaceThermoelectricPotential", &toolbox_t::spaceThermoelectricPotential, "get the potential function space")
        .def( "fieldThermoelectricPotential", static_cast<element_thermoelectricpotential_t const& (toolbox_t::*)() const>(&toolbox_t::fieldThermoelectricPotential), "returns the thermoelectric potential field" )
        .def( "fieldThermoelectricPotentialPtr", static_cast<element_thermoelectricpotential_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldThermoelectricPotentialPtr), "returns the thermoelectric potential field shared_ptr" )
        .def( "spaceThermoelectricField", &toolbox_t::spaceThermoelectricField, "get the field function space")
        .def( "fieldThermoelectricField", static_cast<element_thermoelectricfield_t const& (toolbox_t::*)() const>(&toolbox_t::fieldThermoelectricField), "returns the thermoelectric field" )
        .def( "fieldThermoelectricFieldPtr", static_cast<element_thermoelectricfield_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldThermoelectricFieldPtr), "returns the thermoelectric field shared_ptr" )
        .def( "fieldCurrentDensity", static_cast<element_thermoelectricfield_t const& (toolbox_t::*)() const>(&toolbox_t::fieldCurrentDensity), "returns the current density" )
        .def( "fieldCurrentDensityPtr", static_cast<element_thermoelectricfield_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldCurrentDensityPtr), "returns the current density shared_ptr" )
#endif
        // solve
        .def("solve",&toolbox_t::solve, "solve the thermoelectric mechanics problem, set boolean to true to update velocity and acceleration")
        .def("exportResults",static_cast<void (toolbox_t::*)()>(&toolbox_t::exportResults), "export the results of the thermoelectric mechanics problem")
        .def("exportResults",static_cast<void (toolbox_t::*)( double )>(&toolbox_t::exportResults), "export the results of the thermoelectric mechanics problem", py::arg("time"))
        ;
        
}
    

PYBIND11_MODULE(_thermoelectric, m )
{
    using namespace Feel;

    py::module::import("feelpp.toolboxes.heat");
    py::module::import("feelpp.toolboxes.electric");
    defToolbox<2>(m);
    defToolbox<3>(m);

}

