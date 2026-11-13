/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2018-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include <feel/feelpython/pybind11/pybind11.h>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/electric/electric.hpp>
#include <stdexcept>

namespace py = pybind11;
using namespace Feel;

template<int nDim>
void defSM(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using toolbox_t = Electric< Simplex<nDim,Dynamic>,
                           Lagrange<Dynamic, Scalar,Continuous,PointSetFekete> >;
    using element_electricpotential_t = typename toolbox_t::element_electricpotential_type;
    using element_electricpotential_ptr_t = typename toolbox_t::element_electricpotential_ptrtype;

    std::string pyclass_name = std::string("Electric_") + std::to_string(nDim) + "DDynamic";
    py::class_<toolbox_t,std::shared_ptr<toolbox_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init([](std::string const& prefix, std::string const& keyword, py::object worldComm, std::string const& subprefix, ModelBaseRepository const& modelRep, int orderPotential, int orderGeometry) {
                 worldcomm_ptr_t wc = worldComm.is_none() ? Environment::worldCommPtr() : worldComm.cast<worldcomm_ptr_t>();
                 return new toolbox_t(prefix, keyword, wc, subprefix, modelRep, RuntimeOrder{orderPotential}, RuntimeOrder{orderGeometry});
             }),
             py::arg("prefix"),
             py::arg("keyword")=std::string("electric"),
             py::arg("worldComm")=py::none(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             py::arg("orderPotential") = 1,
             py::arg("orderGeometry") = 1,
             "Initialize the electric mechanics toolbox"
             )
        .def("init",&toolbox_t::init, "initialize the electric mechanics toolbox",py::arg("buildModelAlgebraicFactory")= true)

        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        .def( "rangeMeshElements", &toolbox_t::rangeMeshElements, "get the range of mesh elements" )
        
        // elements
        .def( "spaceElectricPotential", &toolbox_t::spaceElectricPotential, "get the potential function space")
        .def( "fieldElectricPotential", static_cast<element_electricpotential_t const& (toolbox_t::*)() const>(&toolbox_t::fieldElectricPotential), "returns the electric potential field" )
        .def( "fieldElectricPotentialPtr", static_cast<element_electricpotential_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldElectricPotentialPtr), "returns the electric potential field shared_ptr" )
        
        // solve
        .def("solve",&toolbox_t::solve, "solve the electric mechanics problem, set boolean to true to update velocity and acceleration")
        .def("exportResults",static_cast<void (toolbox_t::*)()>(&toolbox_t::exportResults), "export the results of the electric mechanics problem")
        .def("exportResults",static_cast<void (toolbox_t::*)( double )>(&toolbox_t::exportResults), "export the results of the electric mechanics problem", py::arg("time"))
        ;
        
}

PYBIND11_MODULE(_electric, m )
{
    using namespace Feel;

    defSM<2>(m);
    defSM<3>(m);

}
