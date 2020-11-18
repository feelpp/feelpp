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
#include <feel/feelmodels/coefficientformpdes/coefficientformpdes.hpp>


namespace py = pybind11;
using namespace Feel;

template<int nDim>
void defSM(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using toolbox_t =  CoefficientFormPDEs< Simplex<nDim,1>,
                                            Lagrange<0,Scalar,Continuous,PointSetFekete>,
                                            Lagrange<1,Scalar,Continuous,PointSetFekete>,
                                            Lagrange<2,Scalar,Continuous,PointSetFekete>,
                                            Lagrange<1,Vectorial,Continuous,PointSetFekete> >;
    //using element_temperature_t = typename toolbox_t::element_temperature_type;
    //using element_temperature_ptr_t = typename toolbox_t::element_temperature_ptrtype;

    std::string pyclass_name = std::string("cfpdes_") + std::to_string(nDim) + std::string("D");
    py::class_<toolbox_t,std::shared_ptr<toolbox_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,std::string const&,worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("keyword")=std::string("cfpdes"),
             py::arg("worldComm")=Environment::worldCommPtr(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize the coefficient form pdes toolbox"
             )
        .def("init",&toolbox_t::init, "initialize the heat  toolbox",py::arg("buildModelAlgebraicFactory")= true)

        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        //.def( "rangeMeshElements", &toolbox_t::rangeMeshElements, "get the range of mesh elements" )

        // elements
        //.def( "spaceTemperature", &toolbox_t::spaceTemperature, "get the temperature function space")
        //.def( "fieldTemperature", static_cast<element_temperature_t const& (toolbox_t::*)() const>(&toolbox_t::fieldTemperature), "returns the temperature field" )
        //.def( "fieldTemperaturePtr", static_cast<element_temperature_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldTemperaturePtr), "returns the temperature field shared_ptr" )

        // time stepping
        .def("timeStepBase",static_cast<std::shared_ptr<TSBase> (toolbox_t::*)() const>(&toolbox_t::timeStepBase), "get time stepping base")
        .def("updateTimeStep",&toolbox_t::updateTimeStep, "update time stepping")

        // solve
        .def("solve",&toolbox_t::solve, "solve the heat mechanics problem, set boolean to true to update velocity and acceleration")
        .def("exportResults",static_cast<void (toolbox_t::*)()>(&toolbox_t::exportResults), "export the results of the heat mechanics problem")
//        .def("exportResults",static_cast<void (toolbox_t::*)( double )>(&toolbox_t::exportResults), "export the results of the heat mechanics problem", py::arg("time"))
        ;
        
}
    

PYBIND11_MODULE(_cfpdes, m )
{
    using namespace Feel;
    
    defSM<2>(m);
    defSM<3>(m);

}

