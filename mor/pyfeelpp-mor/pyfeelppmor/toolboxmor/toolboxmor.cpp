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
//!
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
// #include <pybind11/eigen.h>

#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/toolboxmor.hpp>

namespace py = pybind11;
using namespace Feel;

template<int nDim>
void defToolboxMor(py::module &m)
{
    using namespace Feel;
    using mor_t = ToolboxMor<nDim>;

    std::string pyclass_name = std::string("ToolboxMor_") + std::to_string(nDim) + std::string("D");
    py::class_<mor_t,std::shared_ptr<mor_t>>(m,pyclass_name.c_str())
        .def(py::init<std::string const&>(),
             py::arg("prefix")=std::string(""),
             "Initialize the toolboxmor"
             )
        .def("initModel", &mor_t::initModel, "initialized the model" )
        .def("postInitModel", &mor_t::postInitModel, "finalize the init of the model")
        .def("setInitialized", &mor_t::setInitialized, "set the model to initialized", py::arg("b") )
        .def("functionSpace", &mor_t::functionSpace, "returns the function space" )
        .def("setFunctionSpaces",&mor_t::setFunctionSpaces, "set the function spaces", py::arg("Vh"))
        .def("setAssembleDEIM", static_cast<void (mor_t::*)(std::function<vector_ptrtype(ParameterSpaceX::Element const&)> const&) >(&mor_t::setAssembleDEIM), "set the function to assemble DEIM", py::arg("fct") )
        .def("setAssembleMDEIM", &mor_t::setAssembleMDEIM, "set the function to assemble MDEIM", py::arg("fct"))
        .def("setOnlineAssembleDEIM", &mor_t::setOnlineAssembleDEIM, "set the function to assemble DEIM for the online model", py::arg("fct"))
        .def("setOnlineAssembleMDEIM", &mor_t::setOnlineAssembleDEIM, "set the function to assemble MDEIM for the online model", py::arg("fct"))
        .def("getDEIMReducedMesh", &mor_t::getDEIMReducedMesh, "get the reduced mesh of DEIM" )
        .def("getMDEIMReducedMesh", &mor_t::getMDEIMReducedMesh, "get the reduced mesh of MDEIM" )
        ;

}


PYBIND11_MODULE(_toolboxmor, m )
{
    using namespace Feel;

    defToolboxMor<2>(m);
    defToolboxMor<3>(m);

    m.def("makeToolboxMorOptions", &makeToolboxMorOptions, "get options for the model" );

}

