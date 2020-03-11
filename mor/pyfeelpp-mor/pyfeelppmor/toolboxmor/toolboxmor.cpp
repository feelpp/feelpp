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
        .def("setOnlineAssembleMDEIM", &mor_t::setOnlineAssembleMDEIM, "set the function to assemble MDEIM for the online model", py::arg("fct"))
        .def("getDEIMReducedMesh", &mor_t::getDEIMReducedMesh, "get the reduced mesh of DEIM" )
        .def("getMDEIMReducedMesh", &mor_t::getMDEIMReducedMesh, "get the reduced mesh of MDEIM" )
        ;
    std::string modelnew_name = std::string("toolboxmor_") + std::to_string(nDim) +std::string("d");
    m.def(modelnew_name.c_str(), []() { return std::make_shared<mor_t>(); }," return a pointer on model");

    std::string crbmodelclass_name = std::string("CRBModel_") + pyclass_name;
    py::class_<CRBModel<mor_t>,std::shared_ptr<CRBModel<mor_t> > >(m, crbmodelclass_name.c_str())
        .def(py::init<>(), "init")
        ;
    std::string crbmodelnew_name = std::string("crbmodel_toolboxmor_") + std::to_string(nDim) +std::string("d");
    m.def(crbmodelnew_name.c_str(), [](std::shared_ptr<mor_t>& m) { return std::make_shared<CRBModel<mor_t> >(m); }," return a pointer on crbmodel");

    std::string crbclass_name = std::string("CRB_") + pyclass_name;
    py::class_<CRB<CRBModel<mor_t> >,std::shared_ptr<CRB<CRBModel<mor_t> > > >(m, crbclass_name.c_str())
        .def(py::init<std::string const&,
             std::shared_ptr<CRBModel<mor_t>> const&,
             crb::stage,
             std::string const&>(),
             py::arg("name"),
             py::arg("model"),
             py::arg("stage")=crb::stage::online,
             py::arg("prefixExt")=std::string(""),
             "init")
        // get rid of the return
        .def("offline", [](CRB<CRBModel<mor_t> >& c) { c.offline(); }, "run the offline step")
        ;
    std::string crbnew_name = std::string("crb_toolboxmor_") + std::to_string(nDim) +std::string("d");
    m.def(crbnew_name.c_str(), [](std::shared_ptr<CRBModel<mor_t> >& m/*, crb::stage stage = crb::stage::online*/) { return std::make_shared<CRB<CRBModel<mor_t> > >("ToolboxMor", m, crb::stage::offline); }," return a pointer on crb");
}


PYBIND11_MODULE(_toolboxmor, m )
{
    using namespace Feel;

    defToolboxMor<2>(m);
    defToolboxMor<3>(m);

    m.def("makeToolboxMorOptions", &makeToolboxMorOptions, "get options for the model" );

}

