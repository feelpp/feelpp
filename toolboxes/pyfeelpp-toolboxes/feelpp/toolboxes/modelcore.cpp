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
#include <feel/feelcore/pybind11_json.hpp>
#include <pybind11/pybind11.h>

#include <feel/feelmodels/modelcore/modelbase.hpp>
#include <feel/feelmodels/modelcore/modelalgebraic.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/options.hpp>


namespace py = pybind11;
using namespace Feel;


PYBIND11_MODULE(_modelcore, m )
{
    using namespace Feel;
    using namespace Feel::FeelModels;

    std::string pyclass_name = std::string("ModelBaseRepository");
    py::class_<ModelBaseRepository>(m,pyclass_name.c_str())
        .def(py::init<std::string const&>(),py::arg("rootDirWithoutNProc")="","initialize root directory without the number of processes")
        .def("root",&ModelBaseRepository::root,"get the root directory")
        .def("rootWithNumProc",&ModelBaseRepository::rootWithNumProc,"get the root directory with the number of processes")
        .def("rootWithoutNumProc",&ModelBaseRepository::rootWithoutNumProc,"get the root directory without the number of processes")
        .def("expr",&ModelBaseRepository::expr,"get the directory where the expressions are stored");

    m.def( "modelbase_options", &modelbase_options, "get modelbase command line options");
    m.def( "modelalgebraic_options", &modelalgebraic_options, "get modelalgebraic command line options");
    m.def( "modelnumerical_options", &modelnumerical_options, "get modelnumerical command line options");
    m.def( "toolboxes_options",py::overload_cast<std::string const&>(&toolboxes_options), "get fluid command line options");
    m.def( "toolboxes_options",py::overload_cast<std::string const&,std::string const&>(&toolboxes_options), "get fluid command line options");
    m.def( "alemesh_options", &alemesh_options, "get alemesh command line options");

    py::class_<ModelBase, std::shared_ptr<ModelBase>>(m,"ModelBase")
        .def(py::init<std::string const&,worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("worldComm")=Environment::worldCommPtr(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize ModelBase base class"
             )
        .def("worldComm", static_cast<worldcomm_ptr_t& (ModelBase::*)()>(&ModelBase::worldCommPtr), "get model WorldComm")
        .def("worldsComm", static_cast<worldscomm_ptr_t& (ModelBase::*)()>(&ModelBase::worldsComm), "get model WorldsComm")

        // prefix
        .def("prefix", &ModelBase::prefix, "get model prefix")
        .def("subprefix", &ModelBase::subPrefix, "get model sub prefix")
        // root repo
        .def("repository", &ModelBase::repository, "get model repository")
        .def("rootRepository", &ModelBase::rootRepository, "get model root repository")
        // verbose
        .def("verbose", &ModelBase::verbose, "return true if model is verbose, false otherwise")
        .def("verboseAllProc", &ModelBase::verboseAllProc, "return true if model is verbose on all processes, false otherwise")
        .def("log", &ModelBase::log, "write message to log", py::arg("className"), py::arg("functionName"), py::arg("msg") )

        // info
        .def("printInfo", py::overload_cast<>(&ModelBase::printInfo,py::const_), "print model info")
        .def("saveInfo", py::overload_cast<>(&ModelBase::saveInfo,py::const_), "save model info")
        .def("printAndSaveInfo", &ModelBase::printAndSaveInfo, "print and save model info")
        ;

    py::class_<ModelAlgebraic, std::shared_ptr<ModelAlgebraic>, ModelBase>(m,"ModelAlgebraic")
        .def(py::init<std::string const&,worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("worldComm")=Environment::worldCommPtr(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize ModelAlgebraic"
             )
        ;

    py::class_<ModelNumerical, std::shared_ptr<ModelNumerical>, ModelAlgebraic>(m,"ModelNumerical")
        .def(py::init<std::string const&,worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("worldComm")=Environment::worldCommPtr(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize ModelNumerical"
             )
        .def("isStationary", &ModelNumerical::isStationary, "return if steady state model, false otherwise")
        .def("setStationary", &ModelNumerical::setStationary, "set model steady state to true or false", py::arg("state") )
        .def("doRestart", &ModelNumerical::doRestart, "return if model is restarted, false otherwise")
        .def("setRestart", &ModelNumerical::setRestart, "set model restart status", py::arg("state"))

        .def("hasModelProperties", &ModelNumerical::hasModelProperties, "returns true if model properties are defined, false otherwise")
        .def("modelProperties", &ModelNumerical::modelPropertiesPtr, "return model properties",py::return_value_policy::reference)
        .def("setModelProperties", static_cast<void (ModelNumerical::*)(std::string const&)>(&ModelNumerical::setModelProperties), "set model properties from filename",py::arg("model"))
        .def("setModelProperties", static_cast<void (ModelNumerical::*)(nl::json const&)>(&ModelNumerical::setModelProperties), "set model properties from json",py::arg("model"))
        .def("addParameterInModelProperties", &ModelNumerical::addParameterInModelProperties, "add new parameter in model properties")
        
        .def("time", &ModelNumerical::time, "get the current time")
        .def("currentTime", &ModelNumerical::currentTime, "get the current time")
        .def("updateTime", &ModelNumerical::updateTime, "update the current time", py::arg("time"))
        .def("timeInitial", &ModelNumerical::timeInitial, "get the initial time")
        .def("setTimeInitial", &ModelNumerical::setTimeInitial, "set the initial time", py::arg("time"))
        .def("timeFinal", &ModelNumerical::timeFinal, "get the final time")
        .def("setTimeFinal", &ModelNumerical::setTimeFinal, "set the final time", py::arg("time"))
        .def("timeStep", &ModelNumerical::timeStep, "get the time step")
        .def("setTimeStep", &ModelNumerical::setTimeStep, "set the time step", py::arg("step"))
        ;
    
}

