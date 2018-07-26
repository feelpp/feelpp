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
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <pybind11/pybind11.h>

#include <feel/feel.hpp>
#include <mpi4py/mpi4py.h>

#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/binding.hpp>
//#include <boost/parameter/python.hpp>
#include <boost/mpl/vector.hpp>

#include<feel/feelcore/environment.hpp>

namespace py = pybind11;
PYBIND11_MODULE(core, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;
    
    py::class_<Environment>(m,"Environment")
        .def( py::init<py::list,po::options_description>(),"Construct a Feel++ Environment")//,py::arg("arg"), py::arg("opts") = feel_nooptions())
        .def( py::init<py::list>(),"Construct a Feel++ Environment")//,py::arg("arg"), py::arg("opts") = feel_nooptions())
        .def_static("initialized",&Feel::Environment::initialized, "return true if MPI is initialized, false otherwise",py::return_value_policy::copy)
        .def_static("finalized",&Feel::Environment::finalized, "return true if MPI is finalized, false otherwise",py::return_value_policy::copy)
        .def_static("numberOfProcessors",&Feel::Environment::numberOfProcessors, "return numberOfProcessors",py::return_value_policy::copy)
        .def_static("rank",&Feel::Environment::rank, "return process rank",py::return_value_policy::copy)
        .def_static("masterRank",&Feel::Environment::masterRank, "return master rank",py::return_value_policy::copy)
        .def_static("isSequential",&Feel::Environment::isSequential, "true if process is sequential, false otherwise",py::return_value_policy::copy)
        .def_static("isParallel",&Feel::Environment::isParallel, "true if process is parallel, false otherwise",py::return_value_policy::copy)
        .def_static("isMasterRank",&Feel::Environment::isMasterRank, "true if rank is 0, false otherwise",py::return_value_policy::copy)
        .def_static("worldComm",&Feel::Environment::worldComm, "get the Environment WorldComm",py::return_value_policy::copy)
        .def_static("rootRepository",&Feel::Environment::rootRepository,"get the root repository for Feel++, default $HOME/feel",py::return_value_policy::move)
        .def_static("findFile",&Feel::Environment::findFile,"find file",py::return_value_policy::move)
        .def_static("expand",&Feel::Environment::expand,"expand variable in string",py::return_value_policy::move)
        ;

    py::class_<Info>(m,"Info")
        //.def( py::init<py::list,po::options_description>(),"Get Information about Feel++")//,py::arg("arg"), py::arg("opts") = feel_nooptions())
        //.def( py::init<py::list>(),"Get Information about Feel++")//,py::arg("arg"), py::arg("opts") = feel_nooptions())
        .def_static("prefix",&Feel::Info::prefix,"prefix directory where Feel++ is installed",py::return_value_policy::copy)
        .def_static("libdir",&Feel::Info::libdir,"directory where libraries and plugins are installed",py::return_value_policy::copy)
        .def_static("datadir",&Feel::Info::datadir,"directory where arch-independent files are installed",py::return_value_policy::copy)
        .def_static("versionMajor",&Feel::Info::versionMajor,"Feel++ major version number",py::return_value_policy::copy)
        .def_static("versionMinor",&Feel::Info::versionMinor,"Feel++ minor version number",py::return_value_policy::copy)
        .def_static("versionMicro",&Feel::Info::versionMicro,"Feel++ micro version number",py::return_value_policy::copy)
        .def_static("version",&Feel::Info::versionString,"Feel++ version string",py::return_value_policy::copy)

        ;
    py::class_<WorldComm,std::shared_ptr<WorldComm>>(m,"WorldComm")
        .def(py::init<>())
        .def("localRank", &Feel::WorldComm::localRank,"returns the rank of the local worldcomm")
        .def("globalRank", &Feel::WorldComm::localRank,"returns the rank of the global worldcomm")
        ;
    
    py::class_<po::options_description>(m,"OptionsDescription");
}
