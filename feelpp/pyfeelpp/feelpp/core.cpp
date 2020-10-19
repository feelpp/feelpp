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
#include <pybind11/stl.h>

#include <feel/feel.hpp>
#include <mpi4py/mpi4py.h>

#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/binding.hpp>
//#include <boost/parameter/python.hpp>
#include <boost/mpl/vector.hpp>

#include<feel/feelcore/environment.hpp>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

void bindRemoteData( py::module& m );

PYBIND11_MAKE_OPAQUE(Feel::worldscomm_ptr_t);

PYBIND11_MODULE(_core, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;

    
    py::class_<po::options_description>(m,"OptionsDescription")
        .def(py::init<>())
        .def("add",static_cast<po::options_description & (po::options_description::*)(po::options_description const& )>(&po::options_description::add), "add options description",py::arg("desc"));

    py::class_<WorldComm,std::shared_ptr<WorldComm>>(m,"WorldComm")
        .def(py::init<>())
        .def("isMasterRank", &Feel::WorldComm::isMasterRank,"returns true if master rank, false otherwise")
        .def("localRank", &Feel::WorldComm::localRank,"returns the rank of the local worldcomm")
        .def("globalRank", &Feel::WorldComm::globalRank,"returns the rank of the global worldcomm")
        .def("masterRank", &Feel::WorldComm::masterRank,"returns the master rank")
        ;
    py::class_<worldscomm_ptr_t>(m,"WorldsComm").def(py::init<>())
        .def("clear", &worldscomm_ptr_t::clear)
        .def("pop_back", &worldscomm_ptr_t::pop_back)
        .def("__len__", [](const worldscomm_ptr_t &v) { return v.size(); })
        .def("__iter__", [](worldscomm_ptr_t &v) {
                return py::make_iterator(v.begin(), v.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */;
    //py::bind_vector<worldscomm_ptr_t>(m, "WorldsComm");
    m.def( "makeWorldsComm", static_cast<worldscomm_ptr_t (*)(int, worldcomm_ptr_t const& )>(&Feel::makeWorldsComm), py::arg("n")=1, py::arg("worldComm"), "create a vector of WorldComm with n entries" );
    
    py::class_<Environment>(m,"Environment")
        .def( py::init<py::list,po::options_description>(),"Construct a Feel++ Environment",py::arg("arg"), py::arg("opts") = feel_nooptions())
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
        .def_static("worldCommPtr",static_cast<worldcomm_ptr_t const& (*)()>(&Feel::Environment::worldCommPtr), "get the Environment WorldComm")
        .def_static("worldsComm",static_cast<worldscomm_ptr_t& (*)(int)>(&Feel::Environment::worldsComm), "get the Environment WorldComm",py::return_value_policy::copy,py::arg("size")=1)
        .def_static("worldsCommSeq",static_cast<worldscomm_ptr_t& (*)(int)>(&Feel::Environment::worldsCommSeq), "get the Environment sequential WorldsComm",py::arg("size")=1)
        .def_static("rootRepository",&Feel::Environment::rootRepository,"get the root repository for Feel++, default $HOME/feel",py::return_value_policy::move)
        .def_static("downloadsRepository",&Feel::Environment::downloadsRepository,"get the downloads repository for Feel++",py::return_value_policy::move)
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
    py::class_<std::vector<bool>>(m,"vector_bool").def(py::init<>());

    
    bindRemoteData( m );
    
}
