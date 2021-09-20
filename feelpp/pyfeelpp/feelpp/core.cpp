//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! Copyright (C) 2017-present Feel++ Consortium
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
//! @copyright 2017 Feel++ Consortiumq
//!
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/pybind11_json.hpp>
#include <feel/feelcore/repository.hpp>
#include <mpi4py/mpi4py.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "mpi_caster.hpp"
#include <boost/parameter/binding.hpp>
#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
//#include <boost/parameter/python.hpp>
#include <boost/mpl/vector.hpp>


namespace py = pybind11;

void bindRemoteData( py::module& m );

PYBIND11_MAKE_OPAQUE(Feel::worldscomm_ptr_t);

PYBIND11_MODULE(_core, m )
{
    using namespace Feel;
    namespace fs = Feel::fs;
    if (import_mpi4py()<0) return ;

    /**
     * @brief bind boost filesystem
     * 
     */
    py::class_<fs::path>(m, "Path")
        .def(py::init<std::string>())
        .def("string",&fs::path::native,"convert path to string")
        ;
    py::implicitly_convertible<std::string, boost::filesystem::path>();
    
    py::class_<po::options_description>(m,"OptionsDescription")
        .def(py::init<>())
        .def("add",static_cast<po::options_description & (po::options_description::*)(po::options_description const& )>(&po::options_description::add), "add options description",py::arg("desc"));

    py::class_<WorldComm, std::shared_ptr<WorldComm>>( m, "WorldComm" )
        .def( py::init<>() )
        .def( "isMasterRank", &Feel::WorldComm::isMasterRank, "returns true if master rank, false otherwise" )
        .def( "localRank", &Feel::WorldComm::localRank, "returns the rank of the local worldcomm" )
        .def( "globalRank", &Feel::WorldComm::globalRank, "returns the rank of the global worldcomm" )
        .def( "masterRank", &Feel::WorldComm::masterRank, "returns the master rank" )
        .def( "localComm", &Feel::WorldComm::localComm, "returns the local communicator" )
        .def( "localSize", &Feel::WorldComm::localSize, "returns the local communicator size" )
        .def( "globalComm", &Feel::WorldComm::globalComm, "returns the global communicator" )
        .def( "globalSize", &Feel::WorldComm::globalSize, "returns the global communicator size" )
        .def( "selfComm", &Feel::WorldComm::selfComm, "returns the self communicator" )
        .def( "barrier", []( std::shared_ptr<WorldComm> const& wc){
                wc->barrier();
            }, "create a barrier" )
        .def( "to_comm", []( std::shared_ptr<WorldComm>& wc ){
            return static_cast<boost::mpi::communicator>( *wc );
            }, "return worldComm as MPI  communicator" );

    py::class_<worldscomm_ptr_t>(m,"WorldsComm").def(py::init<>())
        .def("clear", &worldscomm_ptr_t::clear)
        .def("pop_back", &worldscomm_ptr_t::pop_back)
        .def("__len__", [](const worldscomm_ptr_t &v) {
        return v.size(); })
        .def("__iter__", [](worldscomm_ptr_t &v) {
                return py::make_iterator(v.begin(), v.end());
            }, py::keep_alive<0, 1>()) /* Keep vector alive while iterator is used */;
    //py::bind_vector<worldscomm_ptr_t>(m, "WorldsComm");
    m.def( "makeWorldsComm", static_cast<worldscomm_ptr_t (*)(int, worldcomm_ptr_t const& )>(&Feel::makeWorldsComm), py::arg("n")=1, py::arg("worldComm"), "create a vector of WorldComm with n entries" );
    
    py::enum_<Location>(m,"Location")
        .value("UKNOWN", Location::unknown )
        .value("GLOBAL", Location::global )
        .value("LOCAL", Location::local )
        .value("GIT", Location::git )
        .value("GIVEN", Location::given )
        .export_values();
    py::class_<Repository::Config>(m,"Config")
        .def(py::init<fs::path,Location>(),"Construct a Config",py::arg("directory"),py::arg("location"))
        .def_readwrite("directory", &Repository::Config::directory)
        .def_readwrite("location", &Repository::Config::location)
        .def_readwrite("data", &Repository::Config::data)
        ;
    m.def( "globalRepository", &Feel::globalRepository,  "create a Global Repository Config", py::arg("directory"), py::arg("data")=nl::json{} );
    m.def( "localRepository", &Feel::localRepository,  "create a Local Repository Config", py::arg("directory"), py::arg("data")=nl::json{} );
    m.def( "unknownRepository", &Feel::unknownRepository,  "create an Unknown Repository Config" );
    py::class_<Repository>(m,"Repository")
        .def(py::init<>())
        .def(py::init<Repository::Config>(),"Construct a repository from a path and location type",py::arg("config"))
        .def("root",&Feel::Repository::root, "return the root directory of the repository")
        .def("globalRoot",&Feel::Repository::globalRoot, "return the global root directory")
        .def("geo",&Feel::Repository::geo, "return the geo directory of the repository")
        .def("directory",&Feel::Repository::directory, "return the directory relative to the root of the repository")
        .def("relativeDirectory",&Feel::Repository::relativeDirectory, "return the directory relative to the root of the repository")
        .def("isLocal",&Feel::Repository::isLocal, "return true if the repository is local")
        .def("isGlobal",&Feel::Repository::isGlobal, "return true if the repository is global")
        .def("isGit",&Feel::Repository::isGit, "return true if the repository is git")
        .def("isGiven",&Feel::Repository::isGiven, "return true if the repository is given")
        .def("userName",&Feel::Repository::userName, "return the user name")
        .def("userEmail",&Feel::Repository::userEmail, "return the user email")
        ;

    py::class_<Environment>(m,"Environment")
        .def( py::init<py::list,po::options_description,Repository::Config>(),"Construct a Feel++ Environment",py::arg("arg"), py::arg("opts") = feel_nooptions(),py::arg("config")=localRepository("."))
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
        .def_static("appRepository",&Feel::Environment::appRepository,"get the application repository",py::return_value_policy::move)
        .def_static("findFile",&Feel::Environment::findFile,"find file",py::return_value_policy::move)
        .def_static("expand",&Feel::Environment::expand,"expand variable in string",py::return_value_policy::move)
        .def_static("setConfigFile",&Feel::Environment::setConfigFile,"set config file and update variable map",py::arg("filename"))
        .def_static("changeRepository", []( std::string const& fmt,  bool subdir )
            { Feel::Environment::changeRepository(_directory=boost::format(fmt),_subdir=subdir); },py::arg("directory"), py::arg("subdir")=true,"change repository")
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
    
    m.def( "feel_options", &Feel::feel_options, py::arg("prefix")="", "create feelpp options with optional prefix" );
    m.def( "feel_nooptions", &Feel::feel_nooptions, "create feelpp options with optional prefix" );
}
