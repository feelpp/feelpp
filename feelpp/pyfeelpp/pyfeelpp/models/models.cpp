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

#include <mpi4py/mpi4py.h>
#include<feel/feelmodels/modelproperties.hpp>

namespace py = pybind11;

using namespace Feel;


PYBIND11_MODULE(_models, m )
{
    using namespace Feel;


    std::string pyclass_name = "ModelParameters";
    py::class_<ModelParameters>(m,pyclass_name.c_str())
        .def(py::init<worldcomm_t const&>(),
             py::arg("worldComm"))
        .def("setParameterValues",&ModelParameters::setParameterValues, "set parameter values from a map of string/double pairs");

    pyclass_name = "ModelProperties";
    py::class_<ModelProperties>(m,pyclass_name.c_str())
        .def(py::init<std::string const&, std::string const&, worldcomm_t const&, std::string const&>(),"initialize ModelProperties",py::arg("filename")="",py::arg("directoryLibExpr")="",py::arg("worldComm"),py::arg("prefix")="")
        .def("parameters",static_cast<ModelParameters& (ModelProperties::*)()>(&ModelProperties::parameters), "get parameters of the model");
    
}
