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
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <feel/feelmor/pbdw.hpp>

namespace py = pybind11;

PYBIND11_MODULE(_pbdw, m )
{
    using namespace Feel;
    using pbdw_t = PBDWOnline;
    using vectorN_type = typename pbdw_t::vectorN_type;
    py::class_<pbdw_t,std::shared_ptr<pbdw_t>>(m, "Pbdw")
        .def(py::init<std::string const&, int, std::string const&, std::string const&>(),
             py::arg("name")=std::string("pbdw"),
             py::arg("dbLoad")=1,
             py::arg("dbFilename")=std::string(""),
             py::arg("dbId")=std::string(""),
             "Initialize PBDW"
             )
        .def("dimensionN", &pbdw_t::dimensionN, "Dimension of the Reduced Basis" )
        .def("dimensionM", &pbdw_t::dimensionM, "Dimension of the sensors" )
        .def("dimension", &pbdw_t::dimension, "Dimension of PBDW" )
        .def("dimensionF", &pbdw_t::dimensionF, "Dimension of outputs" )
        .def("matrix", &pbdw_t::matrix, "Matrix of PBDW" )
        .def("matrixF", &pbdw_t::matrixF, "Matrix of outputs" )
        .def("sensorNames", &pbdw_t::sensorNames, "Name of the sensors" )
        .def("outputNames", &pbdw_t::outputNames, "Name of the outputs" )
        .def("online", py::overload_cast<vectorN_type const&>(&pbdw_t::online, py::const_), "returns the coefficients of the solution", py::arg("yobs") )
        .def("outputs", py::overload_cast<vectorN_type const&>(&pbdw_t::outputs, py::const_), "returns the outputs of the mode", py::arg("yobs") )
        .def("online", py::overload_cast<vectorN_type const&, std::vector<int> const&, bool>(&pbdw_t::online, py::const_), "returns the coefficients of the solution", py::arg("yobs"), py::arg("sensors"), py::arg("toComplete")=true )
        .def("outputs", py::overload_cast<vectorN_type const&, std::vector<int> const&, bool>(&pbdw_t::outputs, py::const_), "returns the outputs of the mode", py::arg("yobs"), py::arg("sensors"), py::arg("toComplete")=true )
        .def("online", py::overload_cast<vectorN_type const&, std::vector<std::string> const&>(&pbdw_t::online, py::const_), "returns the coefficients of the solution", py::arg("yobs"), py::arg("sensors") )
        .def("outputs", py::overload_cast<vectorN_type const&, std::vector<std::string> const&>(&pbdw_t::outputs, py::const_), "returns the outputs of the mode", py::arg("yobs"), py::arg("sensors") )
        .def("onlineWithout", py::overload_cast<vectorN_type const&, std::vector<int> const&>(&pbdw_t::onlineWithout, py::const_), "returns the coefficients of the solution", py::arg("yobs"), py::arg("sensors") )
        .def("outputsWithout", py::overload_cast<vectorN_type const&, std::vector<int> const&>(&pbdw_t::outputsWithout, py::const_), "returns the outputs of the mode", py::arg("yobs"), py::arg("sensors") )
        .def("onlineWithout", py::overload_cast<vectorN_type const&, std::vector<std::string> const&>(&pbdw_t::onlineWithout, py::const_), "returns the coefficients of the solution", py::arg("yobs"), py::arg("sensors") )
        .def("outputsWithout", py::overload_cast<vectorN_type const&, std::vector<std::string> const&>(&pbdw_t::outputsWithout, py::const_), "returns the outputs of the mode", py::arg("yobs"), py::arg("sensors") )
        ;

    m.def("makePbdwOptions", &pbdw_options, "make PBDW options" );
}
