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
//! @copyright 2017 Feel++ Consortium
//!
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <mpi4py/mpi4py.h>

#include <boost/parameter/keyword.hpp>
#include <boost/parameter/preprocessor.hpp>
#include <boost/parameter/binding.hpp>
//#include <boost/parameter/python.hpp>
#include <boost/mpl/vector.hpp>

#include<feel/feelalg/backend.hpp>
#include<feel/feelalg/vectorublas.hpp>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

PYBIND11_MODULE(_alg, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;

    m.def( "backend_options", &Feel::backend_options, py::arg("prefix"), "create a backend options descriptions with prefix" );
    py::class_<VectorUblas<double>,std::shared_ptr<VectorUblas<double>>> vublas(m,"VectorUBlas<double,ublas::vector<double>>");
    vublas.def(py::init<>())
        .def(py::self + py::self )
        .def(py::self - py::self)
        .def(double() + py::self)
        .def(double() - py::self)
        .def(py::self + double())
        .def(py::self - double())
        .def(py::self * double())
        .def(double() * py::self)
        //.def(double() / py::self)
        .def(py::self += py::self)
        .def(py::self -= py::self)
        .def(py::self *= double())
        .def(-py::self)
//        .def(py::self /= double())
//        .def(py::self *= py::self)
//        .def(py::self /= py::self)
        ;
}
