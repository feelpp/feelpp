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
//! @date 25 Jun 2021
//! @copyright 2018-2021 Feel++ Consortium
//!
#include <feel/feelcore/pybind11_json.hpp>
#include <pybind11/pybind11.h>

#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>

namespace py = pybind11;
using namespace Feel;

PYBIND11_MODULE( _modelalg, m )
{
    using namespace Feel;
    using namespace Feel::FeelModels;

    py::class_<ModelAlgebraicFactory, std::shared_ptr<ModelAlgebraicFactory>>( m, "ModelAlgebraicFactory" )
        //.def( py::init<std::string, backend_ptr_t<double>>(),"Initialize ModelAlgebraicFactory" )
        .def( "matrix", &ModelAlgebraicFactory::matrix, "get the matrix" )
        .def( "rhs", &ModelAlgebraicFactory::rhs, "get the right hand side" )
        .def( "backend", &ModelAlgebraicFactory::backend, "get the backend" );
}
