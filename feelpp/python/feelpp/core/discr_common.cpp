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
#include <feel/feelpython/pybind11/stl_bind.h>
#include "discr_bindings.hpp"

void
bindDiscrCommon( py::module& m )
{
    std::string pyclass_name = std::string("ComponentType");
    py::enum_<ComponentType>(m,pyclass_name.c_str())
        .value("NO_COMPONENT", ComponentType::NO_COMPONENT )
        .value("X", ComponentType::X )
        .value("Y", ComponentType::Y )
        .value("Z", ComponentType::Z )
        .value("NX", ComponentType::NX )
        .value("NY", ComponentType::NY )
        .value("NZ", ComponentType::NZ )
        .value("TX", ComponentType::TX )
        .value("TY", ComponentType::TY )
        .value("TZ", ComponentType::TZ )
        .export_values();

    pyclass_name = std::string("DofTableExtendedType");
    py::enum_<DofTableExtendedType>(m,pyclass_name.c_str())
        .value("NONE", DofTableExtendedType::NONE )
        .value("VERTICES", DofTableExtendedType::VERTICES )
        .value("DEFAULT", DofTableExtendedType::DEFAULT )
        .export_values();

    py::bind_vector<std::vector<DofTableExtendedType>>(m, "VectorDofTableExtendedType");

    pyclass_name = std::string("Periodic");
    py::class_<Periodic<double>>(m,pyclass_name.c_str()).def(py::init<>());
    pyclass_name = std::string("PeriodicityPeriodic");
    py::class_<Periodicity<Periodic<double>>>(m,pyclass_name.c_str()).def(py::init<>());
    pyclass_name = std::string("NoPeriodicity");
    py::class_<NoPeriodicity>(m,pyclass_name.c_str()).def(py::init<>());
    pyclass_name = std::string("PeriodicityNoPeriodicity");
    py::class_<Periodicity<NoPeriodicity>>(m,pyclass_name.c_str()).def(py::init<>());
}
