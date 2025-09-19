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
//! @date 23 Sept 2021
//! @copyright 2021 Feel++ Consortium
//!
#include <feel/feelpython/pybind11/pybind11.h>

#include <fmt/core.h>
#include <mpi4py/mpi4py.h>

#include <feel/feeltiming/tic.hpp>

namespace py = pybind11;


PYBIND11_MODULE(_timing, m )
{
    using namespace Feel;
    m.def(
        "tic", &tic,
        fmt::format( "Record internal time at its execution. To be used with toc()." ).c_str() );
    m.def(
        "toc", &toc,
        py::arg( "msg" ) = std::string( "" ), py::arg( "display" ) = display, py::arg( "uiname" ) = std::string( "" ),
        fmt::format( "toc returns the time elapsed since the last tic" ).c_str() );
}
