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

#include<feel/feelcore/environment.hpp>
#include <feel/feelts/tsbase.hpp>
#include <mpi4py/mpi4py.h>

#include <boost/shared_ptr.hpp>

namespace py = pybind11;

using namespace Feel;


PYBIND11_MODULE(_ts, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;

    py::class_<TSBase, std::shared_ptr<TSBase>>(m,"TSBase")
        .def(py::init<>())
        .def("isFinished", &TSBase::isFinished, "return if time stepping is finished, false otherwise")
        ;
}
