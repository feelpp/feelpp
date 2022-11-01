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

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/quality.hpp>

namespace py = pybind11;

template<int Dim>
void
quality_inst( py::module &m )
{
    using namespace Feel;using namespace Feel::vf;
    using mesh_t = Mesh < Simplex<Dim, 1>>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    m.def(
        "etaQ", []( mesh_ptr_t const& mesh )
        { return etaQ( mesh ); },
        py::return_value_policy::copy,
        py::arg( "mesh" ),
        fmt::format("compute the piecewise constant field containing the element quality according to etaQ indicator {}D",Dim).c_str() );
    m.def(
        "nsrQ", []( mesh_ptr_t const& mesh )
        { return nsrQ( mesh ); },
        py::return_value_policy::copy,
        py::arg( "mesh" ),
        fmt::format("compute the piecewise constant field containing the element quality according to etaQ indicator {}D",Dim).c_str() );
}
PYBIND11_MODULE(_minmax, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;
    quality_inst<2>(m);
    quality_inst<3>(m);
                                                      
}
