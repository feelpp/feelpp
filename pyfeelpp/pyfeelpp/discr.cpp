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
#include <pybind11/pybind11.h>

#include <feel/feeldiscr/pch.hpp>
#include <mpi4py/mpi4py.h>

namespace py = pybind11;
using namespace Feel;

template<typename MeshT, int Order = 1>
void defDiscr(py::module &m)
{
    using namespace Feel;
    std::string pyclass_name = std::string("Pch_") + std::to_string(MeshT::nDim)+std::string("D_P") + std::to_string(Order);
    using pch_t = Pch_type<MeshT,Order>;
    using mesh_support_vector_t = typename pch_t::mesh_support_vector_type;
    using periodicity_t = typename pch_t::periodicity_type;
    py::class_<pch_t,std::shared_ptr<pch_t>>(m,pyclass_name.c_str())
        .def(py::init<std::shared_ptr<MeshT> const&,mesh_support_vector_t const&, size_type, periodicity_t, std::vector<WorldComm> const&, std::vector<bool>>(),
             py::arg("mesh"),
             py::arg("support")=mesh_support_vector_t(),
             py::arg("components")=MESH_RENUMBER | MESH_CHECK,
             py::arg("periodicity")=periodicity_t(),
             py::arg("worldsComm") = Environment::worldsComm(1),
             py::arg("extendedDofTable") = std::vector<bool>(1,false) );

    
}
    
PYBIND11_MODULE(discr, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;
    
    defDiscr<Mesh<Simplex<1>>,1>( m );
    defDiscr<Mesh<Simplex<2>>,1>( m );
    defDiscr<Mesh<Simplex<2>>,2>( m );
}

