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
//! @file meshmover.cpp
//! @author Luca Berti 
//! @date 03 Mar 2023
//! @copyright 2023 Feel++ Consortium
//!
#include <boost/shared_ptr.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelpython/pybind11/pybind11.h>
#include <mpi4py/mpi4py.h>

namespace py = pybind11;


template<typename MeshT, typename DisplType>
void defMeshmover(py::module &m)
{
    using mesh_t = MeshT;
    using displ_t = DisplType;
    using mesh_ptr_t = std::shared_ptr<MeshT>;

    std::string pyclass_name = fmt::format("MeshMover{}DG{}", MeshT::nDim, MeshT::nOrder );

    m.def("meshMove",[]( mesh_ptr_t & mesh , displ_t const& u ) 
        {return meshMove( mesh,u);}, 
        py::arg("mesh"), py::arg("disp"),"move mesh according to displacement u" );
           
}

PYBIND11_MODULE(_meshmover, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;

    auto ordert = hana::make_tuple( 1_c,2_c );
    hana::for_each( ordert, [&m](auto const& o ){
        constexpr int _order = std::decay_t<decltype(o)>::value;
        defMeshmover<Mesh<Simplex<2>>,typename Pch_type< Mesh<Simplex<2>>,_order>::element_type >( m );
        defMeshmover<Mesh<Simplex<2>>,typename Pchv_type< Mesh<Simplex<2>>,_order>::element_type >( m );
        defMeshmover<Mesh<Simplex<3>>,typename Pch_type< Mesh<Simplex<3>>,_order>::element_type >( m );
        defMeshmover<Mesh<Simplex<3>>,typename Pchv_type< Mesh<Simplex<3>>,_order>::element_type >( m );
    });
}