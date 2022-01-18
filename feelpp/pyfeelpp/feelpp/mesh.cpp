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
#include "mesh.hpp"

void mesh1( py::module& m );
void mesh2( py::module& m );
void mesh3( py::module& m );

PYBIND11_MODULE( _mesh, m )
{
    using namespace Feel;

    if ( import_mpi4py() < 0 ) return;

    py::enum_<GeomapStrategyType>( m, "GeomapStrategyType" )
        .value( "GEOMAP_OPT", GeomapStrategyType::GEOMAP_OPT )
        .value( "GEOMAP_O1", GeomapStrategyType::GEOMAP_O1 )
        .value( "GEOMAP_HO", GeomapStrategyType::GEOMAP_HO )
        .export_values();

    mesh1( m );
    mesh2( m );
    mesh3( m );

    defMesh<MeshStructured<Hypercube<2>>>(m);
}
