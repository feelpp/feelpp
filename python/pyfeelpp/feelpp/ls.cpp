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

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feells/distancetorange.hpp>
#include <feel/feelmesh/metric.hpp>
#include <feel/feelvf/vf.hpp>
#include <pybind11/pybind11.h>

#include <fmt/core.h>
#include <mpi4py/mpi4py.h>
namespace py = pybind11;

template<int Dim>
void
dist2range_inst( py::module &m )
{
    using namespace Feel;using namespace Feel::vf;
    using mesh_t = Mesh < Simplex<Dim, 1>>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    m.def(
        "distanceToRange", []( Pch_ptrtype<mesh_t, 1> const& Xh, faces_reference_wrapper_t<mesh_ptr_t> const& facets, double maxDistance, double fmStride )
        { return distanceToRange( _space=Xh, _range=facets, _max_distance=maxDistance, _fm_stride=fmStride ); },
        py::return_value_policy::copy,
        py::arg( "space" ),
        py::arg( "faces" ),
        py::arg( "maxDistance" ) = -1.,
        py::arg( "fmStride" ) = -1.,
        fmt::format("compute the distance field in space to the range of faces in {}D",Dim).c_str() );
    m.def(
        "gradedls", []( Pch_ptrtype<mesh_t, 1> const& Xh, faces_reference_wrapper_t<mesh_ptr_t> const& facets, double hclose, double hfar )
        { 
            return gradedfromls( Xh, facets, hclose, hfar);
        },
        py::return_value_policy::copy,
        py::arg( "space" ),
        py::arg( "faces" ),
        py::arg( "hclose" ),
        py::arg( "hfar" ),
        fmt::format("compute the graded metric field with metric set to hclose on facets and hfar the farthest away from facets in {}D",Dim).c_str() );
    m.def(
        "expr", []( Pch_ptrtype<mesh_t, 1> const& Xh, std::string const& e )
        { 
            return expr( Xh, expr(e) );
        },
        py::return_value_policy::copy,
        py::arg( "space" ),
        py::arg( "expr" ),
        fmt::format("compute the  metric field with metric set from expression in {}D",Dim).c_str() );
}
PYBIND11_MODULE(_ls, m )
{
    using namespace Feel;
    using namespace Feel::vf;

    if (import_mpi4py()<0) return ;
    dist2range_inst<2>(m);
    dist2range_inst<3>(m);
                                                      
}
