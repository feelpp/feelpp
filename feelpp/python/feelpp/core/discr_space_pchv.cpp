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
#include "discr_bindings.hpp"

namespace
{
template<int Dim, int Geo, int Order>
void bindPchvOrder( py::module& m )
{
    defDiscr<Pchv_type<Mesh<Simplex<Dim, Geo>>, Order>, Order>( m );
}

template<int Dim, int Geo>
void bindPchvGeo( py::module& m )
{
    bindPchvOrder<Dim, Geo, 0>( m );
    bindPchvOrder<Dim, Geo, 1>( m );
    bindPchvOrder<Dim, Geo, 2>( m );
    bindPchvOrder<Dim, Geo, 3>( m );
    bindPchvOrder<Dim, Geo, Dynamic>( m );
}

template<int Dim>
void bindPchvDim( py::module& m )
{
    bindPchvGeo<Dim, 1>( m );
    bindPchvGeo<Dim, 2>( m );
}
}

void
bindDiscrPchv( py::module& m )
{
    bindPchvDim<1>( m );
    bindPchvDim<2>( m );
    bindPchvDim<3>( m );
}
