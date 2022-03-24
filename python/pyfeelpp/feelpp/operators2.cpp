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
//! @date 23 Mar 2022
//! @copyright 2022 Feel++ Consortium
//!
#include "operators.hpp"

void operators2( py::module& m )
{
    using namespace Feel;
    using namespace hana::literals;
    if ( import_mpi4py() < 0 ) return;

    auto dimt = hana::make_tuple( hana::int_c<2> );
    auto ordert = hana::make_tuple( hana::int_c<0>, hana::int_c<1>, hana::int_c<2>, hana::int_c<3> );
    auto geot = hana::make_tuple( hana::int_c<1>, hana::int_c<2> );
    hana::for_each(hana::cartesian_product(hana::make_tuple(dimt, ordert, geot)),
                   [&m]( auto const& d )
                       {
                           constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                           constexpr int _order = std::decay_t<decltype(hana::at_c<1>(d))>::value;
                           constexpr int _geo = std::decay_t<decltype( hana::at_c<2>( d ) )>::value;
                           defOperator<Pch_type<Mesh<Simplex<_dim, _geo>>, _order>>( m );
                           defOperator<Pdh_type<Mesh<Simplex<_dim, _geo>>, _order>>( m );
                           defOperator<Pchv_type<Mesh<Simplex<_dim, _geo>>, _order>>( m );
                           defOperator<Pdhv_type<Mesh<Simplex<_dim, _geo>>, _order>>( m );
                       });
}
