/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2021-06-02

  Copyright (C) 2021 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef FEELPP_DISCR_DETAIL_SPACECONTEXT_HPP
#define FEELPP_DISCR_DETAIL_SPACECONTEXT_HPP

#include <feel/feelcore/feel.hpp>

namespace Feel
{
namespace detail
{

template <typename MeshType, typename ... CTX>
auto selectGeomapContextFromSpaceContext( std::shared_ptr<MeshType> const& mesh, CTX const& ... ctx )
{
    using geoelement_type = typename MeshType::element_type;
    using gmc_type = typename geoelement_type::gm_type::template Context<geoelement_type>;
    std::shared_ptr<gmc_type> res;
    hana::for_each( hana::make_tuple( ctx... ), [&mesh,&res]( auto const& e )
                    {
                        if constexpr ( std::is_same_v<std::decay_t<decltype(e->gmContext())>, std::shared_ptr<gmc_type> > )
                        {
                            if ( e->gmContext()->element().mesh()->isSameMesh( mesh ) )
                                res = e->gmContext();
                        }
                    } );
    return res;
}


} // namespace detail
} // namespace Feel

#endif
