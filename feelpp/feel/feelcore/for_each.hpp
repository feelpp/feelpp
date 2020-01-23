/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Thibaut Metivet <thibaut.metivet@inria.fr>
       Date: 2020-01-23

  Copyright (C) 2020 Inria

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file for_each.hpp
   \author Thibaut Metivet <thibaut.metivet@inria.fr>
   \date 2020-01-23
 */

#ifndef FEELPP_FOREACH_HPP
#define FEELPP_FOREACH_HPP 1

#include <feel/feelcore/traits.hpp>
#include <boost/hana/for_each.hpp>
#include <algorithm>

namespace Feel {

template<typename Xs, typename F,
    std::enable_if_t<boost::hana::Foldable<Xs>::value, int> = 0 >
void for_each( Xs && xs, F && f )
{
    boost::hana::for_each( xs, f );
}

template<typename Xs, typename F,
    std::enable_if_t<Feel::is_iterable<Xs>::value, int> = 0 >
void for_each( Xs && xs, F && f )
{
    std::for_each( begin(xs), end(xs), f );
}

} // namespace Feel
#endif
