/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-12-04

  Copyright (C) 2013 Université de Strasbourg

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
/**
   \file isp0continuous.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-12-04
*/
#if !defined(FEELPP_IS_P0_CONTINUOUS_HPP)
#define FEELPP_IS_P0_CONTINUOUS_HPP 1

namespace Feel
{
template<typename PSet>
struct isP0Continuous
{
    static constexpr bool result = mpl::and_< boost::is_same<mpl::int_<PSet::nOrder>,mpl::int_<0> >,
                                              boost::is_same<typename PSet::continuity_type, Continuous > >::type::value;
};
}

#endif /* FEELPP_IS_P0_CONTINUOUS_HPP */
