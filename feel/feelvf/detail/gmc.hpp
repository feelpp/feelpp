/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013 Feel++ Consortium

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
   \file gmc.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_DETAIL_GMC_HPP)
#define FEELPP_DETAIL_GMC_HPP 1

#include <boost/fusion/sequence.hpp>

namespace Feel {

//namespace blas = boost::numeric::bindings::blas;
//namespace traits = boost::numeric::bindings::traits;
namespace fusion = boost::fusion;


namespace vf { namespace detail {

/// \cond detail
template<int Index> struct gmc
{
    static const int value = Index;
    typedef mpl::void_ reference_element_type;
} ;

template<typename Geo_t>
struct ExtractGm
{
    typedef typename mpl::if_<fusion::result_of::has_key<Geo_t,vf::detail::gmc<0> >,mpl::identity<vf::detail::gmc<0> >,mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
    typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
    typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

    static gmc_ptrtype get( Geo_t const& geom )
    {
        return fusion::at_key<key_type>( geom ).get();
    }
    static Geo_t clone( Geo_t const& geom )
    {
        Geo_t geom2( geom );
        fusion::at_key<key_type>( geom2 )  = fusion::at_key<key_type>( geom )->clone();
        return geom2;
    }
};
/// \endcond
} } }
#endif /* FEELPP_DETAIL_GMC_HPP */
