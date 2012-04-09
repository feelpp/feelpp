/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-01-21

  Copyright (C) 2009-2011 Universite Joseph Fourier (Grenoble I)

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
   \file bases.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-01-21
 */
#ifndef __Bases_H
#define __Bases_H 1

#include <boost/mpl/at.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/support/is_sequence.hpp>
#include <boost/fusion/sequence.hpp>

namespace Feel
{
namespace mpl = boost::mpl;
using boost::mpl::_;

namespace detail
{
struct bases_base {};
struct meshes_base {};


/**
 * \class bases
 * \brief classes that store sequences of basis functions to define function spaces
 *
 * @author Christophe Prud'homme
 * @see
 */



template <class A0=mpl::void_, class A1=mpl::void_, class A2=mpl::void_, class A3=mpl::void_>
struct bases
        :
    public detail::bases_base,
    public mpl::if_<boost::is_same<A1,mpl::void_>,
        boost::fusion::vector<A0>,
        typename mpl::if_<boost::is_same<A2,mpl::void_>,
        boost::fusion::vector<A0,A1>,
        typename mpl::if_<boost::is_same<A3,mpl::void_>,
        boost::fusion::vector<A0,A1,A2>,
        boost::fusion::vector<A0,A1,A2,A3> >::type>::type>::type
{
};

} // namespace detail


struct void_basis : public mpl::void_
{
    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef mpl::void_ type;
    };
};

template <class A0=void_basis, class A1=void_basis, class A2=void_basis, class A3=void_basis>
struct bases
        :
    public detail::bases_base,
    public mpl::if_<boost::is_same<A1,void_basis>,
        boost::fusion::vector<typename A0::template ChangeTag<0>::type >,
      typename mpl::if_<boost::is_same<A2,void_basis>,
      boost::fusion::vector<typename A0::template ChangeTag<0>::type,
            typename A1::template ChangeTag<1>::type >,
                     typename mpl::if_<boost::is_same<A3,void_basis>,
                     boost::fusion::vector<typename A0::template ChangeTag<0>::type,
                           typename A1::template ChangeTag<1>::type,
                                    typename A2::template ChangeTag<2>::type >,
                                             boost::fusion::vector<typename A0::template ChangeTag<0>::type,
                                                   typename A1::template ChangeTag<1>::type,
                                                            typename A2::template ChangeTag<2>::type,
                                                                     typename A3::template ChangeTag<3>::type > >::type>::type>::type
{
};


template <class A0=mpl::void_, class A1=mpl::void_, class A2=mpl::void_, class A3=mpl::void_>
struct meshes
        :
    public detail::meshes_base,
    public mpl::if_<boost::is_same<A1,mpl::void_>,
        boost::fusion::vector<A0>,
        typename mpl::if_<boost::is_same<A2,mpl::void_>,
        boost::fusion::vector<A0,A1>,
        typename mpl::if_<boost::is_same<A3,mpl::void_>,
        boost::fusion::vector<A0,A1,A2>,
        boost::fusion::vector<A0,A1,A2,A3> >::type>::type>::type
{
};


}// Feel

#endif /* __Bases_H */
