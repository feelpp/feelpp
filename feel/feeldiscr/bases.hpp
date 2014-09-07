/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-21
 */
#ifndef __Bases_H
#define __Bases_H 1

#include <boost/mpl/at.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/support/is_sequence.hpp>
#include <boost/fusion/sequence.hpp>
#include <feel/feeldiscr/mortar.hpp>

namespace Feel
{
namespace mpl = boost::mpl;
using boost::mpl::_;
    /// @cond DETAIL
namespace detail
{
struct bases_base {};
struct meshes_base {};
struct mortars_base {};
struct periodic_base {};
/**
 *
 * \brief classes that store sequences of basis functions to define function spaces
 *
 * @author Christophe Prud'homme
 * @see
 */
#if FEELPP_CLANG_AT_LEAST(3,1) || FEELPP_GNUC_AT_LEAST(4,7)

template<typename... Args>
struct bases
    :
        public Feel::detail::bases_base,
        public boost::fusion::vector<Args...>
{};

#else

template <class A0=mpl::void_, class A1=mpl::void_, class A2=mpl::void_, class A3=mpl::void_, class A4=mpl::void_>
struct bases
        :
        public Feel::detail::bases_base,
        public mpl::if_<boost::is_same<A1,mpl::void_>,
                        boost::fusion::vector<A0>,
                        typename mpl::if_<boost::is_same<A2,mpl::void_>,
                                          boost::fusion::vector<A0,A1>,
                                          typename mpl::if_<boost::is_same<A3,mpl::void_>,
                                                            boost::fusion::vector<A0,A1,A2>,
                                                            typename mpl::if_<boost::is_same<A4,mpl::void_>,
                                                                              boost::fusion::vector<A0,A1,A2,A3>,
                                                                              boost::fusion::vector<A0,A1,A2,A3,A4> >::type>::type>::type>::type
{
};

#endif

template <class BasisFusionVectorType>
struct bases2
    :
        public Feel::detail::bases_base,
        public BasisFusionVectorType
{
};
} // namespace detail
/// @endcond
#if FEELPP_CLANG_AT_LEAST(3,1) || FEELPP_GNUC_AT_LEAST(4,7)

struct ChangeBasisTag
{
public:
    template<typename Sig>
    struct result;

    template<typename Lhs, typename Rhs>
    struct result<ChangeBasisTag( Lhs,Rhs )>
    {
	    typedef typename boost::remove_const<typename boost::remove_reference<Lhs>::type>::type lhs_noref_type;
	    typedef typename boost::remove_const<typename boost::remove_reference<Rhs>::type>::type rhs_noref_type;

	    typedef typename fusion::result_of::size<lhs_noref_type>::type index;
	    typedef typename fusion::result_of::push_back<lhs_noref_type, typename rhs_noref_type::template ChangeTag<index::value>::type>::type type;
    };

};
template<typename... Args>
struct bases
    :
        public Feel::detail::bases_base,
        public fusion::result_of::as_vector<typename fusion::result_of::accumulate<fusion::vector<Args...>, fusion::vector<>, ChangeBasisTag >::type>::type
{};


template<typename... Args>
struct meshes
    :
        public Feel::detail::meshes_base,
        public boost::fusion::vector<Args...>
{
    typedef boost::fusion::vector<Args...> super;
    typedef meshes<Args...> this_type;
	static const int s = sizeof...(Args);
    meshes( super const& m) : super( m ) {}
};

template<typename... Args>
struct mortars
    :
        public Feel::detail::mortars_base,
        public Feel::detail::mortar_base,
        public boost::fusion::vector<Args...>
{
    typedef boost::fusion::vector<Args...> super;
    typedef mortars<Args...> this_type;
	static const int s = sizeof...(Args);
    mortars( super const& m) : super( m ) {}
};

template<typename... Args>
struct Periodicity
    :
        public Feel::detail::periodic_base,
        public Feel::detail::periodicity_base,
        public boost::fusion::vector<Args...>
{
    typedef boost::fusion::vector<Args...> super;
    typedef Periodicity<Args...> this_type;
	static const int s = sizeof...(Args);
    Periodicity() : super() {}
    Periodicity( super const& m) : super( m ) {}
    Periodicity( Args... args ) : super( fusion::make_vector(args...) ) {}

    uint16_type tag1() const { return fusion::at_c<0>(*this).tag1(); }
    uint16_type tag2() const { return fusion::at_c<0>(*this).tag2(); }
};

template<typename... Args>
Periodicity<Args...> periodicity( Args... args )
{ return Periodicity<Args...>( args... ); }

#else



struct void_basis : public mpl::void_
{
    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef mpl::void_ type;
    };
};

template <class A0=void_basis, class A1=void_basis, class A2=void_basis, class A3=void_basis, class A4=void_basis>
struct bases
        :
        public Feel::detail::bases_base,
        public mpl::if_<boost::is_same<A1,void_basis>,
                        boost::fusion::vector<typename A0::template ChangeTag<0>::type >,
                        typename mpl::if_<boost::is_same<A2,void_basis>,
                                          boost::fusion::vector<typename A0::template ChangeTag<0>::type,
                                                                typename A1::template ChangeTag<1>::type >,
                                          typename mpl::if_<boost::is_same<A3,void_basis>,
                                                            boost::fusion::vector<typename A0::template ChangeTag<0>::type,
                                                                                  typename A1::template ChangeTag<1>::type,
                                                                                  typename A2::template ChangeTag<2>::type >,
                                                            typename mpl::if_<boost::is_same<A4,void_basis>,
                                                                              boost::fusion::vector<typename A0::template ChangeTag<0>::type,
                                                                                                    typename A1::template ChangeTag<1>::type,
                                                                                                    typename A2::template ChangeTag<2>::type,
                                                                                                    typename A3::template ChangeTag<3>::type >,
                                                                              boost::fusion::vector<typename A0::template ChangeTag<0>::type,
                                                                                                    typename A1::template ChangeTag<1>::type,
                                                                                                    typename A2::template ChangeTag<2>::type,
                                                                                                    typename A3::template ChangeTag<3>::type,
                                                                                                    typename A4::template ChangeTag<4>::type > >::type>::type>::type>::type
{
};

template <class A0=mpl::void_, class A1=mpl::void_, class A2=mpl::void_, class A3=mpl::void_, class A4=mpl::void_>
struct meshes
        :
        public Feel::detail::meshes_base,
        public mpl::if_<boost::is_same<A0,mpl::void_>,
                        boost::fusion::vector<>,
                        typename mpl::if_<boost::is_same<A1,mpl::void_>,
                                          boost::fusion::vector<A0>,
                                          typename mpl::if_<boost::is_same<A2,mpl::void_>,
                                                            boost::fusion::vector<A0,A1>,
                                                            typename mpl::if_<boost::is_same<A3,mpl::void_>,
                                                                              boost::fusion::vector<A0,A1,A2>,
                                                                              typename mpl::if_<boost::is_same<A4,mpl::void_>,
                                                                                                boost::fusion::vector<A0,A1,A2,A3>,
                                                                                                boost::fusion::vector<A0,A1,A2,A3,A4> >::type>::type>::type>::type>::type

{
    typedef typename mpl::if_<boost::is_same<A0,mpl::void_>,
                              boost::fusion::vector<>,
                              typename mpl::if_<boost::is_same<A1,mpl::void_>,
                                                boost::fusion::vector<A0>,
                                                typename mpl::if_<boost::is_same<A2,mpl::void_>,
                                                                  boost::fusion::vector<A0,A1>,
                                                                  typename mpl::if_<boost::is_same<A3,mpl::void_>,
                                                                                    boost::fusion::vector<A0,A1,A2>,
                                                                                    typename mpl::if_<boost::is_same<A4,mpl::void_>,
                                                                                                      boost::fusion::vector<A0,A1,A2,A3>,
                                                                                                      boost::fusion::vector<A0,A1,A2,A3,A4> >::type>::type>::type>::type>::type super;

    typedef meshes<A0,A1,A2,A3,A4> this_type;
    meshes( super const& m ) : super( m ) {}
};

template <class A0=mpl::void_, class A1=mpl::void_, class A2=mpl::void_, class A3=mpl::void_, class A4=mpl::void_>
struct mortars
    :
        public Feel::detail::mortars_base,
        public mpl::if_<boost::is_same<A0,mpl::void_>,
                        boost::fusion::vector<>,
                        typename mpl::if_<boost::is_same<A1,mpl::void_>,
                                          boost::fusion::vector<A0>,
                                          typename mpl::if_<boost::is_same<A2,mpl::void_>,
                                                            boost::fusion::vector<A0,A1>,
                                                            typename mpl::if_<boost::is_same<A3,mpl::void_>,
                                                                              boost::fusion::vector<A0,A1,A2>,
                                                                              typename mpl::if_<boost::is_same<A4,mpl::void_>,
                                                                                                boost::fusion::vector<A0,A1,A2,A3>,
                                                                                                boost::fusion::vector<A0,A1,A2,A3,A4> >::type>::type>::type>::type>::type

{
    typedef typename mpl::if_<boost::is_same<A0,mpl::void_>,
                              boost::fusion::vector<>,
                              typename mpl::if_<boost::is_same<A1,mpl::void_>,
                                                boost::fusion::vector<A0>,
                                                typename mpl::if_<boost::is_same<A2,mpl::void_>,
                                                                  boost::fusion::vector<A0,A1>,
                                                                  typename mpl::if_<boost::is_same<A3,mpl::void_>,
                                                                                    boost::fusion::vector<A0,A1,A2>,
                                                                                    typename mpl::if_<boost::is_same<A4,mpl::void_>,
                                                                                                      boost::fusion::vector<A0,A1,A2,A3>,
                                                                                                      boost::fusion::vector<A0,A1,A2,A3,A4> >::type>::type>::type>::type>::type super;

    typedef mortars<A0,A1,A2,A3,A4> this_type;
    mortars( super const& m ) : super( m ) {}
};

template <class A0=mpl::void_, class A1=mpl::void_, class A2=mpl::void_, class A3=mpl::void_, class A4=mpl::void_>
struct Periodicity
    :
        public Feel::detail::periodic_base,
        public Feel::detail::periodicity_base,
        public mpl::if_<boost::is_same<A0,mpl::void_>,
                        boost::fusion::vector<>,
                        typename mpl::if_<boost::is_same<A1,mpl::void_>,
                                          boost::fusion::vector<A0>,
                                          typename mpl::if_<boost::is_same<A2,mpl::void_>,
                                                            boost::fusion::vector<A0,A1>,
                                                            typename mpl::if_<boost::is_same<A3,mpl::void_>,
                                                                              boost::fusion::vector<A0,A1,A2>,
                                                                              typename mpl::if_<boost::is_same<A4,mpl::void_>,
                                                                                                boost::fusion::vector<A0,A1,A2,A3>,
                                                                                                boost::fusion::vector<A0,A1,A2,A3,A4> >::type>::type>::type>::type>::type

{
    typedef typename mpl::if_<boost::is_same<A0,mpl::void_>,
                              boost::fusion::vector<>,
                              typename mpl::if_<boost::is_same<A1,mpl::void_>,
                                                boost::fusion::vector<A0>,
                                                typename mpl::if_<boost::is_same<A2,mpl::void_>,
                                                                  boost::fusion::vector<A0,A1>,
                                                                  typename mpl::if_<boost::is_same<A3,mpl::void_>,
                                                                                    boost::fusion::vector<A0,A1,A2>,
                                                                                    typename mpl::if_<boost::is_same<A4,mpl::void_>,
                                                                                                      boost::fusion::vector<A0,A1,A2,A3>,
                                                                                                      boost::fusion::vector<A0,A1,A2,A3,A4> >::type>::type>::type>::type>::type super;

    typedef Periodicity<A0,A1,A2,A3,A4> this_type;
    Periodicity() : super() {}
    Periodicity( super const& m ) : super( m ) {}
    Periodicity( A0 a0 ) : super( fusion::make_vector(a0) ) {}
    Periodicity( A0 a0, A1 a1 ) : super( fusion::make_vector(a0,a1) ) {}
    Periodicity( A0 a0, A1 a1, A2 a2 ) : super( fusion::make_vector(a0,a1,a2) ) {}
    Periodicity( A0 a0, A1 a1, A2 a2, A3 a3 ) : super( fusion::make_vector(a0,a1,a2,a3) ) {}

    uint16_type tag1() const { return fusion::at_c<0>(*this).tag1(); }
    uint16_type tag2() const { return fusion::at_c<0>(*this).tag2(); }
};

template<typename A0> Periodicity<A0> periodicity( A0 const& a0 ) { return Periodicity<A0>( a0 ); }
template<typename A0,typename A1> Periodicity<A0,A1> periodicity( A0 const& a0, A1 const& a1 ) { return Periodicity<A0,A1>( a0,a1 ); }
template<typename A0,typename A1,typename A2> Periodicity<A0,A1,A2> periodicity( A0 const& a0, A1 const& a1, A2 const& a2 ) { return Periodicity<A0,A1,A2>( a0,a1,a2 ); }

#endif

}// Feel

#endif /* __Bases_H */
