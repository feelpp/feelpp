/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-03-28

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file test_fusion.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-03-28
*/
#include <cassert>

#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/container.hpp>
#include <boost/fusion/container.hpp>

template<typename T>
class F2
{
public:
    template<typename Sig>
    struct result;

    template<typename Lhs, typename Rhs>
    struct result<F2( Lhs,Rhs )>
    {
#if 1
        typedef typename boost::remove_const<typename boost::remove_reference<Lhs>::type>::type lhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<Rhs>::type>::type rhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<T>::type>::type T_noref_type;

        typedef typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::push_back<lhs_noref_type,
        typename boost::fusion::result_of::make_pair<T_noref_type,rhs_noref_type>::type>::type>::type type;
#else
        typedef typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::push_back<Lhs,
        typename boost::fusion::result_of::make_pair<T,Rhs>::type>::type>::type type;
        //typedef typename boost::fusion::vector< typename boost::fusion::result_of::make_pair<T_noref_type,rhs_noref_type>::type> type;
        //typedef typename boost::fusion::result_of::as_vector< lhs_noref_type,typename boost::fusion::result_of::make_pair<T_noref_type,rhs_noref_type>::type> type;

#endif
    };
    T m_T;
    F2( T const& t ) : m_T( t ) {}

    template<typename Lhs, typename Rhs>
    typename result<F2( Lhs,Rhs )>::type
    operator()( Lhs const&  lhs, Rhs const& rhs ) const
    {
#if 1
        typedef typename boost::remove_const<typename boost::remove_reference<Lhs>::type>::type lhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<Rhs>::type>::type rhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<T>::type>::type T_noref_type;
        //return boost::fusion::push_back( lhs, boost::fusion::make_pair<T>( rhs ) ) ;
        //return boost::fusion::make_vector(boost::fusion::make_pair<T>( rhs ));
        return boost::fusion::as_vector( boost::fusion::push_back( lhs,boost::fusion::make_pair<T>( rhs ) ) );
        //return boost::fusion::push_back( lhs, boost::fusion::make_pair<T_noref_type>( rhs ) ) ;
#else
        return boost::fusion::as_vector( boost::fusion::push_back( lhs,boost::fusion::make_pair<T>( rhs ) ) );
#endif
    }

};
template<typename T>
class F1
{
public:
    template<typename Sig>
    struct result;

    template<typename Lhs, typename Rhs>
    struct result<F1( Lhs,Rhs )>
    {
        typedef typename boost::remove_const<typename boost::remove_reference<Lhs>::type>::type lhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<Rhs>::type>::type rhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<T>::type>::type T_noref_type;
        typedef typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::join<lhs_noref_type,
        typename boost::fusion::result_of::accumulate<T_noref_type,
        const boost::fusion::vector<>,
        F2<rhs_noref_type> >::type>::type>::type type;
    };

    T m_T;
    F1( T const& t ) : m_T( t ) {}

    template<typename Lhs, typename Rhs>
    typename result<F1( Lhs,Rhs )>::type
    operator()( Lhs const& lhs, Rhs const&   rhs ) const
    {
        typedef typename boost::remove_const<typename boost::remove_reference<Lhs>::type>::type lhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<Rhs>::type>::type rhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<T>::type>::type T_noref_type;
        return boost::fusion::as_vector( boost::fusion::join( lhs, boost::fusion::accumulate( m_T, boost::fusion::vector<>(), F2<Rhs>( rhs ) ) ) );
        ;
    }
};

template<typename T>
class F0
{
public:
    template<typename Sig>
    struct result;

    template<typename Lhs, typename Rhs>
    struct result<F0( Lhs,Rhs )>
    {

        typedef typename boost::remove_reference<Lhs>::type lhs_noref_type;
        typedef typename boost::remove_reference<Rhs>::type rhs_noref_type;
        typedef typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::push_back<lhs_noref_type,rhs_noref_type>::type>::type type;

    };

    T m_T;
    F0( T const& t ) : m_T( t ) {}

    template<typename Lhs, typename Rhs>
    typename result<F0( Lhs,Rhs )>::type
    operator()( Lhs  lhs, Rhs    rhs )  const
    {
        return boost::fusion::as_vector( boost::fusion::push_back( lhs, rhs ) );
    }
};

int main()
{

    //boost::fusion::vector<boost::mpl::int_<1>,boost::mpl::int_<2> > v1;
    //boost::fusion::vector<boost::mpl::int_<5>,boost::mpl::int_<6>, boost::mpl::int_<7> > v2;
    boost::fusion::vector2<boost::mpl::int_<1>,boost::mpl::int_<2> > v1;
    boost::fusion::vector<boost::mpl::int_<5>, boost::mpl::int_<6> > v2;
    auto v3 = boost::fusion::accumulate( v1, boost::fusion::vector<>(), F0<decltype( v2 )>( v2 ) );
    BOOST_MPL_ASSERT( ( boost::is_same<decltype( v3 ),decltype( v1 )> ) );
    assert( ( boost::is_same<decltype( v3 ),decltype( v1 )>::value ) );

    //const boost::fusion::vector<int,int> vec(1,2);
    //assert( boost::is_same<boost::fusion::insert(vec,  boost::fusion::next(begin(vec)), 3),boost::fusion::make_vector(1,3,2));

    boost::fusion::vector<> v;
    auto v4 = boost::fusion::accumulate( v1, v, F1<decltype( v2 )>( v2 ) );
    typedef boost::fusion::vector4<boost::fusion::pair<boost::mpl::int_<1>,boost::mpl::int_<5> >, boost::fusion::pair<boost::mpl::int_<1>,boost::mpl::int_<6> >, boost::fusion::pair<boost::mpl::int_<2>,boost::mpl::int_<5> >, boost::fusion::pair<boost::mpl::int_<2>,boost::mpl::int_<6> > > res_type;
    BOOST_MPL_ASSERT( ( boost::is_same<decltype( v4 ),res_type> ) );
    assert( ( boost::is_same<decltype( v4 ),res_type >::value ) );


}
#if 0
F1<boost::fusion::vector<mpl_::int_<5>, mpl_::int_<6>, mpl_::int_<7> > >::result<
F1<boost::fusion::vector<mpl_::int_<5>, mpl_::int_<6>, mpl_::int_<7> > >( boost::fusion::vector<>, mpl_::int_<1> )>::type

{
    aka
    boost::fusion::joint_view<boost::fusion::vector<>,
    boost::fusion::joint_view<
    const boost::fusion::joint_view<
    const boost::fusion::joint_view<const boost::fusion::vector<>, const boost::fusion::single_view<boost::fusion::pair<mpl_::int_<1>, mpl_::int_<5> > > >,
    const boost::fusion::single_view<boost::fusion::pair<mpl_::int_<1>, mpl_::int_<6> > > >,
    const boost::fusion::single_view<boost::fusion::pair<mpl_::int_<1>, mpl_::int_<7> > > > >


    //to non-scalar type const State1 {aka

    const boost::fusion::joint_view<const boost::fusion::vector<>,
    boost::fusion::joint_view<
    const boost::fusion::joint_view<
    const boost::fusion::joint_view<const boost::fusion::vector<>,const boost::fusion::single_view<boost::fusion::pair<mpl_::int_<1>, mpl_::int_<5> > > >,
    const boost::fusion::single_view<boost::fusion::pair<mpl_::int_<1>, mpl_::int_<6> > > >,
    const boost::fusion::single_view<boost::fusion::pair<mpl_::int_<1>, mpl_::int_<7> > > > >


    boost::fusion::joint_view<boost::fusion::vector<>, const boost::fusion::single_view<boost::fusion::pair<mpl_::int_<1>, mpl_::int_<5> > > >


    const boost::fusion::joint_view<const boost::fusion::vector<>, const boost::fusion::single_view<boost::fusion::pair<mpl_::int_<1>, const mpl_::int_<5> > > >
#endif
