/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date     : Tue Feb 25 06:36:40 2014

   Copyright (C) 2014 Feel++ Consortium

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
#ifndef FEELPP_CREATE_ELEMENT_VECTOR_HPP
#define FEELPP_CREATE_ELEMENT_VECTOR_HPP 1

namespace Feel { namespace detail {

template<typename ElementType>
struct CreateElementVector
{
public:
    template<typename Sig>
    struct result;

    template<typename Lhs, typename Rhs>
    struct result<CreateElementVector( Lhs,Rhs )>
    {
	    typedef typename boost::remove_const<typename boost::remove_reference<Lhs>::type>::type lhs_noref_type;
	    typedef typename boost::remove_const<typename boost::remove_reference<Rhs>::type>::type::element_type rhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<ElementType>::type>::type ElementType_noref_type;
	    typedef typename boost::fusion::result_of::size<lhs_noref_type>::type index;
        typedef typename ElementType_noref_type::template sub_element<index::value>::type elt_type;
        BOOST_MPL_ASSERT( ( boost::is_same<typename elt_type::functionspace_type,rhs_noref_type> ) );
        typedef typename boost::fusion::result_of::make_vector<elt_type>::type v_elt_type;

	    typedef typename boost::fusion::result_of::push_back<lhs_noref_type, elt_type>::type ptype;
        typedef typename boost::fusion::result_of::as_vector<ptype>::type type;
    };
    CreateElementVector( ElementType const& e ) : M_e( e ), M_names() {}
    CreateElementVector( ElementType const& e, std::vector<std::string> const& names ) : M_e( e ), M_names( names ) {}
    ElementType const& M_e;
    std::vector<std::string> M_names;

    template<typename Lhs, typename Rhs>
    typename result<CreateElementVector( Lhs,Rhs )>::type
    operator()( Lhs const&  lhs, Rhs const& rhs ) const
    {
        typedef typename boost::remove_const<typename boost::remove_reference<Lhs>::type>::type lhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<Rhs>::type>::type rhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<ElementType>::type>::type ElementType_noref_type;
	    typedef typename boost::fusion::result_of::size<lhs_noref_type>::type index;
        typename ElementType_noref_type::template sub_element<index::value>::type elt = M_e.template element<index::value>();
        static const int s = mpl::size<typename ElementType::functionspace_type::bases_list>::type::value;
        BOOST_STATIC_ASSERT( (boost::is_same<decltype(elt), typename ElementType::template sub_element<index::value>::type>::value ) );
        if ( !M_names.empty() && M_names.size() > index::value )
        {

            FEELPP_ASSERT( M_names.size() == s  )
                ( M_names.size() )( s ).error( "incompatible number of function names and functions");
            elt.setName( M_names[index::value] );
        }
        else if  ( ( M_names.size() == 1 )  && s > 1 )
        {
            elt.setName( (boost::format( "%1%-%2%" ) % M_names[0] % index::value ).str() );
        }
        else
        {
            elt.setName( (boost::format( "%1%-%2%" ) % M_e.name() % index::value ).str() );
        }
        return boost::fusion::as_vector( boost::fusion::push_back( lhs, elt ) );
    }
};


} } // Feel / detail
#endif // FEELPP_CREATE_ELEMENT_VECTOR_HPP
