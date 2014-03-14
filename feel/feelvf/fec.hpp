/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-03-13

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
   \file fec.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-03-13
 */
#ifndef __FEC_H
#define __FEC_H 1

namespace Feel
{
namespace vf
{
/// \cond detail
namespace detail
{
template<uint16_type Type, typename FormContextType>
struct FEContextInit
{
    // 0 : test, 1 : trial
    static const uint16_type type = Type;
    typedef typename FormContextType::test_fe_type test_fe_type;
    typedef typename FormContextType::trial_fe_type trial_fe_type;


    typedef typename mpl::if_<mpl::equal_to<mpl::int_<type>, mpl::int_<0> >, mpl::identity<test_fe_type>, mpl::identity<trial_fe_type> >::type::type fe_type;
    typedef boost::shared_ptr<fe_type> fe_ptrtype;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<type>, mpl::int_<0> >, mpl::identity<typename FormContextType::test_fecontext_type>,mpl::identity<typename FormContextType::trial_fecontext_type> >::type::type fecontext_type;
    typedef boost::shared_ptr<fecontext_type> fecontext_ptrtype;

    typedef typename mpl::if_<mpl::equal_to<mpl::int_<type>, mpl::int_<0> >,
            mpl::identity<typename FormContextType::test_geometric_mapping_context_ptrtype>,
            mpl::identity<typename FormContextType::trial_geometric_mapping_context_ptrtype> >::type::type geometric_mapping_context_ptrtype;

    //typedef typename FormContextType::test_geometric_mapping_context_ptrtype geometric_mapping_context_ptrtype;

    //typedef typename FormContextType::form_type form_type;
    typedef FormContextType form_type;
    typedef boost::shared_ptr<form_type> form_ptrtype;

    template<typename Sig>
    struct result;

    template<typename T>
    struct result<FEContextInit( T )>
    {
        typedef fusion::pair<typename boost::remove_reference<T>::type::first_type,fecontext_ptrtype> type;
    };

    FEContextInit( fe_ptrtype const& fe, FormContextType const& form )
        :
        M_fe( fe ),
        M_form( form )
    {}
    template<typename T>
    fusion::pair<typename boost::remove_reference<T>::type::first_type,fecontext_ptrtype>
    operator()( T const& t ) const
    {
        return operator()( t, mpl::int_<type>() );
    }

private:

    // Test FE context
    template<typename T>
    fusion::pair<typename boost::remove_reference<T>::type::first_type,fecontext_ptrtype>
    operator()( T const& t, mpl::int_<0> ) const
    {
        geometric_mapping_context_ptrtype gmcptr( t.second );
        typedef typename boost::remove_reference<T>::type::first_type first_type;
        return fusion::make_pair<first_type>( fecontext_ptrtype( new fecontext_type( M_fe,
                                              gmcptr,
                                              M_form.testPc( gmcptr->faceId(), gmcptr->permutation() ) ) ) );
    }

    // Trial FE context
    template<typename T>
    fusion::pair<typename boost::remove_reference<T>::type::first_type,fecontext_ptrtype>
    operator()( T const& t, mpl::int_<1> ) const
    {
        geometric_mapping_context_ptrtype gmcptr( t.second );
        typedef typename boost::remove_reference<T>::type::first_type first_type;
        return fusion::make_pair<first_type>( fecontext_ptrtype( new fecontext_type( M_fe,
                                              gmcptr,
                                              M_form.trialPc( gmcptr->faceId(), gmcptr->permutation() ) ) ) );
    }

    fe_ptrtype const& M_fe;
    form_type const& M_form;
};


template<uint16_type Type, typename FormContextType>
struct FEContextUpdate
{
    // 0 : test, 1 : trial
    static const uint16_type type = Type;

    //typedef typename FormContextType::test_geometric_mapping_context_ptrtype geometric_mapping_context_ptrtype;
    //typedef typename FormContextType::map_test_geometric_mapping_context_type map_geometric_mapping_context_type;
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<type>, mpl::int_<0> >,
            mpl::identity<typename FormContextType::test_geometric_mapping_context_ptrtype>,
            mpl::identity<typename FormContextType::trial_geometric_mapping_context_ptrtype> >::type::type geometric_mapping_context_ptrtype;

    typedef typename mpl::if_<mpl::equal_to<mpl::int_<type>, mpl::int_<0> >,
            mpl::identity<typename FormContextType::map_test_geometric_mapping_context_type>,
            mpl::identity<typename FormContextType::map_trial_geometric_mapping_context_type> >::type::type map_geometric_mapping_context_type;


    //typedef typename FormContextType::form_type form_type;
    typedef FormContextType form_type;
    typedef boost::shared_ptr<form_type> form_ptrtype;

    FEContextUpdate( map_geometric_mapping_context_type const& mapgmc,
                     form_type const& form )
        :
        M_mapgmc( mapgmc ),
        M_form( form )
    {}

    template<typename T>
    void operator()( T& t ) const
    {
        return operator()( t, mpl::int_<type>() );
    }
    template<typename T>
    void operator()( T& t, mpl::int_<0> ) const
    {
        typedef typename boost::remove_reference<T>::type::first_type first_type;
        geometric_mapping_context_ptrtype gmcptr( fusion::at_key<first_type>( M_mapgmc ) );
        t.second->update( gmcptr, M_form.testPc( gmcptr->faceId(), gmcptr->permutation()  ) );
    }
    template<typename T>
    void operator()( T& t, mpl::int_<1> ) const
    {
        typedef typename boost::remove_reference<T>::type::first_type first_type;
        geometric_mapping_context_ptrtype gmcptr( fusion::at_key<first_type>( M_mapgmc ) );
        t.second->update( gmcptr, M_form.trialPc( gmcptr->faceId(), gmcptr->permutation()  ) );
    }

    map_geometric_mapping_context_type const& M_mapgmc;
    form_type const& M_form;
};


}
/// \endcond
} // vf
}
#endif /* __FEC_H */
