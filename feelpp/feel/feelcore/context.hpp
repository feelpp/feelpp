/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-03-08

  Copyright (C) 2009 Université de Grenoble 1
  Copyright (C) 2005,2006 EPFL

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
   \file context.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-03-08
 */
#ifndef __Context_H
#define __Context_H 1

#include <feel/feelcore/feel.hpp>

namespace Feel
{

template<size_type Contextv, size_type Value>
using has_value = mpl::bool_<( Contextv & Value ) != 0>;

template<size_type Contextv, size_type Value>
constexpr bool has_value_v = has_value<Contextv,Value>::value;

template<size_type Contextv, size_type Value>
using set_value = mpl::size_t<( Contextv | Value )>;

template<size_type Contextv, size_type Value>
constexpr size_type set_value_v = set_value<Contextv,Value>::value;

template<size_type Contextv, size_type Value>
using clear_value = mpl::size_t<Contextv & ( ~Value )>;

template<size_type Contextv, size_type Value>
constexpr size_type clear_value_v = clear_value<Contextv,Value>::value;

template<size_type Contextv, size_type Value1, size_type... Values>
struct clear_values{
    using type = typename clear_values< clear_value_v<Contextv,Value1>, Values...>::type;
};
template<size_type Contextv, size_type Value1>
struct clear_values<Contextv,Value1> {
    using type = clear_value<Contextv,Value1>;
};

template<size_type Contextv, size_type... Values>
constexpr size_type clear_values_v = clear_values<Contextv,Values...>::type::value;


template<size_type Contextv>
inline bool hasValue( size_type c ) { return (c & Contextv) != 0; }

inline bool hasValue( size_type c, size_type v ) { return (c & v) != 0; }

namespace meta
{
/*!
  \class Context
 *\ingroup Core
 *\brief Context class

  @author Christophe Prud'homme
*/
template <typename StorageType>
class Context
{
public:
    typedef StorageType storage_type;

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     *
     * @param c context
     */
    explicit Context( storage_type c )
        :
        M_context( c )
    {}

    /**
     * copy constructor
     *
     */
    Context( Context const & c )
        :
        M_context( c.M_context )
    {}

    /**
     * destructor
     *
     */
    ~Context()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     * assignment operator
     *
     * @return the context
     */
    Context& operator=( Context const& __c )
    {
        if (  this != &__c )
        {
            M_context = __c.M_context;
        }

        return *this;
    }

    /**
     * assignment operator
     *
     * @return the context
     */
    Context& operator=( storage_type __c )
    {
        M_context = __c;
        return *this;
    }

    /**
     * get the context
     *
     * @return the context
     */
    storage_type operator()() const
    {
        return M_context;
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * get context value
     *
     * @return the context value
     */
    storage_type context() const
    {
        return M_context;
    }


    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set context value
     *
     * @param __v context value
     */
    void setContext( storage_type __v )
    {
        M_context = __v ;
    }


    //@}

    /** @name  Methods
     */
    //@{


    bool test( storage_type b ) const
    {
        return ( M_context&b )!=0;
    }
    template<typename T> bool test( T b ) const
    {
        return ( M_context&storage_type( b ) )!=0;
    }
    void set( storage_type b )
    {
        M_context |= b;
    }
    void set( storage_type b, bool v );
    void clear( storage_type b )
    {
        M_context &= ( uint )( ~b );
    }
    void reset()
    {
        M_context = 0;
    }

    //@}

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & M_context;
        }

private:

    storage_type M_context;

};

} // namespace meta

using Context = meta::Context<size_type>;
}
#endif /* __Context_H */
