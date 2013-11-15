/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
struct has_value
{
    static const bool value = ( Contextv & Value ) != 0;
};

template<size_type Contextv, size_type Value>
struct set_value
{
    static const bool value = ( Contextv | Value );
};

template<size_type Contextv, size_type Value>
struct clear_value
{
    static const bool value = Contextv & ( ~Value );
};

/*!
  \class Context
 *\ingroup Core
 *\brief Context class

  @author Christophe Prud'homme
*/
class Context
{
public:

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     *
     * @param c context
     */
    explicit Context( size_type c )
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
    Context& operator=( size_type __c )
    {
        M_context = __c;
        return *this;
    }

    /**
     * get the context
     *
     * @return the context
     */
    size_type operator()() const
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
    size_type context() const
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
    void setContext( size_type __v )
    {
        M_context = __v ;
    }


    //@}

    /** @name  Methods
     */
    //@{


    bool test( size_type b ) const
    {
        return ( M_context&b )!=0;
    }
    template<typename T> bool test( T b ) const
    {
        return ( M_context&size_type( b ) )!=0;
    }
    void set( size_type b )
    {
        M_context |= b;
    }
    void set( size_type b, bool v );
    void clear( size_type b )
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

    size_type M_context;

};
}
#endif /* __Context_H */
