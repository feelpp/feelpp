/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-17

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file functordomain.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-17
 */
#ifndef __Functordomain_H
#define __Functordomain_H 1

#include <feel/feelcore/feel.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
/*!
  \class Functordomain
  \brief Concept to define the domain of definition of a functor

  @author Christophe Prud'homme
*/
template<typename T = double>
class FunctorDomain
{
public:
    typedef T value_type;

    FunctorDomain()
    {}

    virtual ~FunctorDomain()
    {}

    virtual bool hasLowerBound() const
    {
        return false;
    }

    virtual value_type lowerBound() const
    {
        FEELPP_ASSERT( true )( "FunctorDomain::lowerBound() called for a domain without "
                               "a lower bound" );
        return 0.0;
    }

    virtual bool hasUpperBound() const
    {
        return false;
    }

    virtual value_type upperBound() const
    {
        FEELPP_ASSERT( true )( "FunctorDomain::upperBound() called for a domain without "
                               "a upper bound" );
        return 0.0;
    }

    virtual bool hasExcludedPoint() const
    {
        return false;
    }

    virtual value_type excludedPoint() const
    {
        FEELPP_ASSERT( true )( "FunctorDomain::excludedPoint() called for a domain without "
                               "an excluded point" );
        return 0.0;
    }

};
template<typename T = double>
class UnboundedDomain : public FunctorDomain<T>
{
    typedef FunctorDomain<T> super;
public:


    typedef typename super::value_type value_type;

    UnboundedDomain()
        :
        super()
    {}
};

template<typename T = double>
class PositiveDomain : public FunctorDomain<T>
{
    typedef FunctorDomain<T> super;

public:

    typedef typename super::value_type value_type;

    PositiveDomain()
        :
        super()
    {}

    virtual bool hasLowerBound() const
    {
        return true;
    }

    virtual value_type lowerBound() const
    {
        return 0.0;
    }
};

template<typename T = double>
class BoundedDomain : public FunctorDomain<T>
{
    typedef FunctorDomain<T> super;

public:

    typedef typename super::value_type value_type;

    BoundedDomain( const value_type& lower, const value_type& upper )
        :
        super(),
        lower_( lower ),
        upper_( upper )
    {}

    virtual bool hasLowerBound() const
    {
        return true;
    }

    virtual value_type lowerBound() const
    {
        return lower_;
    }

    virtual bool hasUpperBound() const
    {
        return true;
    }

    virtual value_type upperBound() const
    {
        return upper_;
    }

private:
    value_type lower_;

    value_type upper_;
};

template<typename T = double>
class LowerBoundedDomain : public FunctorDomain<T>
{
    typedef FunctorDomain<T> super;

public:

    typedef typename super::value_type value_type;

    LowerBoundedDomain( const value_type& lower )
        :
        super(),
        lower_( lower )
    {}

    virtual bool hasLowerBound() const
    {
        return true;
    }

    virtual value_type lowerBound() const
    {
        return lower_;
    }

private:
    value_type lower_;
};

template<typename T = double>
class NonzeroDomain : public FunctorDomain<T>
{
    typedef FunctorDomain<T> super;

public:

    typedef typename super::value_type value_type;

    NonzeroDomain()
        :
        super()
    {}

    virtual bool hasExcludedPoint() const
    {
        return true;
    }

    virtual value_type excludedPoint() const
    {
        return 0.0;
    }
};
/// \endcond
}
}
#endif /* __Functordomain_H */
