/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-03-20

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
   \file timer.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-03-20
 */
#if !defined( FEELPP_TIMING_TIMER_HPP)
#define FEELPP_TIMING_TIMER_HPP 1


#include <stack>
#include <boost/assert.hpp>

namespace Feel
{
namespace details
{
//////////////////////////////////////////////////////////////////////////////
// counter<R,T> is an implementation detail that gather and store cycles or
// seconds measures between tic and toc calls.
//////////////////////////////////////////////////////////////////////////////
template<class R,class T> class counter
{
public :
    typedef T timer_type;
    typedef R type;

    void  tic() const
    {
        times().push( time() );
    }

    type toc( std::string const& msg, bool display ) const
    {
        BOOST_ASSERT_MSG( !empty(), "Unbalanced timing calls" );
        type t = time()-times().top();
        times().pop();

        if ( display ) timer_type::print( msg, t );

        return t;
    }

    bool  empty() const
    {
        return times().empty();
    }
    type  time()  const
    {
        return timer_type::time();
    }

    std::stack<type>& times() const
    {
        static std::stack<type> local;
        return local;
    }
};
}
}

#endif /* FEELPP_TIMING_TIMER_HPP */
