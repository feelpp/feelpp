/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2004-09-10

   Copyright (C) 2006,2009 Universite Joseph Fourier (Grenoble I)
   Copyright (C) 2004 EPFL
   Copyright (C) 2011-2016 Feel++ Consortium

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
   \file policy.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-09-10
 */
#ifndef FEELPP_FEELCORE_POLICY_HPP
#define FEELPP_FEELCORE_POLICY_HPP 1

namespace Feel
{
/**
 * \class PolicyCreationUsingNew
 *  default creation policy
 *
 *  \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
template <class T>
struct PolicyCreationUsingNew
{
    static T* create()
    {
        return new T;
    }

    static void destroy( T* p )
    {
        delete p;
    }
};

/**
 * \class PolicyFeelTimeDefault
 *  default feel time policy
 *
 *  \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
template <class T>
struct PolicyFeelTimeDefault
{
    static void scheduleDestruction( T*, void ( *pFun )() )
    {
        std::atexit( pFun );
    }

    static void onDeadReference()
    {
        throw std::logic_error( "Dead Reference Detected" );
    }
};
}

#endif /* FEELPP_FEELCORE_POLICY_HPP */
