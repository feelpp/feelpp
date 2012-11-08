/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <prudhomm@mit.edu>
   Date: 2004-09-10

   Copyright (C) 2006,2009 Université Joseph Fourier (Grenoble I)
   Copyright (C) 2004 EPFL

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
   \file singleton.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-09-10
 */
#ifndef __Singleton_H
#define __Singleton_H 1

#include <cstdlib>
#include <cassert>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <new>

#include <feel/feelcore/policy.hpp>

namespace Feel
{
/**
   \class Singleton
  *\ingroup Core
  *\brief implement the Singleton pattern

   A Singleton pattern implementation using the ideas
   from Alexandrescu's book "modern C++ design"
   http://www.moderncppdesign.com/


   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
*/
template <typename T>
class Singleton
{
public:

    typedef T singleton_type;

    typedef PolicyFeelTimeDefault<singleton_type> feeltime_policy;
    typedef PolicyCreationUsingNew<singleton_type> creation_policy;


    /**
     * return the instance of the singleton
     */
    static singleton_type& instance();

private:
    // Helpers
    static void makeInstance();

    /**
       Singleton::makeInstance (helper for Instance)
    */
    static void destroySingleton();

    // Protection
    Singleton();

    // Data
    typedef singleton_type* instance_type;
    static instance_type _S_instance;
    static bool _S_destroyed;
};

//
// Instantiate Singleton's static data
//
template <class T>
typename Singleton<T>::instance_type Singleton<T>::_S_instance;

template <class T>
bool Singleton<T>::_S_destroyed;


template <class T>
inline
T&
Singleton<T>::instance()
{
    if ( !_S_instance )
    {
        makeInstance();
    }

    return *_S_instance;
}


template <class T>
void
Singleton<T>::makeInstance()
{
    if ( !_S_instance )
    {
        if ( _S_destroyed )
        {
            feeltime_policy::onDeadReference();
            _S_destroyed = false;
        }

        _S_instance = creation_policy::create();
        feeltime_policy::scheduleDestruction( _S_instance, &destroySingleton );
    }
}

template <class T>
void
Singleton<T>::destroySingleton()
{
    assert( !_S_destroyed );
    creation_policy::destroy( _S_instance );
    _S_instance = 0;
    _S_destroyed = true;
}
}
#endif /* __singleton_HPP */
