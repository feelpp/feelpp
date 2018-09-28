/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 18 July 2017

 Copyright (C) 2017 Feel++ Consortium

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

#ifndef FEELPP_HASHTABLES_HPP
#define FEELPP_HASHTABLES_HPP 1

#include <boost/functional/hash.hpp>

namespace Feel {

namespace HashTables
{
template <typename T>
struct HasherContainers
{
    size_t operator()(const std::set<T>& v) const
        {
            return boost::hash_range( v.begin(),v.end() );
        }
    size_t operator()(const std::vector<T>& v) const
        {
            return boost::hash_range( v.begin(),v.end() );
        }
};

}

}

#endif
