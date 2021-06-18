/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2020-05-04

   Copyright (C) 2020 Université de Strasbourg

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
#pragma once
#include <map>
#include <memory>
#include <utility>
#include <optional>
namespace Feel
{
/**
 * class bimap
 */
template<typename L, typename R>
class bimap
{
private:
    /* data */
public:
    bimap(/* args */);
    ~bimap();
    using left_t = L;
    using right_t = R;
    using value_type = std::pair<left_t, right_t>;
    using relation = value_type;

    bimap& operator=( std::map<left_t, right_t> const& );

    //! @return size of the bimap
    int size() const;

    //! insert a new element in the bimap
    void insert( std::pair<left_t, right_t>&& p );

    //! given an right index, get the left index
    std::optional<right_t> leftFind( left_t ) const;
    //! given an left index, get the right index
    std::optional<left_t> rightFind( right_t ) const;

    //! given an right index, get the left index
    right_t leftAt( left_t ) const;
    //! given an left index, get the right index
    left_t rightAt( right_t ) const;

private:
    class Pimpl;
    std::unique_ptr<Pimpl> pimpl_;
};

}