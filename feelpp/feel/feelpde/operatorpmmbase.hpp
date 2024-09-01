/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date:18 Dec 2017

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

#ifndef FEELPP_OPERATORPMMBASE_HPP
#define FEELPP_OPERATORPMMBASE_HPP 1

#include <feel/feelalg/operator.hpp>

namespace Feel
{

template<typename T>
class OperatorPMMBase : public OperatorBase<T>
{
    typedef OperatorBase<T> super;
public :
    using datamap_type = typename super::datamap_type;
    using datamap_ptrtype = typename super::datamap_ptrtype;
    OperatorPMMBase( datamap_ptrtype const& map, std::string label, bool use_transpose, bool has_norminf )
        :
        super( map, label, use_transpose, has_norminf )
        {}

    virtual sparse_matrix_ptrtype pressureMassMatrix() const = 0;

};

}

#endif
