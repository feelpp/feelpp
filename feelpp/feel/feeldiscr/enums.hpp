/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2020-02-04

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2020 Feel++ Consortium

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
#ifndef FEELPP_DISCR_ENUMS_HPP
#define FEELPP_DISCR_ENUMS_HPP 1

namespace  Feel
{

enum class ComponentType
{
    NO_COMPONENT = -1,
    X = 0,
    Y,
    Z,
    NX,
    NY,
    NZ,
    TX,
    TY,
    TZ
};
using Component = ComponentType;
enum FunctionSpaceType
{
    SCALAR = 0,
    VECTORIAL = 1,
    TENSOR2,
    TENSOR2_SYMM
};

}

#endif
