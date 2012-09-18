/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-12-14

  Copyright (C) 2006 Universit� Joseph Fourier (Grenoble)

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
   \file entities.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-12-14
 */
#include <boost/detail/identifier.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelmesh/entities.hpp>

namespace Feel
{
/// \cond detail
const uint16_type no_permutation::NO_PERMUTATION;
const uint16_type no_permutation::IDENTITY;
const uint16_type no_permutation::N_PERMUTATIONS;

const uint16_type line_permutations::NO_PERMUTATION;
const uint16_type line_permutations::IDENTITY;
const uint16_type line_permutations::REVERSE_PERMUTATION;
const uint16_type line_permutations::N_PERMUTATIONS;

const uint16_type triangular_faces_type::NO_PERMUTATION;
const uint16_type triangular_faces_type::IDENTITY;
const uint16_type triangular_faces_type::REVERSE_HEIGHT;
const uint16_type triangular_faces_type::REVERSE_BASE;
const uint16_type triangular_faces_type::REVERSE_HYPOTENUSE;
const uint16_type triangular_faces_type::ROTATION_ANTICLOCK;
const uint16_type triangular_faces_type::ROTATION_CLOCKWISE;
const uint16_type triangular_faces_type::N_PERMUTATIONS;
const uint16_type triangular_faces_type::PRINCIPAL_DIAGONAL;
const uint16_type triangular_faces_type::SECOND_DIAGONAL;
const uint16_type triangular_faces_type::ROTATION_TWICE_CLOCKWISE;

const uint16_type quadrangular_faces::NO_PERMUTATION;
const uint16_type quadrangular_faces::IDENTITY;
const uint16_type quadrangular_faces::SECOND_DIAGONAL;
const uint16_type quadrangular_faces::REVERSE_BASE;
const uint16_type quadrangular_faces::ROTATION_ANTICLOCK;
const uint16_type quadrangular_faces::PRINCIPAL_DIAGONAL;
const uint16_type quadrangular_faces::ROTATION_TWICE_CLOCKWISE;
const uint16_type quadrangular_faces::ROTATION_CLOCKWISE;
const uint16_type quadrangular_faces::REVERSE_HEIGHT;
const uint16_type quadrangular_faces::N_PERMUTATIONS;

/// \endcond
}
