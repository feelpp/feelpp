/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013 Feel++ Consortium

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
   \file crbenums.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_CRBENUMS_HPP)
#define FEELPP_CRBENUMS_HPP 1

namespace Feel {

/**
 * CRBErrorType
 * Determine the type of error estimation used
 * - CRB_RESIDUAL : use the residual error estimation without algorithm SCM
 * - CRB_RESIDUAL_SCM : use the residual error estimation and also algorithm SCM
 * - CRB_NO_RESIDUAL : in this case we don't compute error estimation
 * - CRB_EMPIRICAL : compute |S_n - S_{n-1}| where S_n is the output obtained by using a reduced basis with n elements
 */
enum CRBErrorType { CRB_RESIDUAL = 0, CRB_RESIDUAL_SCM=1, CRB_NO_RESIDUAL=2 , CRB_EMPIRICAL=3};


}

#endif /* FEELPP_CRBENUMS_HPP */
