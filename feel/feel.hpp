//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! Includes all the Feel++ library components
//!
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 05 Feb 2017
//! @copyright 2017 Feel++ Consortium
//!
#if !defined(FEELPP_FEEL_HPP)
#define FEELPP_FEEL_HPP 1

#include <boost/math/constants/constants.hpp>

#include <feel/options.hpp>

#include <feel/feelcore/core.hpp>

#include <feel/feelalg/alg.hpp>

#include <feel/feelpoly/poly.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelts/ts.hpp>

#include <feel/feeldiscr/discr.hpp>

#include <feel/feelpde/pde.hpp>

#include <feel/feeltiming/tic.hpp>

#include <ginac/ginac.h>
namespace Feel
{
using GiNaC::symbol;
using GiNaC::ex;
}

#include <feel/feelfilters/filters.hpp>

#endif /* FEELPP_FEEL_HPP */
