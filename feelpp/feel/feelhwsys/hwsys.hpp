// This file is part of the Feel library
//
// Author(s): Feel++ Contortium
//      Date: 2017-07-10
//
// @copyright (C) 2017 University of Strasbourg
// @copyright (C) 2012-2017 Feel++ Consortium
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef FEELPP_HWSYS_HPP
#define FEELPP_HWSYS_HPP 1

#include <feel/feelhwsys/hwsysbase.hpp>

#if defined(FEELPP_HAS_KWSYS )
#include <feel/feelhwsys/kwsys.hpp>
#endif

namespace Feel
{
namespace Sys
{

//! Feel++ dedicated hardware and OS information class
//! (no external libraries).
class HWSys
    : public HwSysBase
{
    HWSys() = default;

    ~HWSys() override = default;
};

} // Sys namespace
} // Feel namespace

#endif // FEELPP_HWSYS_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
