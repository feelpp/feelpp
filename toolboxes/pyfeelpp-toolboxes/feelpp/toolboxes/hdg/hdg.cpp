//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
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
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 25 Jul 2018
//! @copyright 2018 Feel++ Consortium
//!
#include <pybind11/pybind11.h>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include "hdg-poisson.hpp"

namespace py = pybind11;
using namespace Feel;

    

PYBIND11_MODULE(_hdg, m )
{
    using namespace Feel;

    m.def( "hdg_poisson_options", &FeelModels::makeMixedPoissonOptions, "get hdg mixed Poisson command line options",
           py::arg("prefix")=std::string(""),py::arg("prefix_toolbox")=std::string("hdg.poisson"));
    
    defHDGPoisson<2,1>(m);
    defHDGPoisson<2,2>(m);
    defHDGPoisson<3,1>(m);
    defHDGPoisson<3,2>(m);
}

