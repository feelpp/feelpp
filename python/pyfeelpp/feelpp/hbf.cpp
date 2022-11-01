//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! Copyright (C) 2017-present Feel++ Consortium
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
//! @date 01 Nov 2020
//! @copyright 2020 Feel++ Consortium
//!
#include <feel/feelpython/pybind11/pybind11.h>
#include <feel/feelpython/pybind11/stl.h>
#include <feel/feelpython/pybind11/eigen.h>

#include<feel/feelcore/environment.hpp>
#include<feel/feelfilters/hbf.hpp>

namespace py = pybind11;
void
defHBF( py::module & m )
{
    using namespace Feel;

    m.def( "readHBF",static_cast<holo3_image<float>(*)(std::string const&)>(&readHBF<float>),"read a HBF file of float", pybind11::return_value_policy::copy, py::arg("filename") );
    m.def( "readHBF",static_cast<holo3_image<float>(*)(std::istream&)>(&readHBF<float>),"read a HBF file of float", pybind11::return_value_policy::copy, py::arg("stream") );
    m.def( "writeHBF",static_cast<void (*)(std::string const&, holo3_image<float> const&, worldcomm_ptr_t)>(&writeHBF<float>),"write a HBF file of float", py::arg("filename"), py::arg("hbf"), py::arg("worldComm") );
    m.def( "writeHBF",static_cast<void (*)(std::ostream&, holo3_image<float> const&)>(&writeHBF<float>),"write a HBF file of float", py::arg("stream"), py::arg("hbf") );
    //m.def( "readHBF",&readHBF<double>,"read a HBF file of double", pybind11::return_value_policy::copy, py::arg("filename") );
    //m.def( "readHBF",&readHBF<int>,"read a HBF file of int", pybind11::return_value_policy::copy, py::arg("filename") );

    
}
