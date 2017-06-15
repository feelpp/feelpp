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
//! @date 14 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
#include <vector>

#include <feel/feelcrb/crbdata.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE( libfcrb )
{
    using namespace Feel;
    class_<CRBResults>("CRBResults")
        .def("output", &CRBResults::output)
        .def("errorBound", &CRBResults::errorbound)
        ;
    class_<ParameterSpaceX::Element>("ParameterSpaceElement")
        ;
    class_<ParameterSpaceX::Sampling>("ParameterSpaceSampling", init<boost::shared_ptr<ParameterSpaceX>,int,boost::shared_ptr<ParameterSpaceX::Sampling>>())
        ;
#if 1
    class_<ParameterSpaceX>("ParameterSpace")
        .def("sampling", &ParameterSpaceX::sampling)
        .def("element", &ParameterSpaceX::element)
        ;
#endif

}
