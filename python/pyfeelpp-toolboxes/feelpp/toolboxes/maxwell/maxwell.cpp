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
#include <feel/feelpython/pybind11/pybind11.h>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/maxwell/maxwell.hpp>

namespace py = pybind11;
using namespace Feel;

template<typename ConvexT, int Order=1>
void defToolbox(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using convex_t = ConvexT;
    using toolbox_t = Maxwell<convex_t>;
    using space_magneticpotential_t = typename toolbox_t::space_magneticpotential_type;
    using space_magneticpotential_ptr_t = typename toolbox_t::space_magneticpotential_ptrtype;
    using element_magneticpotential_t = typename toolbox_t::element_magneticpotential_type;
    using element_magneticpotential_ptr_t = typename toolbox_t::element_magneticpotential_ptrtype;
    using element_magneticfield_t = typename toolbox_t::element_magneticfield_type;
    using element_magneticfield_ptr_t = typename toolbox_t::element_magneticfield_ptrtype;

    std::string pyclass_name = std::string("Maxwell_") + std::to_string(nDim) + std::string("DP") + std::to_string(OrderPotential);
    py::class_<toolbox_t,std::shared_ptr<toolbox_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,bool,worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("buildmesh")=true,
             py::arg("worldComm")=Environment::worldCommPtr(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize the maxwell mechanics toolbox"
             )
        .def("init",&toolbox_t::init, py::arg("buildModelAlgebraicFactory")= true, "initialize the maxwell mechanics toolbox")

        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        .def( "rangeMeshElements", &toolbox_t::rangeMeshElements, "get the range of mesh elements" )
        
        // elements
        .def( "spaceMagneticpotential", &toolbox_t::spaceMagneticpotential, "get the potential function space")
        .def( "fieldMagneticpotential", static_cast<element_magneticpotential_t const& (toolbox_t::*)() const>(&toolbox_t::fieldMagneticpotential), "returns the maxwell potential field" )
        .def( "fieldMagneticpotentialPtr", static_cast<element_magneticpotential_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldMagneticpotentialPtr), "returns the maxwell potential field shared_ptr" )
        .def( "spaceMagneticfield", &toolbox_t::spaceMagneticfield, "get the field function space")
        .def( "fieldMagneticfield", static_cast<element_magneticfield_t const& (toolbox_t::*)() const>(&toolbox_t::fieldMagneticfield), "returns the maxwell field" )
        .def( "fieldMagneticfieldPtr", static_cast<element_magneticfield_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldMagneticfieldPtr), "returns the maxwell field shared_ptr" )

        // solve
        .def("solve",&toolbox_t::solve, "solve the maxwell mechanics problem, set boolean to true to update velocity and acceleration")
        .def("exportResults",static_cast<void (toolbox_t::*)()>(&toolbox_t::exportResults), "export the results of the maxwell mechanics problem")
        .def("exportResults",static_cast<void (toolbox_t::*)( double )>(&toolbox_t::exportResults), "export the results of the maxwell mechanics problem", py::arg("time"))
        ;
        
}
    

PYBIND11_MODULE(_maxwell, m )
{
    using namespace Feel;

    defToolbox<2,1>(m);
    defToolbox<3,1>(m);
}

