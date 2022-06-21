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
#include <feel/feelmodels/fsi/fsi.hpp>
#include <mpi4py/mpi4py.h>

namespace py = pybind11;
using namespace Feel;

template<int nDim, int OrderPotential>
void defToolbox(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using toolbox_t = Fsi< Simplex<nDim,1>,
                           Lagrange<OrderPotential, Scalar,Continuous,PointSetFekete> >;
    using element_fsipotential_t = typename toolbox_t::element_fsipotential_type;
    using element_fsipotential_ptr_t = typename toolbox_t::element_fsipotential_ptrtype;
    using element_fsifield_t = typename toolbox_t::element_fsifield_type;
    using element_fsifield_ptr_t = typename toolbox_t::element_fsifield_ptrtype;

    std::string pyclass_name = std::string("Fsi_") + std::to_string(nDim) + std::string("DP") + std::to_string(OrderPotential);
    py::class_<toolbox_t,std::shared_ptr<toolbox_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,bool,worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("buildmesh")=true,
             py::arg("worldComm")=Environment::worldCommPtr(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize the fsi mechanics toolbox"
             )
        .def("init",&toolbox_t::init, "initialize the fsi mechanics toolbox",py::arg("buildModelAlgebraicFactory")= true)

        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        .def( "rangeMeshElements", &toolbox_t::rangeMeshElements, "get the range of mesh elements" )
        
        // elements
        .def( "spaceFsiPotential", &toolbox_t::spaceFsiPotential, "get the potential function space")
        .def( "fieldFsiPotential", static_cast<element_fsipotential_t const& (toolbox_t::*)() const>(&toolbox_t::fieldFsiPotential), "returns the fsi potential field" )
        .def( "fieldFsiPotentialPtr", static_cast<element_fsipotential_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldFsiPotentialPtr), "returns the fsi potential field shared_ptr" )
        .def( "spaceFsiField", &toolbox_t::spaceFsiField, "get the field function space")
        .def( "fieldFsiField", static_cast<element_fsifield_t const& (toolbox_t::*)() const>(&toolbox_t::fieldFsiField), "returns the fsi field" )
        .def( "fieldFsiFieldPtr", static_cast<element_fsifield_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldFsiFieldPtr), "returns the fsi field shared_ptr" )
        .def( "fieldCurrentDensity", static_cast<element_fsifield_t const& (toolbox_t::*)() const>(&toolbox_t::fieldCurrentDensity), "returns the current density" )
        .def( "fieldCurrentDensityPtr", static_cast<element_fsifield_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldCurrentDensityPtr), "returns the current density shared_ptr" )
        
        // solve
        .def("solve",&toolbox_t::solve, "solve the fsi mechanics problem, set boolean to true to update velocity and acceleration")
        .def("exportResults",static_cast<void (toolbox_t::*)()>(&toolbox_t::exportResults), "export the results of the fsi mechanics problem")
        .def("exportResults",static_cast<void (toolbox_t::*)( double )>(&toolbox_t::exportResults), "export the results of the fsi mechanics problem", py::arg("time"))
        ;
        
}
    

PYBIND11_MODULE(_fsi, m )
{
    using namespace Feel;

    if (import_mpi4py()<0) return ;

    
    
    defToolbox<2,1>(m);
    defToolbox<2,2>(m);
    defToolbox<3,1>(m);
    defToolbox<3,2>(m);

}

