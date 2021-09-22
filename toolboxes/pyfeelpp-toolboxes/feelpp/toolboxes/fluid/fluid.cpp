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
#include <feel/feelmodels/fluid/fluidmechanics.hpp>

namespace py = pybind11;
using namespace Feel;

template<typename MeshT>
std::shared_ptr<MeshT>
loadmesh( std::string const& n, double h )
{
    return loadMesh( _mesh=new MeshT, _filename=n, _h=h );
}
template<int nDim, int OrderVelocity=2, int OrderPressure = 1, int OrderGeo = 1>
void defFM(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using fm_t = FluidMechanics< Simplex<nDim,OrderGeo>,
                                 Lagrange<OrderVelocity, Vectorial,Continuous,PointSetFekete>,
                                 Lagrange<OrderPressure, Scalar,Continuous,PointSetFekete> > ;
    std::string pyclass_name = std::string("Fluid_") + std::to_string(fm_t::mesh_type::nDim) + std::string("D") + std::string("P") + std::to_string(OrderVelocity) + std::string("P") + std::to_string(OrderPressure) + std::string("G") + std::to_string(OrderGeo);
    py::class_<fm_t,std::shared_ptr<fm_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,std::string const&,worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("keyword")=std::string("fluid"),
             py::arg("worldComm"),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize the fluid mechanics toolbox"
             )
        .def("init",static_cast<void (fm_t::*)(bool)>(&fm_t::init), "initialize the fluid mechanics toolbox",py::arg("buildModelAlgebraicFactory")= true)
        self_ptrtype const& other_toolbox, std::vector<std::string> const& markersInterpolate,
               bool buildModelAlgebraicFactory = true ); 
        .def("init",static_cast<void (fm_t::*)(std::shared_ptr<fm_t> const&, std::vector<std::string> const&, bool)>(&fm_t::init), "initialize another toolbox and a set of markers",
             py::arg("toolbox"),py::arg("markers"),py::arg("buildModelAlgebraicFactory")= true);
        // function spaces and elements
        .def("functionSpaceVelocity",&fm_t::functionSpaceVelocity, "get the velocity function space")
        .def("fieldVelocity",static_cast<typename fm_t::element_velocity_type& (fm_t::*)()>(&fm_t::fieldVelocity), "get the velocity field")

        .def("functionSpacePressure",&fm_t::functionSpacePressure, "get the pressure function space")
        .def("fieldPressure",static_cast<typename fm_t::element_pressure_type& (fm_t::*)()>(&fm_t::fieldPressure), "get the pressure field")

        // time stepping
        .def("timeStepBase",static_cast<std::shared_ptr<TSBase> (fm_t::*)() const>(&fm_t::timeStepBase), "get time stepping base")
        .def("updateTimeStep",&fm_t::updateTimeStep, "update time stepping")

        // solve
        .def("solve",&fm_t::solve, "solve the fluid mechanics problem")
        .def("exportResults",static_cast<void (fm_t::*)()>(&fm_t::exportResults), "export the results of the fluid mechanics problem")
        .def("exportResults",static_cast<void (fm_t::*)( double )>(&fm_t::exportResults), "export the results of the fluid mechanics problem", py::arg("time"))
        ;
        
}
    

PYBIND11_MODULE(_fluid, m )
{
    using namespace Feel;
    
    defFM<2,2,1,1>(m);
    defFM<2,3,2,1>(m);
    defFM<3,2,1,1>(m);
    defFM<3,3,2,1>(m);

}

