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
#include <feel/feelmodels/solid/solidmechanics.hpp>

namespace py = pybind11;
using namespace Feel;

template<typename MeshT>
std::shared_ptr<MeshT>
loadmesh( std::string const& n, double h )
{
    return loadMesh( _mesh=new MeshT, _filename=n, _h=h );
}
template<int nDim, int OrderDisp>
void defSM(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using sm_t = SolidMechanics< Simplex<nDim,1>,
                                 Lagrange<OrderDisp, Vectorial,Continuous,PointSetFekete> >;
    using element_displacement_t = typename sm_t::element_displacement_type;
    using element_displacement_ptr_t = typename sm_t::element_displacement_ptrtype;
    using element_pressure_t = typename sm_t::element_pressure_type;
    using element_pressure_ptr_t = typename sm_t::element_pressure_ptrtype;
    
    std::string pyclass_name = std::string("Solid_") + std::to_string(nDim) + std::string("DP") + std::to_string(OrderDisp);
    py::class_<sm_t,std::shared_ptr<sm_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init<std::string const&,std::string const&,worldcomm_ptr_t const&,std::string const&, ModelBaseRepository const&>(),
             py::arg("prefix"),
             py::arg("keyword")=std::string("solid"),
             py::arg("worldComm"),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize the solid mechanics toolbox"
             )
        .def_static( "create", &sm_t::New,
                     py::arg("prefix"),
                     py::arg("keyword")=std::string("solid"),
                     py::arg("worldComm"),
                     py::arg("subprefix")=std::string(""),
                     py::arg("modelRep") = ModelBaseRepository(),
                     "Initialize the solid mechanics toolbox"
                     )
        .def("init",&sm_t::init, "initialize the solid mechanics toolbox",py::arg("buildModelAlgebraicFactory")= true)
        .def("mesh",static_cast<typename sm_t::mesh_ptrtype  (sm_t::*)() const>(&sm_t::mesh), "get the mesh")
        .def("setMesh",static_cast<void (sm_t::*)(typename sm_t::mesh_ptrtype const&)>(&sm_t::setMesh), "set the mesh of the toolbox",py::arg("mesh"))
        .def("updateParameterValues", &sm_t::updateParameterValues, "update parameter values" )

        .def("is1dReducedModel",&sm_t::is1dReducedModel, "returns true if 1D reduced model, false otherwise")
        .def("isStandardModel",&sm_t::isStandardModel, "returns true if standard model, false otherwise")
        
        // time stepping
        .def("timeStepBase",static_cast<std::shared_ptr<TSBase> (sm_t::*)() const>(&sm_t::timeStepBase), "get time stepping base")
        .def("startTimeStep", &sm_t::startTimeStep, "start time stepping" )
        .def("updateTimeStep",&sm_t::updateTimeStep, "update time stepping")

        // elements
        .def( "fieldDisplacement", static_cast<element_displacement_t& (sm_t::*)()>(&sm_t::fieldDisplacement), "returns the displacement field" )
        .def( "fieldDisplacementPtr", static_cast<element_displacement_ptr_t& (sm_t::*)()>(&sm_t::fieldDisplacementPtr), "returns the displacement field shared_ptr" )
        .def( "fieldPressure", static_cast<element_pressure_t& (sm_t::*)()>(&sm_t::fieldPressure), "returns the pressure field" )
        .def( "fieldPressurePtr", static_cast<element_pressure_ptr_t& (sm_t::*)()>(&sm_t::fieldPressurePtr), "returns the pressure field shared_ptr" )
        
        // solve
        .def("solve",&sm_t::solve, py::arg("upVelAcc")=true, "solve the solid mechanics problem, set boolean to true to update velocity and acceleration")
        .def("exportResults",static_cast<void (sm_t::*)()>(&sm_t::exportResults), "export the results of the solid mechanics problem")
        .def("exportResults",static_cast<void (sm_t::*)( double )>(&sm_t::exportResults), "export the results of the solid mechanics problem", py::arg("time"))
        ;
        
}
    

PYBIND11_MODULE(_solid, m )
{
    using namespace Feel;

    defSM<2,1>(m);
    defSM<2,2>(m);
    defSM<3,1>(m);
    defSM<3,2>(m);

}

