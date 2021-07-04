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
#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feeldiscr/pch.hpp>

namespace py = pybind11;
using namespace Feel;

template<int nDim, int Order>
void defSM(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using toolbox_t = Heat< Simplex<nDim,1,nDim>,
                            Lagrange<Order, Scalar,Continuous,PointSetFekete,0> >;
    using element_temperature_t = typename toolbox_t::element_temperature_type;
    using element_temperature_ptr_t = typename toolbox_t::element_temperature_ptrtype;
    // using element_velocityconvection_t = typename toolbox_t::element_velocityconvection_type;
    // using element_velocityconvection_ptr_t = typename toolbox_t::element_velocityconvection_ptrtype;
    using space_temperature_ptr_t = typename toolbox_t::space_temperature_ptrtype;

    std::string pyclass_name = std::string("Heat_") + std::to_string(nDim) + std::string("DP") + std::to_string(Order);
    py::class_<toolbox_t, std::shared_ptr<toolbox_t>, ModelNumerical>( m, pyclass_name.c_str() )
        .def( py::init<std::string const&, std::string const&, worldcomm_ptr_t const&, std::string const&, ModelBaseRepository const&>(),
              py::arg( "prefix" ),
              py::arg( "keyword" ) = std::string( "heat" ),
              py::arg( "worldComm" ) = Environment::worldCommPtr(),
              py::arg( "subprefix" ) = std::string( "" ),
              py::arg( "modelRep" ) = ModelBaseRepository(),
              "Initialize the heat mechanics toolbox" )
        .def( "init", &toolbox_t::init, "initialize the heat  toolbox", py::arg( "buildModelAlgebraicFactory" ) = true )

        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        .def( "rangeMeshElements", &toolbox_t::rangeMeshElements, "get the range of mesh elements" )

        // Algrebrzic factory
        .def( "algebraicFactory", static_cast<typename toolbox_t::model_algebraic_factory_ptrtype const&  (toolbox_t::*)() const>(&toolbox_t::algebraicFactory), "get the algebraic factory" )
#if 0
        .def( "matrix", &toolbox_t::matrix, "matrix used for assembly linear system or the jacobian" )
        .def( "rhs", &toolbox_t::rhs, "vector used for assembly rhs in linear system or the residual" )
        .def( "backend", &ModelAlgebraic::backend, "get the backend" )
        .def( "model", &ModelAlgebraic::model, "get the model" )
#endif        
        // elements
        .def( "spaceTemperature", &toolbox_t::spaceTemperature, "get the temperature function space" )
        .def( "fieldTemperature", static_cast<element_temperature_t const& (toolbox_t::*)() const>( &toolbox_t::fieldTemperature ), "returns the temperature field" )
        .def( "fieldTemperaturePtr", static_cast<element_temperature_ptr_t const& (toolbox_t::*)() const>( &toolbox_t::fieldTemperaturePtr ), "returns the temperature field shared_ptr" )
        //.def( "spaceVelocityConvection", &toolbox_t::spaceVelocityConvection, "get the field function space")
        // .def( "fieldVelocityConvection", static_cast<element_velocityconvection_t const& (toolbox_t::*)() const>(&toolbox_t::fieldVelocityConvection), "returns the convection velocity field" )
        // .def( "fieldVelocityConvectionPtr", static_cast<element_velocityconvection_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldVelocityConvectionPtr), "returns the convection velocity field shared_ptr" )

        // time stepping
        .def( "timeStepBase", static_cast<std::shared_ptr<TSBase> ( toolbox_t::* )() const>( &toolbox_t::timeStepBase ), "get time stepping base" )
        .def( "updateTimeStep", &toolbox_t::updateTimeStep, "update time stepping" )

        // solve
        .def( "solve", &toolbox_t::solve, "solve the heat mechanics problem, set boolean to true to update velocity and acceleration" )
        .def( "exportResults", static_cast<void ( toolbox_t::* )()>( &toolbox_t::exportResults ), "export the results of the heat mechanics problem" )
        .def( "exportResults", static_cast<void ( toolbox_t::* )( double )>( &toolbox_t::exportResults ), "export the results of the heat mechanics problem", py::arg( "time" ) )

        .def( "setMesh", &toolbox_t::setMesh, "set the mesh", py::arg( "mesh" ) )
        .def( "updateParameterValues", &toolbox_t::updateParameterValues, "update parameter values" )
        .def(
            "assembleMatrix", []( const toolbox_t& t )
            {
                auto mat = t.algebraicFactory()->matrix(); //->clone();
                mat->zero();
                auto rhsTMP = t.algebraicFactory()->rhs()->clone();
                t.algebraicFactory()->applyAssemblyLinear( t.algebraicBlockVectorSolution()->vectorMonolithic(), mat, rhsTMP, { "ignore-assembly.rhs" } );
                return mat;
            },
            "returns the assembled matrix" )
        .def(
            "assembleRhs", []( const toolbox_t& t )
            {
                auto rhs = t.algebraicFactory()->rhs()->clone();
                rhs->zero();
                auto matTMP = t.algebraicFactory()->matrix()->clone();
                t.algebraicFactory()->applyAssemblyLinear( t.algebraicBlockVectorSolution()->vectorMonolithic(), matTMP, rhs, { "ignore-assembly.lhs" } );
                return rhs;
            },
            "returns the assembled rhs" )
        // .def( "updateFieldVelocityConvection", static_cast<void (toolbox_t::*)(bool)>(&toolbox_t::updateFieldVelocityConvection), "update field velocity convection", py::arg("onlyExprWithTimeSymbol")=false )
        // .def( "assembleLinear", &toolbox_t::assembleLinear, "assemble linear matrix and vector" )
        // .def( "rhs", &toolbox_t::rhs, "returns the right hand side" )
        // .def( "matrix", &toolbox_t::matrix, "returns the matrix" )
        ;
}


PYBIND11_MODULE(_heat, m )
{
    using namespace Feel;

    defSM<2,1>(m);
    defSM<2,2>(m);
    defSM<3,1>(m);
    defSM<3,2>(m);

}

