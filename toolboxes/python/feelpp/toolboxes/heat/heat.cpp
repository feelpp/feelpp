/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Thomas Saigre <saigre@math.unistra.fr>
 * SPDX-FileContributor: Romain Hild <romain.hild@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2018-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#include <feel/feelpython/pybind11/pybind11.h>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feeldiscr/pch.hpp>

namespace py = pybind11;
using namespace Feel;

template<int nDim>
void defSM(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using toolbox_t = Heat< Simplex<nDim,Dynamic,nDim>,
                            Lagrange<Dynamic, Scalar,Continuous,PointSetFekete,0> >;
    using element_temperature_t = typename toolbox_t::element_temperature_type;
    using element_temperature_ptr_t = typename toolbox_t::element_temperature_ptrtype;
    // using element_velocityconvection_t = typename toolbox_t::element_velocityconvection_type;
    // using element_velocityconvection_ptr_t = typename toolbox_t::element_velocityconvection_ptrtype;
    using space_temperature_ptr_t = typename toolbox_t::space_temperature_ptrtype;

    std::string pyclass_name = std::string("Heat_") + std::to_string(nDim) + "DDynamic";
    py::class_<toolbox_t, std::shared_ptr<toolbox_t>, ModelNumerical>( m, pyclass_name.c_str() )
        .def( py::init([](std::string const& prefix, std::string const& keyword, py::object worldComm, std::string const& subprefix, ModelBaseRepository const& modelRep, int order, int orderGeometry) {
                  worldcomm_ptr_t wc = worldComm.is_none() ? Environment::worldCommPtr() : worldComm.cast<worldcomm_ptr_t>();
                  return new toolbox_t(prefix, keyword, wc, subprefix, modelRep, ModelBaseCommandLineOptions{}, RuntimeOrder{order}, RuntimeOrder{orderGeometry});
              }),
              py::arg( "prefix" ),
              py::arg( "keyword" ) = std::string( "heat" ),
              py::arg( "worldComm" ) = py::none(),
              py::arg( "subprefix" ) = std::string( "" ),
              py::arg( "modelRep" ) = ModelBaseRepository(),
              py::arg( "order" ) = 1,
              py::arg( "orderGeometry" ) = 1,
              "Initialize the heat mechanics toolbox" )
        .def( "init", &toolbox_t::init, "initialize the heat  toolbox", py::arg( "buildModelAlgebraicFactory" ) = true )
        .def( "stabilizationGLS", &toolbox_t::stabilizationGLS, "whether GLS stabilization is enabled" )

        // mesh
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        .def( "rangeMeshElements", &toolbox_t::rangeMeshElements, "get the range of mesh elements" )

        // elements
        .def( "spaceTemperature", &toolbox_t::spaceTemperature, "get the temperature function space" )
        .def( "fieldTemperature", static_cast<element_temperature_t const& (toolbox_t::*)() const>( &toolbox_t::fieldTemperature ), "returns the temperature field" )
        .def( "fieldTemperaturePtr", static_cast<element_temperature_ptr_t const& (toolbox_t::*)() const>( &toolbox_t::fieldTemperaturePtr ), "returns the temperature field shared_ptr" )
        //.def( "spaceVelocityConvection", &toolbox_t::spaceVelocityConvection, "get the field function space")
        // .def( "fieldVelocityConvection", static_cast<element_velocityconvection_t const& (toolbox_t::*)() const>(&toolbox_t::fieldVelocityConvection), "returns the convection velocity field" )
        // .def( "fieldVelocityConvectionPtr", static_cast<element_velocityconvection_ptr_t const& (toolbox_t::*)() const>(&toolbox_t::fieldVelocityConvectionPtr), "returns the convection velocity field shared_ptr" )

        // time stepping
        .def( "timeStepBase", static_cast<std::shared_ptr<TSBase> ( toolbox_t::* )() const>( &toolbox_t::timeStepBase ), "get time stepping base" )
        .def( "startTimeStep", &toolbox_t::startTimeStep, "start time stepping" )
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
                return mat; },
            "returns the assembled matrix" )
        .def(
            "assembleRhs", []( const toolbox_t& t )
            {
                auto rhs = t.algebraicFactory()->rhs()->clone();
                rhs->zero();
                auto matTMP = t.algebraicFactory()->matrix()->clone();
                t.algebraicFactory()->applyAssemblyLinear( t.algebraicBlockVectorSolution()->vectorMonolithic(), matTMP, rhs, { "ignore-assembly.lhs" } );
                return rhs; },
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

    defSM<2>(m);
    defSM<3>(m);

}

