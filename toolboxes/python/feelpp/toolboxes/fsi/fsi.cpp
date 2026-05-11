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
#include <feel/feelpython/pybind11/stl.h>
#include <feel/feelpython/pybind11/functional.h>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/fsi/fsi.hpp>
#include <feel/feelcore/pybind11_json.hpp>
#include "magneto.hpp"



namespace py = pybind11;
using namespace Feel;

template<int nDim, int OrderV, int OrderP,int OrderG>
void defToolbox(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;

    using model_solid_type = FeelModels::SolidMechanics< Simplex<nDim,1>, 
                                                        Lagrange<OrderG, Vectorial,Continuous,PointSetFekete> >;
    using model_fluid_type = FeelModels::FluidMechanics< Simplex<nDim,1>,
                                                        Lagrange<OrderV, Vectorial,Continuous,PointSetFekete>,
                                                        Lagrange<OrderP, Scalar,Continuous,PointSetFekete> >;

    using toolbox_t = FeelModels::FSI<model_fluid_type, model_solid_type>;
    
    using toolbox_ptr_t = std::shared_ptr<toolbox_t>;

    //using space_displacement_t = typename toolbox_t::model_solid_type::space_displacement_type;
    //using element_displacement_t = typename toolbox_t::model_solid_type::element_displacement_type;
    //using element_displacement_ptr_t = typename toolbox_t::model_solid_type::element_displacement_ptrtype;
//    using element_solidfluidpotential_ptr_t = typename toolbox_t::electric_model_type::element_electricpotential_ptrtype;

    std::string pyclass_name = fmt::format("Fsi_{}DP{}P{}G{}",nDim,OrderV,OrderP,OrderG);

    auto add_magneto_torque = [](toolbox_t& FSImodel)
    {
        auto add_torque = [&FSImodel](FeelModels::ModelAlgebraic::DataUpdateLinear & data)
        {
            auto const& t = unwrap_ptr(FSImodel.fluidModel());
            magnetoTorqueModelFSI<nDim,0,toolbox_t>(t, data);
        };

        FSImodel.fluidModel()->algebraicFactory()->addFunctionLinearAssembly( add_torque );
    };

    auto add_magneto_torque_residual = [](toolbox_t& FSImodel)
    {
        auto add_torque_residual = [&FSImodel](FeelModels::ModelAlgebraic::DataUpdateResidual & data)
        {
            auto const& t = unwrap_ptr(FSImodel.fluidModel());
            magnetoTorqueModelFSI<nDim,1, toolbox_t>(t, data);
        };

        FSImodel.fluidModel()->algebraicFactory()->addFunctionResidualAssembly( add_torque_residual );
    };

    py::class_<toolbox_t,std::shared_ptr<toolbox_t>,ModelNumerical>(m,pyclass_name.c_str())
        .def(py::init([](std::string const& prefix, std::string const& keyword, py::object worldComm, std::string const& subprefix, ModelBaseRepository const& modelRep) {
                 worldcomm_ptr_t wc = worldComm.is_none() ? Environment::worldCommPtr() : py::cast<worldcomm_ptr_t>(worldComm);
                 return new toolbox_t(prefix, keyword, wc, subprefix, modelRep);
             }),
             py::arg("prefix"),
             py::arg("keyword")=std::string("fsi"),
             py::arg("worldComm")=py::none(),
             py::arg("subprefix")=std::string(""),
             py::arg("modelRep") = ModelBaseRepository(),
             "Initialize the FSI toolbox"
             )
        .def("init", []( toolbox_t& t, bool /*buildModelAlgebraicFactory*/ ) { t.init(); },
             "initialize the FSI toolbox",
             py::arg("buildModelAlgebraicFactory")=true)

        // mesh
        //.def( "mesh", &toolbox_t::mesh, "get the mesh" ) //TODO
        //.def( "setMesh", &toolbox_t::setMesh, "set the mesh", py::arg( "mesh" ) ) //TODO
        .def( "updateParameterValues", &toolbox_t::updateParameterValues, "update parameter values" )
        .def( "setParameterValues", &toolbox_t::setParameterValues, "set parameter values", py::arg( "paramValues" ) )
        .def( "meshSize", &toolbox_t::meshSize, "get the FSI mesh size" )
        //.def( "rangeMeshElements", &toolbox_t::rangeMeshElements, "get the range of mesh elements" )

        // FSI coupling state
        .def( "fsiCouplingType", &toolbox_t::fsiCouplingType, "get the FSI coupling type" )
        .def( "fsiCouplingBoundaryCondition", &toolbox_t::fsiCouplingBoundaryCondition, "get the FSI coupling boundary condition" )
        .def( "useFSISemiImplicitScheme", &toolbox_t::useFSISemiImplicitScheme, "return true if FSI uses a semi-implicit scheme" )
        .def( "interfaceFSIisConforme", &toolbox_t::interfaceFSIisConforme, "return true if the FSI interface is conforming" )
        .def( "fixPointTolerance", &toolbox_t::fixPointTolerance, "get the FSI fix-point tolerance" )
        .def( "fixPointInitialTheta", &toolbox_t::fixPointInitialTheta, "get the FSI fix-point initial theta" )
        .def( "fixPointMinTheta", &toolbox_t::fixPointMinTheta, "get the FSI fix-point minimum theta" )
        .def( "fixPointMaxIt", &toolbox_t::fixPointMaxIt, "get the FSI fix-point maximum iterations" )
        .def( "fixPointMinItConvergence", &toolbox_t::fixPointMinItConvergence, "get the FSI fix-point minimum convergence iterations" )

        // temperature space and field
        .def( "modelSolid", []( toolbox_ptr_t& t ) { return t->solidModel(); } , "get the solid model" )
        //.def( "spaceDisplacement", []( toolbox_ptr_t& t ) { return t->solidModel()->spaceDisplacement(); } , "get the Displacement function space") //TODO
        .def( "fieldDisplacement", []( toolbox_ptr_t& t ) { return t->solidModel()->fieldDisplacement(); } , "get the Displacement function space")
        .def( "fieldDisplacementPtr", []( toolbox_ptr_t& t ) { return t->solidModel()->fieldDisplacementPtr(); }, "returns the Displacement field shared_ptr" )

        // fluid space and fields
        .def( "modelFluid", []( toolbox_ptr_t& t ) { return t->fluidModel(); } , "get the fluid model" )
        .def( "spaceVelocity", []( toolbox_ptr_t& t ) { return t->fluidModel()->functionSpaceVelocity(); } , "get the velocity function space" )
        .def( "spacePressure", []( toolbox_ptr_t& t ) { return t->fluidModel()->functionSpacePressure(); } , "get the pressure function space" )
        .def( "fieldVelocity", []( toolbox_ptr_t& t ) { return t->fluidModel()->fieldVelocity(); } , "get the velocity field" )
        .def( "fieldPressure", []( toolbox_ptr_t& t ) { return t->fluidModel()->fieldPressure(); } , "get the pressure field" )

        // solve
        .def("solve",&toolbox_t::solve, "solve the FSI problem")
        .def("exportResults",static_cast<void (toolbox_t::*)()>(&toolbox_t::exportResults), "export the results of the FSI problem")
        .def("exportResults",static_cast<void (toolbox_t::*)( double )>(&toolbox_t::exportResults), "export the results of the FSI problem", py::arg("time"))

        //time
        .def("timeStepBase",static_cast<std::shared_ptr<TSBase> (toolbox_t::*)() const>(&toolbox_t::timeStepBase), "get time stepping base")
        .def("fluidTimeStepBase",static_cast<std::shared_ptr<TSBase> (toolbox_t::*)() const>(&toolbox_t::fluidTimeStepBase), "get fluid time stepping base")
        .def("solidTimeStepBase",static_cast<std::shared_ptr<TSBase> (toolbox_t::*)() const>(&toolbox_t::solidTimeStepBase), "get solid time stepping base")
        .def("updateTime",static_cast<void (toolbox_t::*)( double )>(&toolbox_t::updateTime), "update FSI, fluid, and solid model time", py::arg("time"))
        //.def("startTimeStep",static_cast<void (toolbox_t::*)( bool )>(&toolbox_t::startTimeStep), "start time stepping", py::arg("preprocess")=true )
        .def("startTimeStep", &toolbox_t::startTimeStep, "start time stepping")
        .def("updateTimeStep",&toolbox_t::updateTimeStep, "update time stepping")

        .def( "addMagnetoTorqueModelFSI", add_magneto_torque, "add function linear assembly" )
        .def( "addMagnetoTroqueResModelFSI", add_magneto_torque_residual, "add function residual assembly" )
        .def( "addMagnetoTorqueResModelFSI", add_magneto_torque_residual, "add function residual assembly" );
}
    

PYBIND11_MODULE(_fsi, m )
{
    using namespace Feel;

    defToolbox<2,2,1,1>(m);
    defToolbox<2,3,2,1>(m);
    defToolbox<3,2,1,1>(m);
    defToolbox<3,3,2,1>(m);

}
