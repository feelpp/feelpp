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
#include <fmt/core.h>
#include <feel/feelmodels/coefficientformpdes/coefficientformpdes.hpp>
#include <feel/feelmodels/coefficientformpdes/coefficientformpdes_registered_type.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelpython/pybind11/pybind11.h>

namespace py = pybind11;
using namespace Feel;

template<int nDim, typename BasisT>
void defcfpde(py::module &m, std::string const& basisname )
{
    using namespace Feel;
    using namespace Feel::FeelModels;
    using toolbox_cfpde_t = FeelModels::CoefficientFormPDE<Simplex<nDim, 1>,BasisT>;
    using exporter_cfpde_t = Exporter<typename toolbox_cfpde_t::mesh_type,1>;
    using exporter_cfpde_ptr_t = std::shared_ptr<exporter_cfpde_t>;
    using step_t = typename exporter_cfpde_t::step_type;
    using step_ptr_t = typename exporter_cfpde_t::step_ptrtype;

    std::string pyclass_name = fmt::format("cfpde_{}D_{}",nDim,basisname);
    py::class_<toolbox_cfpde_t, std::shared_ptr<toolbox_cfpde_t>, ModelNumerical>( m, pyclass_name.c_str() )
#if 0    
        .def( py::init<std::string const&, std::string const&, worldcomm_ptr_t const&, std::string const&, ModelBaseRepository const&>(),
              py::arg( "prefix" ),
              py::arg( "keyword" ) = std::string( "cfpde" ),
              py::arg( "worldComm" ) = Environment::worldCommPtr(),
              py::arg( "subprefix" ) = std::string( "" ),
              py::arg( "modelRep" ) = ModelBaseRepository(),
              "Initialize the coefficient form pdes toolbox" )

        .def( "init", &toolbox_cfpde_t::init, "initialize the heat  toolbox", py::arg( "buildModelAlgebraicFactory" ) = true )
#endif
        // mesh
        //.def( "mesh", &toolbox_cfpde_t::mesh, "get the mesh", py::arg("mesh") )
        //.def( "setMesh", &toolbox_cfpde_t::setMesh, "set the mesh", py::arg("mesh") )
        .def( "unknownIsScalar", &toolbox_cfpde_t::unknownIsScalar, "return true if the unknown is scalar")
        .def( "spaceUnknown", &toolbox_cfpde_t::spaceUnknown, "return the space of the unknown")
        .def( "fieldUnknown", &toolbox_cfpde_t::fieldUnknownPtr, "return the field of the unknown");

}
template<int nDim>
void defcfpdes(py::module &m)
{
    using namespace Feel;
    using namespace Feel::FeelModels;

    using toolbox_cfpdes_t = FeelModels::coefficient_form_PDEs_t<Simplex<nDim, 1>>;
    using exporter_cfpdes_t = Exporter<typename toolbox_cfpdes_t::mesh_type,1>;
    using exporter_cfpdes_ptr_t = std::shared_ptr<exporter_cfpdes_t>;
    using step_t = typename exporter_cfpdes_t::step_type;
    using step_ptr_t = typename exporter_cfpdes_t::step_ptrtype;

    std::string pyclass_name = std::string("cfpdes_") + std::to_string(nDim) + std::string("D");
    py::class_<toolbox_cfpdes_t, std::shared_ptr<toolbox_cfpdes_t>, ModelNumerical>( m, pyclass_name.c_str() )
        .def( py::init<std::string const&, std::string const&, worldcomm_ptr_t const&, std::string const&, ModelBaseRepository const&>(),
              py::arg( "prefix" ),
              py::arg( "keyword" ) = std::string( "cfpdes" ),
              py::arg( "worldComm" ) = Environment::worldCommPtr(),
              py::arg( "subprefix" ) = std::string( "" ),
              py::arg( "modelRep" ) = ModelBaseRepository(),
              "Initialize the coefficient form pdes toolbox" )
        .def( "init", &toolbox_cfpdes_t::init, "initialize the heat  toolbox", py::arg( "buildModelAlgebraicFactory" ) = true )

        // mesh
        .def( "mesh", &toolbox_cfpdes_t::mesh, "get the mesh" )
        .def( "setMesh", &toolbox_cfpdes_t::setMesh, "set the mesh" )
        .def(
            "exportSolutionToStep", []( std::shared_ptr<toolbox_cfpdes_t> const& t, step_ptr_t& s )
            {
                t->apply( [&s]( auto const& cfpde )
                          {
                              std::cout << fmt::format( "[cfpde] exporting {}_{} ...", cfpde->equationName(), cfpde->unknownName() ) << std::endl;
                              s->add( fmt::format( "{}_{}", cfpde->equationName(), cfpde->unknownName() ), cfpde->fieldUnknown() );
                          } );
            },
            "apply external exporter on solution" )

        //.def( "rangeMeshElements", &toolbox_cfpdes_t::rangeMeshElements, "get the range of mesh elements" )

        // elements
        //.def( "spaceTemperature", &toolbox_cfpdes_t::spaceTemperature, "get the temperature function space")
        //.def( "fieldTemperature", static_cast<element_temperature_t const& (toolbox_cfpdes_t::*)() const>(&toolbox_cfpdes_t::fieldTemperature), "returns the temperature field" )
        //.def( "fieldTemperaturePtr", static_cast<element_temperature_ptr_t const& (toolbox_cfpdes_t::*)() const>(&toolbox_cfpdes_t::fieldTemperaturePtr), "returns the temperature field shared_ptr" )

        // time stepping
        .def( "timeStepBase", static_cast<std::shared_ptr<TSBase> ( toolbox_cfpdes_t::* )() const>( &toolbox_cfpdes_t::timeStepBase ), "get time stepping base" )
        .def( "startTimeStep", &toolbox_cfpdes_t::startTimeStep, "start time stepping" )
        .def( "updateTimeStep", &toolbox_cfpdes_t::updateTimeStep, "update time stepping" )
        // solve
        .def( "solve", &toolbox_cfpdes_t::solve, "solve the cfpde problem" )
        .def( "exportResults", static_cast<void ( toolbox_cfpdes_t::* )()>( &toolbox_cfpdes_t::exportResults ), "export the results of the cfpde problem" )
        .def( "exportResults", static_cast<void ( toolbox_cfpdes_t::* )( double )>( &toolbox_cfpdes_t::exportResults ), py::arg("time"), "export the results of the cfpde problem at time 'time'" )
        .def( "checkResults", static_cast<bool ( toolbox_cfpdes_t::* )() const>( &toolbox_cfpdes_t::checkResults ), "check the results of the cfpde problem" )
        //        .def("exportResults",static_cast<void (toolbox_cfpdes_t::*)( double )>(&toolbox_cfpdes_t::exportResults), "export the results of the heat mechanics problem", py::arg("time"))

        .def( "setMesh", &toolbox_cfpdes_t::setMesh, "set the mesh", py::arg( "mesh" ) )
        .def( "updateParameterValues", &toolbox_cfpdes_t::updateParameterValues, "update parameter values" )
        .def( "pdePch1", 
              []( toolbox_cfpdes_t const& toolbox, std::string& nameeq ){
                using t = hana::type<Lagrange<1, Scalar, Continuous, PointSetFekete>>;
                return toolbox.template coefficientFormPDE<t>(nameeq, t() );                
              }, "get the coefficients form pdes for a Pch1", py::arg( "name" ) )
        .def( "pdePch2", 
              []( toolbox_cfpdes_t const& toolbox, std::string& nameeq ){
                using t = hana::type<Lagrange<2, Scalar, Continuous, PointSetFekete>>;
                return toolbox.template coefficientFormPDE<t>(nameeq, t() ); 
              }, "get the coefficients form pdes for a Pch2", py::arg( "name" ) )
        .def( "pdePchv1", 
              []( toolbox_cfpdes_t const& toolbox, std::string& nameeq ){
                using t = hana::type<Lagrange<1, Vectorial, Continuous, PointSetFekete>>;
                return toolbox.template coefficientFormPDE<t>(nameeq, t() ); 
              }, "get the coefficients form pdes for a Pchv1", py::arg( "name" ) )
        .def( "pdePchv2", 
              []( toolbox_cfpdes_t const& toolbox, std::string& nameeq ){
                using t = hana::type<Lagrange<2, Vectorial, Continuous, PointSetFekete>>;
                return toolbox.template coefficientFormPDE<t>(nameeq, t() );
              }, "get the coefficients form pdes for a Pchv2", py::arg( "name" ) )
        ;
}


PYBIND11_MODULE(_cfpdes, m )
{
    using namespace Feel;
    
    defcfpde<2, Lagrange<1, Scalar, Continuous, PointSetFekete>>( m, "Pch1" );
    defcfpde<3, Lagrange<1, Scalar, Continuous, PointSetFekete>>( m, "Pch1" );
    defcfpde<2, Lagrange<2, Scalar, Continuous, PointSetFekete>>( m, "Pch2" );
    defcfpde<3, Lagrange<2, Scalar, Continuous, PointSetFekete>>( m, "Pch2" );
    defcfpde<2, Lagrange<1, Vectorial, Continuous, PointSetFekete>>( m, "Pchv1" );
    defcfpde<3, Lagrange<1, Vectorial, Continuous, PointSetFekete>>( m, "Pchv1" );
    defcfpde<2, Lagrange<2, Vectorial, Continuous, PointSetFekete>>( m, "Pchv2" );
    defcfpde<3, Lagrange<2, Vectorial, Continuous, PointSetFekete>>( m, "Pchv2" );

    defcfpdes<2>(m);
    defcfpdes<3>(m);

}

