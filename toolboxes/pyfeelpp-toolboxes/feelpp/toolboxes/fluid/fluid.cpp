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
#include <pybind11/stl.h>
#include <pybind11/functional.h>
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
        .def("mesh",static_cast<typename fm_t::mesh_ptrtype  (fm_t::*)() const>(&fm_t::mesh), "get the mesh")
        .def("setMesh",static_cast<void (fm_t::*)(typename fm_t::mesh_ptrtype const&)>(&fm_t::setMesh), "set the mesh of the toolbox",py::arg("mesh"))
        // function spaces and elements
        .def("functionSpaceVelocity",&fm_t::functionSpaceVelocity, "get the velocity function space")
        //.def("fieldVelocity",static_cast<typename fm_t::element_velocity_ptrtype& (fm_t::*)()>(&fm_t::fieldVelocityPtr), "get the velocity field")
        .def("fieldVelocity",[]( std::shared_ptr<fm_t>& self ) {
            self->fieldVelocityPtr()->printMatlab("velocityptr.m");
            return self->fieldVelocityPtr();
        } )
        .def("setFieldVelocity",
            []( std::shared_ptr<fm_t>& self, typename fm_t::element_velocity_ptrtype& v ) {
                v->printMatlab("v.m");
                self->fieldVelocity() = *v;
                self->fieldVelocityPtr()->printMatlab("velocity.m");
            }, "set the velocity field", py::arg("field"))
        .def("setFieldPressure",
            []( std::shared_ptr<fm_t>& self, typename fm_t::element_pressure_ptrtype& p ) {
                self->fieldPressure() = *p;
            }, "set the pressure field", py::arg("field"))
        .def("functionSpacePressure",&fm_t::functionSpacePressure, "get the pressure function space")
        .def("fieldPressure",static_cast<typename fm_t::element_pressure_ptrtype const& (fm_t::*)() const>(&fm_t::fieldPressurePtr), "get the pressure field")

        // time stepping
        .def("timeStepBase",static_cast<std::shared_ptr<TSBase> (fm_t::*)() const>(&fm_t::timeStepBase), "get time stepping base")
        .def("startTimeStep",static_cast<void (fm_t::*)( bool )>(&fm_t::startTimeStep), "start time stepping", py::arg("preprocess")=true )
        .def("updateTimeStep",&fm_t::updateTimeStep, "update time stepping")

        // solve
        .def("solve",&fm_t::solve, "solve the fluid mechanics problem")
        .def("exportResults",static_cast<void (fm_t::*)()>(&fm_t::exportResults), "export the results of the fluid mechanics problem")
        .def("exportResults",static_cast<void (fm_t::*)( double )>(&fm_t::exportResults), "export the results of the fluid mechanics problem", py::arg("time"))

        // remesh
        .def("applyRemesh",&fm_t::applyRemesh, "apply remesh to toolbox and regenerate the necessary data structure")

        .def(
            "contactForce",[]( const fm_t& t, FeelModels::ModelAlgebraic::DataUpdateLinear & data)
            {
                //auto data = t.algebraicFactory();
                
                //bool buildCstPart = data->buildCstPart();
                //if( buildCstPart )
                //    return;

                double r = 0.125;
                double rho = 0.12;
                double eps = 0.005;
                double c = 0.001; 

                // Get center of mass
                std::vector<double> massCenters;
                auto nbr_bodies = 0;

                for ( auto const& [bpname,bpbc] : t.bodySetBC() )
                {
                    for (auto &coord_massCenter : bpbc.body().massCenter())
                    {
                        massCenters.push_back(coord_massCenter);
                    }   
                    nbr_bodies ++;
                }

                for (auto it = massCenters.begin(); it != massCenters.end(); ++it)
                {
                    std::cout << "Coordinates : " << *it << std::endl;
                }


                // Define repulsive force
                std::vector<std::vector<double> > rep_forces(nDim, std::vector<double>(nbr_bodies));

                for (int b_i = 0;b_i < nbr_bodies; b_i++)
                {
                    for (int b_j=b_i+1; b_j < nbr_bodies; b_j++)
                    {
                        auto dist_ij = sqrt(std::pow(massCenters[2*b_j]- massCenters[2*b_i],2) + std::pow(massCenters[2*b_j + 1]-massCenters[2*b_i + 1],2));
                        std::cout << "Distance : " << dist_ij << std::endl;
                        
                        if (-(dist_ij-r-r-rho)/rho > 0)
                        {
                            for (int d = 0; d<nDim; d++)
                            {
                                auto G_ij_d = (massCenters[2*b_j+d]- massCenters[2*b_i+d])/dist_ij;
                                auto F_ij_d = c/eps * std::pow(-(dist_ij-r-r-rho)/rho,2)*G_ij_d;
                                rep_forces[d][b_i] += F_ij_d;
                                rep_forces[d][b_j] += - F_ij_d;
                            }
                        }
                    }
                }

                // Affichage
                for(auto& row: rep_forces)
                {
                    for(auto& col : row)
                    { 
                        std::cout << col << std::endl;
                    }
                }
                            
                // Ajout de la force
                auto rhs = data.rhs();
                auto rowStartInVector = t.rowStartInVector();

                int B = 0;
                for ( auto const& [bpname,bpbc] : t.bodySetBC() )
                {
                    size_type startBlockIndexTranslationalVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
                    rhs->setIsClosed( false );

                    auto const& basisToContainerGpTranslationalVelocityVector = rhs->map().dofIdToContainerId( rowStartInVector+startBlockIndexTranslationalVelocity );
                    
                    for (int d=0;d<nDim;++d)
                    {   
                        std::cout << "Force ajoutée : " << rep_forces[d][B] << " dim : " << d << " body : " << B << std::endl;
                        
                        rhs->add( basisToContainerGpTranslationalVelocityVector[d],
                                 rep_forces[d][B]);       

                    }
                    B++;   
                }
            },
            "contact force"
        )
  
        .def(
            "contactForceRes",[](const fm_t& t,FeelModels::ModelAlgebraic::DataUpdateResidual & data)
            {
                //auto data = t.algebraicFactory();

                //bool buildCstPart = data->buildCstPart();
                //if( buildCstPart )
                //    return;

                double r = 0.125;
                double rho = 0.12;
                double eps = 0.00005;
                double c = 0.001; 

                // Get center of mass
                std::vector<double> massCenters;
                auto nbr_bodies = 0;

                for ( auto const& [bpname,bpbc] : t.bodySetBC() )
                {

                    for (auto &coord_massCenter : bpbc.body().massCenter())
                    {
                        massCenters.push_back(coord_massCenter);
                    }   
                    nbr_bodies ++;
                }

                
                for (auto it = massCenters.begin(); it != massCenters.end(); ++it)
                {
                    std::cout << "Coordinates : " << *it << std::endl;
                }


                // Define repulsive force
                std::vector<std::vector<double> > rep_forces(nDim, std::vector<double>(nbr_bodies));

                for (int b_i = 0;b_i < nbr_bodies; b_i++)
                {
                    for (int b_j=b_i+1; b_j < nbr_bodies; b_j++)
                    {
                        auto dist_ij = sqrt(std::pow(massCenters[2*b_j]- massCenters[2*b_i],2) + std::pow(massCenters[2*b_j + 1]-massCenters[2*b_i + 1],2));
                        std::cout << "Distance : " << dist_ij << std::endl;


                        std::cout << "Distance - rho" << - (dist_ij-r-r-rho) << std::endl;
                        
                        if (-(dist_ij-r-r-rho)/rho > 0)
                        {
                            for (int d = 0; d<nDim; d++)
                            {
                                auto G_ij_d = (massCenters[2*b_j+d]- massCenters[2*b_i+d])/dist_ij;
                                auto F_ij_d = c/eps * std::pow(-(dist_ij-r-r-rho)/rho,2)*G_ij_d;
                                rep_forces[d][b_i] += F_ij_d;
                                rep_forces[d][b_j] += - F_ij_d;
                            }
                        }                       
                    }
                }

                // Affichage
                for(auto& row: rep_forces)
                {
                    for(auto& col : row)
                    { 
                        std::cout << col << std::endl;
                    }
                }
                            
                // Ajout de la force
                //auto rhs = data->residual();
                auto rhs = data.residual();
                auto rowStartInVector = t.rowStartInVector();

                int B = 0;
                for ( auto const& [bpname,bpbc] : t.bodySetBC() )
                {
                    size_type startBlockIndexTranslationalVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
                    
                    rhs->setIsClosed( false );

                    auto const& basisToContainerGpTranslationalVelocityVector = rhs->map().dofIdToContainerId( rowStartInVector+startBlockIndexTranslationalVelocity );
                    
                    for (int d=0;d<nDim;++d)
                    {   
                        std::cout << "Force ajoutée : " << rep_forces[d][B] << " dim : " << d << " body : " << B << std::endl;
                        
                        rhs->add( basisToContainerGpTranslationalVelocityVector[d],
                                 rep_forces[d][B]);       

                    }
                    B++;   
                }
            },
            "contact res force"
        )

        .def(
            "linAssemblyPython",[](const fm_t& t, std::function<void (FeelModels::ModelAlgebraic::DataUpdateLinear& )> const& func)
            {
                t.algebraicFactory()->addFunctionLinearAssembly(func);

            },
            "add function linear assembly"
        )

        .def(
            "ResAssemblyPython",[](const fm_t& t, std::function<void(FeelModels::ModelAlgebraic::DataUpdateResidual& )> const& func)
            {
                t.algebraicFactory()->addFunctionResidualAssembly(func);

            },
            "add function residual assembly"
        )
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

