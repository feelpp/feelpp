/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2016-12-12

  Copyright (C) 2016 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file thermoelectric.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2016-12-12
 */

#ifndef FEELPP_TOOLBOXES_THERMOELECTRIC_HPP
#define FEELPP_TOOLBOXES_THERMOELECTRIC_HPP 1

#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feelmodels/electric/electric.hpp>


namespace Feel
{
namespace FeelModels
{

template< typename HeatType, typename ElectricType>
class ThermoElectric : public ModelNumerical,
                       public std::enable_shared_from_this< ThermoElectric<HeatType,ElectricType> >
{

public:
    typedef ModelNumerical super_type;
    typedef ThermoElectric<HeatType,ElectricType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef HeatType heat_model_type;
    typedef std::shared_ptr<heat_model_type> heat_model_ptrtype;

    typedef ElectricType electric_model_type;
    typedef std::shared_ptr<electric_model_type> electric_model_ptrtype;

    // mesh
    typedef typename heat_model_type::mesh_type mesh_heat_type;
    typedef typename electric_model_type::mesh_type mesh_electric_type;
    typedef mesh_heat_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    // exporter
    typedef Exporter<mesh_type,mesh_type::nOrder> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    //___________________________________________________________________________________//
    // constructor
    ThermoElectric( std::string const& prefix,
                    bool buildMesh = true,
                    worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
                    std::string const& subPrefix = "",
                    ModelBaseRepository const& modelRep = ModelBaseRepository() );
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"ThermoElectricMesh.path"); }
    std::shared_ptr<std::ostringstream> getInfo() const;

private :
    void loadParameterFromOptionsVm();
    void createMesh();
    void initPostProcess();
public :
    // update for use
    void init( bool buildModelAlgebraicFactory = true );

    BlocksBaseGraphCSR buildBlockMatrixGraph() const;
    int nBlockMatrixGraph() const;

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );

    void updateParameterValues();

    //___________________________________________________________________________________//

    mesh_ptrtype const& mesh() const { return M_mesh; }
    void setMesh( mesh_ptrtype const& mesh ) { M_mesh = mesh; }
    //elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

    heat_model_ptrtype const& heatModel() const { return M_heatModel; }
    heat_model_ptrtype heatModel() { return M_heatModel; }

    electric_model_ptrtype const& electricModel() const { return M_electricModel; }
    electric_model_ptrtype electricModel() { return M_electricModel; }


    backend_ptrtype const& backend() const { return M_backendMonolithic; }
    BlocksBaseVector<double> const& blockVectorSolutionMonolithic() const { return M_blockVectorSolutionMonolithic; }
    BlocksBaseVector<double> & blockVectorSolutionMonolithic() { return M_blockVectorSolutionMonolithic; }

    //___________________________________________________________________________________//

    std::shared_ptr<TSBase> timeStepBase() { return this->heatModel()->timeStepBase(); }
    std::shared_ptr<TSBase> timeStepBase() const { return this->heatModel()->timeStepBase(); }
    void updateTimeStep() {  this->heatModel()->updateTimeStep(); }

    //___________________________________________________________________________________//
    // apply assembly and solver
    void solve();

    void updateLinearPreAssemblyJouleLaw( DataUpdateLinear & data ) const;
    void updateResidualPreAssemblyJouleLaw( DataUpdateResidual & data ) const;
    void updateGenericPreAssemblyJouleLaw( vector_ptrtype& F, bool applyOnResidual ) const;

    void updateLinearElectricDependingOnTemperature( DataUpdateLinear & data ) const;

    void updateLinearPDE( DataUpdateLinear & data ) const;

    void updateNewtonInitialGuess( vector_ptrtype& U ) const;
    void updateJacobian( DataUpdateJacobian & data ) const;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const;
    void updateResidual( DataUpdateResidual & data ) const;
    void updateResidualDofElimination( DataUpdateResidual & data ) const;

    //___________________________________________________________________________________//
    void updateCurrentDensity();

private :
    heat_model_ptrtype M_heatModel;
    electric_model_ptrtype M_electricModel;

    bool M_hasBuildFromMesh, M_isUpdatedForUse;

    mesh_ptrtype M_mesh;
    //elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;
    // materials range
    std::map<std::string, elements_reference_wrapper_t<mesh_type> > M_rangeMeshElementsByMaterial;

    // physical parameter
    std::string M_modelName;
    bool M_modelUseJouleEffect;

    // solver
    std::string M_solverName;
    bool M_solverNewtonInitialGuessUseLinearThermoElectric,M_solverNewtonInitialGuessUseLinearHeat,M_solverNewtonInitialGuessUseLinearElectric;

    // algebraic data/tools
    backend_ptrtype M_backendMonolithic;
    model_algebraic_factory_ptrtype M_algebraicFactoryMonolithic;
    BlocksBaseVector<double> M_blockVectorSolutionMonolithic;

    // post-process
    export_ptrtype M_exporter;
    std::set<std::string> M_postProcessFieldExportedHeat, M_postProcessFieldExportedElectric;
};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_TOOLBOXES_THERMOELECTRIC_HPP
