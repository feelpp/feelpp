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

#ifndef FEELPP_THERMOELECTRIC_HPP
#define FEELPP_THERMOELECTRIC_HPP 1

#include <feel/feelmodels/thermodyn/thermodynamics.hpp>
#include <feel/feelmodels/thermoelectric/electricpropertiesdescription.hpp>


namespace Feel
{
namespace FeelModels
{

enum class ThermoElectricPostProcessFieldExported
{
    Temperature = 0, ElectricPotential, ElectricField, Pid
};

template< typename ConvexType, typename BasisTemperatureType>
class ThermoElectric : public ModelNumerical,
                       public MarkerManagementDirichletBC,
                       public MarkerManagementNeumannBC,
                       public MarkerManagementRobinBC,
                       public boost::enable_shared_from_this< ThermoElectric<ConvexType,BasisTemperatureType> >
{

public:
    typedef ModelNumerical super_type;
    typedef ThermoElectric<ConvexType,BasisTemperatureType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef ThermoDynamics<ConvexType,BasisTemperatureType> thermodyn_model_type;
    typedef boost::shared_ptr<thermodyn_model_type> thermodyn_model_ptrtype;

    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef mesh_type mesh_electric_type;

    // function space electric-potential
    typedef BasisTemperatureType basis_electricpotential_type;
    static const uint16_type nOrderPolyElectricPotential = basis_electricpotential_type::nOrder;
    typedef FunctionSpace<mesh_type, bases<basis_electricpotential_type> > space_electricpotential_type;
    typedef boost::shared_ptr<space_electricpotential_type> space_electricpotential_ptrtype;
    typedef typename space_electricpotential_type::element_type element_electricpotential_type;
    typedef boost::shared_ptr<element_electricpotential_type> element_electricpotential_ptrtype;
    typedef typename space_electricpotential_type::element_external_storage_type element_electricpotential_external_storage_type;
    // function space electric-field
    typedef Lagrange<nOrderPolyElectricPotential, Vectorial,Discontinuous/*Continuous*/,PointSetFekete> basis_electricfield_type;
    typedef FunctionSpace<mesh_electric_type, bases<basis_electricfield_type> > space_electricfield_type;
    typedef boost::shared_ptr<space_electricfield_type> space_electricfield_ptrtype;
    typedef typename space_electricfield_type::element_type element_electricfield_type;
    typedef boost::shared_ptr<element_electricfield_type> element_electricfield_ptrtype;

    // mechanical properties desc
    typedef bases<Lagrange<0, Scalar,Discontinuous> > basis_scalar_P0_type;
    typedef FunctionSpace<mesh_type, basis_scalar_P0_type> space_scalar_P0_type;
    typedef boost::shared_ptr<space_scalar_P0_type> space_scalar_P0_ptrtype;
    typedef ElectricPropertiesDescription<space_scalar_P0_type> electricproperties_type;
    typedef boost::shared_ptr<electricproperties_type> electricproperties_ptrtype;

    // exporter
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef boost::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    // context for evaluation
    typedef typename space_electricpotential_type::Context context_electricpotential_type;
    typedef boost::shared_ptr<context_electricpotential_type> context_electricpotential_ptrtype;


    //___________________________________________________________________________________//
    // constructor
    ThermoElectric( std::string const& prefix,
                    bool buildMesh = true,
                    WorldComm const& _worldComm = Environment::worldComm(),
                    std::string const& subPrefix = "",
                    std::string const& appliShortRepository = ModelBase::rootRepositoryByDefault() );
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"ThermoDynamicsMesh.path"); }
    boost::shared_ptr<std::ostringstream> getInfo() const;

    // load config
    void loadConfigBCFile();
    void loadConfigMeshFile( std::string const& geofilename );
    void loadParameterFromOptionsVm();

    void createMesh();
    // update for use
    void init( bool buildModelAlgebraicFactory = true );
    BlocksBaseGraphCSR buildBlockMatrixGraph() const;
    int nBlockMatrixGraph() const;

    void initPostProcess();
    void restartExporters();
    // void exportMeasures( double time );
    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );
    void setDoExportResults( bool b ) { if (M_exporter) M_exporter->setDoExport( b ); }
    bool hasPostProcessFieldExported( ThermoElectricPostProcessFieldExported const& key ) const { return M_postProcessFieldExported.find( key ) != M_postProcessFieldExported.end(); }

    void updateParameterValues();

    //___________________________________________________________________________________//

    mesh_ptrtype const& mesh() const { return M_mesh; }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

    thermodyn_model_ptrtype const& thermodynModel() const { return M_thermodynModel; }
    thermodyn_model_ptrtype thermodynModel() { return M_thermodynModel; }

    space_electricpotential_ptrtype const& spaceElectricPotential() const { return M_XhElectricPotential; }
    element_electricpotential_ptrtype const& fieldElectricPotentialPtr() const { return M_fieldElectricPotential; }
    element_electricpotential_type const& fieldElectricPotential() const { return *M_fieldElectricPotential; }

    space_electricfield_ptrtype const& spaceElectricField() const { return M_XhElectricField; }
    element_electricfield_ptrtype const& fieldElectricFieldPtr() const { return M_fieldElectricField; }
    element_electricfield_type const& fieldElectricField() const { return *M_fieldElectricField; }

    electricproperties_ptrtype const& electricProperties() const { return M_electricProperties; }

    backend_ptrtype const& backend() const { return M_backendMonolithic; }
    BlocksBaseVector<double> const& blockVectorSolutionMonolithic() const { return M_blockVectorSolutionMonolithic; }
    BlocksBaseVector<double> & blockVectorSolutionMonolithic() { return M_blockVectorSolutionMonolithic; }

    //___________________________________________________________________________________//

    boost::shared_ptr<TSBase> timeStepBase() { return this->thermodynModel()->timeStepBase(); }
    boost::shared_ptr<TSBase> timeStepBase() const { return this->thermodynModel()->timeStepBase(); }
    void updateTimeStep() {  this->thermodynModel()->updateTimeStep(); }

    //___________________________________________________________________________________//
    // apply assembly and solver
    void solve();

    void updateLinearPDE( DataUpdateLinear & data ) const;
    void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const;
    void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;
    void updateSourceTermLinearPDE(vector_ptrtype& F, bool buildCstPart) const;
    void updateLinearPreAssemblyJouleLaw( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const;

    void updateNewtonInitialGuess( vector_ptrtype& U ) const;
    void updateJacobian( DataUpdateJacobian & data ) const;
    void updateBCStrongDirichletJacobian(sparse_matrix_ptrtype& J,vector_ptrtype& RBis ) const;
    void updateBCWeakJacobian( element_electricpotential_external_storage_type const& v, sparse_matrix_ptrtype& J, bool buildCstPart ) const;
    void updateResidual( DataUpdateResidual & data ) const;
    void updateBCDirichletStrongResidual( vector_ptrtype& R ) const;
    void updateBCWeakResidual( element_electricpotential_external_storage_type const& v, vector_ptrtype& R, bool buildCstPart ) const;
    void updateSourceTermResidual( vector_ptrtype& R, bool buildCstPart ) const;

    //___________________________________________________________________________________//
    void updateElectricField();

private :
    void updateBoundaryConditionsForUse();

private :
    thermodyn_model_ptrtype M_thermodynModel;

    bool M_hasBuildFromMesh, M_isUpdatedForUse;

    mesh_ptrtype M_mesh;
    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

    space_electricpotential_ptrtype M_XhElectricPotential;
    element_electricpotential_ptrtype M_fieldElectricPotential;
    space_electricfield_ptrtype M_XhElectricField;
    element_electricfield_ptrtype M_fieldElectricField;
    // physical parameter
    std::string M_modelName;
    std::string M_solverName;
    electricproperties_ptrtype M_electricProperties;
    // boundary conditions
    map_scalar_field<2> M_bcDirichlet;
    map_scalar_field<2> M_bcNeumann;
    map_scalar_fields<2> M_bcRobin;
    map_scalar_field<2> M_volumicForcesProperties;

    // algebraic data/tools
    backend_ptrtype M_backendMonolithic;
    model_algebraic_factory_ptrtype M_algebraicFactoryMonolithic;
    BlocksBaseVector<double> M_blockVectorSolutionMonolithic;
    std::map<std::string,std::set<size_type> > M_dofsWithValueImposed;
    // start dof index fields in matrix (temperature,electric-potential,...)
    std::map<std::string,size_type> M_startBlockIndexFieldsInMatrix;

    backend_ptrtype M_backendElectricModel;
    model_algebraic_factory_ptrtype M_algebraicFactoryElectricModel;

    // post-process
    export_ptrtype M_exporter;
    std::set<ThermoElectricPostProcessFieldExported> M_postProcessFieldExported;
    std::set<std::string> M_postProcessUserFieldExported;


};

} // namespace FeelModels
} // namespace Feel

#endif // INCLUDE_THERMOELECTRIC_HPP
