/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2014-06-04

  Copyright (C) 2014 Université Joseph Fourier (Grenoble I)

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
   \file thermodyn.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-06-04
 */

#ifndef FEELPP_THERMODYNAMICSBASE_HPP
#define FEELPP_THERMODYNAMICSBASE_HPP 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
//#include <feel/feelvf/vf.hpp>
#include <feel/feelts/bdf.hpp>

#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/options.hpp>

#include <feel/feelmodels/thermodyn/thermalpropertiesdescription.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename ConvexType, typename BasisTemperatureType>
class ThermoDynamicsBase : public ModelNumerical,
                           public MarkerManagementDirichletBC,
                           public MarkerManagementNeumannBC,
                           public MarkerManagementRobinBC
{
  public:
    typedef ModelNumerical super_type;
    typedef ThermoDynamicsBase<ConvexType, BasisTemperatureType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    // basis
    static const uint16_type nOrderPoly = BasisTemperatureType::nOrder;
    typedef BasisTemperatureType basis_temperature_type;
    typedef Lagrange<nOrderPoly, Vectorial, Continuous, PointSetFekete> basis_velocityconvection_type;
    // function space temperature
    typedef FunctionSpace<mesh_type, bases<basis_temperature_type>> space_temperature_type;
    typedef boost::shared_ptr<space_temperature_type> space_temperature_ptrtype;
    typedef typename space_temperature_type::element_type element_temperature_type;
    typedef boost::shared_ptr<element_temperature_type> element_temperature_ptrtype;
    typedef typename space_temperature_type::element_external_storage_type element_temperature_external_storage_type;
    // function space velocity convection
    typedef FunctionSpace<mesh_type, bases<basis_velocityconvection_type>> space_velocityconvection_type;
    typedef boost::shared_ptr<space_velocityconvection_type> space_velocityconvection_ptrtype;
    typedef typename space_velocityconvection_type::element_type element_velocityconvection_type;
    typedef boost::shared_ptr<element_velocityconvection_type> element_velocityconvection_ptrtype;
    // mechanical properties desc
    typedef bases<Lagrange<0, Scalar, Discontinuous>> basis_scalar_P0_type;
    typedef FunctionSpace<mesh_type, basis_scalar_P0_type> space_scalar_P0_type;
    typedef boost::shared_ptr<space_scalar_P0_type> space_scalar_P0_ptrtype;
    typedef ThermalPropertiesDescription<space_scalar_P0_type> thermalproperties_type;
    typedef boost::shared_ptr<thermalproperties_type> thermalproperties_ptrtype;
    // time scheme
    typedef Bdf<space_temperature_type> bdf_temperature_type;
    typedef boost::shared_ptr<bdf_temperature_type> bdf_temperature_ptrtype;
    // exporter
    typedef Exporter<mesh_type, nOrderGeo> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef boost::shared_ptr<model_algebraic_factory_type> model_algebraic_factory_ptrtype;

    // context for evaluation
    typedef typename space_temperature_type::Context context_temperature_type;
    typedef boost::shared_ptr<context_temperature_type> context_temperature_ptrtype;

    ThermoDynamicsBase( std::string const& prefix,
                        bool buildMesh,
                        WorldComm const& worldComm,
                        std::string const& subPrefix,
                        std::string const& rootRepository );

    std::string fileNameMeshPath() const { return prefixvm( this->prefix(), "ThermoDynamicsMesh.path" ); }
    //___________________________________________________________________________________//
    // mesh, space, element temperature
    mesh_ptrtype const& mesh() const { return M_mesh; }
    space_temperature_ptrtype const& spaceTemperature() const { return M_Xh; }
    element_temperature_ptrtype const& fieldTemperaturePtr() const { return M_fieldTemperature; }
    element_temperature_type const& fieldTemperature() const { return *M_fieldTemperature; }
    element_velocityconvection_ptrtype const& fieldVelocityConvectionPtr() const { return M_fieldVelocityConvection; }
    element_velocityconvection_ptrtype& fieldVelocityConvectionPtr() { return M_fieldVelocityConvection; }
    element_velocityconvection_type const& fieldVelocityConvection() const { return *M_fieldVelocityConvection; }
    bool fieldVelocityConvectionIsUsed() const { return M_fieldVelocityConvectionIsUsed; }
    bool fieldVelocityConvectionIsIncompressible() const { return M_fieldVelocityConvectionIsIncompressible; }
    void setFieldVelocityConvectionIsUsed( bool b ) { M_fieldVelocityConvectionIsUsed = b; }
    bool fieldVelocityConvectionIsOperational() const { return ( M_fieldVelocityConvection.use_count() > 0 ); }
    bool fieldVelocityConvectionIsUsedAndOperational() const { return this->fieldVelocityConvectionIsUsed() && this->fieldVelocityConvectionIsOperational(); }
    void setFieldVelocityConvectionIsIncompressible( bool b ) { M_fieldVelocityConvectionIsIncompressible = b; }
    //___________________________________________________________________________________//
    // physical parameters
    thermalproperties_ptrtype const& thermalProperties() const { return M_thermalProperties; }
    thermalproperties_ptrtype& thermalProperties() { return M_thermalProperties; }

    //___________________________________________________________________________________//
    // algebraic data and solver
    backend_ptrtype const& backend() const { return M_backend; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    BlocksBaseVector<double>& blockVectorSolution() { return M_blockVectorSolution; }
    size_type nLocalDof() const;
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype& algebraicFactory() { return M_algebraicFactory; }
    //___________________________________________________________________________________//
    // time step scheme
    bdf_temperature_ptrtype const& timeStepBdfTemperature() const { return M_bdfTemperature; }
    boost::shared_ptr<TSBase> timeStepBase() { return this->timeStepBdfTemperature(); }
    boost::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBdfTemperature(); }
    void updateBdf();
    void updateTimeStep() { this->updateBdf(); }
    //___________________________________________________________________________________//

    boost::shared_ptr<std::ostringstream> getInfo() const;

    virtual void loadConfigBCFile() = 0;
    virtual void loadConfigMeshFile( std::string const& geofilename ) = 0;

    void loadParameterFromOptionsVm();
    void createMesh();
    void createFunctionSpaces();
    void createTimeDiscretisation();
    void createExporters();
    BlocksBaseGraphCSR buildBlockMatrixGraph() const;
    int nBlockMatrixGraph() const { return 1; }
    void init( bool buildModelAlgebraicFactory, model_algebraic_factory_type::appli_ptrtype const& app );
    void updateForUseFunctionSpacesVelocityConvection();

    void initPostProcess();
    void restartExporters();
    void exportMeasures( double time );
    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );
    void setDoExportResults( bool b )
    {
        if ( M_exporter ) M_exporter->setDoExport( b );
    }

    void build();
    void loadMesh( mesh_ptrtype mesh );

    void updateParameterValues();

    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // apply assembly and solver
    /*virtual*/ void solve();

    void updateLinearPDE( DataUpdateLinear& data ) const;
    virtual void updateWeakBCLinearPDE( sparse_matrix_ptrtype& A, vector_ptrtype& F, bool buildCstPart ) const = 0;
    virtual void updateBCStrongDirichletLinearPDE( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const = 0;
    virtual void updateSourceTermLinearPDE( vector_ptrtype& F, bool buildCstPart ) const = 0;

    // non linear (newton)
    void updateNewtonInitialGuess( vector_ptrtype& U ) const;
    void updateJacobian( DataUpdateJacobian& data ) const;
    void updateResidual( DataUpdateResidual& data ) const;
    void updateBCDirichletStrongResidual( vector_ptrtype& R ) const;
    void updateBCNeumannResidual( vector_ptrtype& R, bool buildCstPart ) const;
    void updateBCRobinResidual( element_temperature_external_storage_type const& u, vector_ptrtype& R, bool buildCstPart ) const;
    void updateSourceTermResidual( vector_ptrtype& R, bool buildCstPart ) const;
    void updateBCStrongDirichletJacobian( sparse_matrix_ptrtype& J, vector_ptrtype& RBis ) const;
    void updateBCRobinJacobian( sparse_matrix_ptrtype& J, bool buildCstPart ) const;

    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // update field from expr
    void updateFieldVelocityConvection( bool onlyExprWithTimeSymbol = false );
    template <typename ExprT>
    void updateFieldVelocityConvection( vf::Expr<ExprT> const& expr )
    {
        if ( !M_fieldVelocityConvection )
            this->updateForUseFunctionSpacesVelocityConvection();
        M_exprVelocityConvection.reset(); // symbolic expression is remove
        M_fieldVelocityConvection->on( _range = elements( this->mesh() ), _expr = expr );
    }

    //private :
  protected:
    bool M_hasBuildFromMesh, M_isUpdatedForUse;

    mesh_ptrtype M_mesh;

    space_temperature_ptrtype M_Xh;
    element_temperature_ptrtype M_fieldTemperature;
    bool M_fieldVelocityConvectionIsUsed, M_fieldVelocityConvectionIsIncompressible;
    space_velocityconvection_ptrtype M_XhVelocityConvection;
    element_velocityconvection_ptrtype M_fieldVelocityConvection; // only define with convection effect
    boost::optional<vector_field_expression<nDim, 1, 2>> M_exprVelocityConvection;

    bdf_temperature_ptrtype M_bdfTemperature;

    // physical parameter
    space_scalar_P0_ptrtype M_XhScalarP0;
    thermalproperties_ptrtype M_thermalProperties;

    // boundary conditions
    map_scalar_field<2> M_bcDirichlet;
    map_scalar_field<2> M_bcNeumann;
    map_scalar_fields<2> M_bcRobin;
    map_scalar_field<2> M_volumicForcesProperties;

    // algebraic data/tools
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;

    // post-process
    export_ptrtype M_exporter;
    bool M_doExportAll, M_doExportVelocityConvection;
    std::vector<ModelMeasuresForces> M_postProcessMeasuresForces;
    context_temperature_ptrtype M_postProcessMeasuresContextTemperature;

    typedef boost::function<void( vector_ptrtype& F, bool buildCstPart )> updateSourceTermLinearPDE_function_type;
    updateSourceTermLinearPDE_function_type M_overwritemethod_updateSourceTermLinearPDE;
};

} // namespace FeelModels
} // namespace Feel

#endif /* FEELPP_THERMODYNAMICSBASE_HPP */
