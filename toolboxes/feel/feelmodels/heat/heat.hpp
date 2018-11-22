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
   \file heat.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-06-04
 */


#ifndef FEELPP_TOOLBOXES_HEAT_HPP
#define FEELPP_TOOLBOXES_HEAT_HPP 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
//#include <feel/feelvf/vf.hpp>
#include <feel/feelts/bdf.hpp>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>

#include <feel/feelmodels/heat/thermalpropertiesdescription.hpp>

#include <feel/feelmodels/modelcore/stabilizationglsparameterbase.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisTemperatureType>
class Heat : public ModelNumerical,
                     public std::enable_shared_from_this< Heat<ConvexType,BasisTemperatureType> >,
                     public MarkerManagementDirichletBC,
                     public MarkerManagementNeumannBC,
                     public MarkerManagementRobinBC
    {
    public:
        typedef ModelNumerical super_type;
        typedef Heat<ConvexType,BasisTemperatureType> self_type;
        typedef std::shared_ptr<self_type> self_ptrtype;
        //___________________________________________________________________________________//
        // mesh
        typedef ConvexType convex_type;
        static const uint16_type nDim = convex_type::nDim;
        static const uint16_type nOrderGeo = convex_type::nOrder;
        typedef Mesh<convex_type> mesh_type;
        typedef std::shared_ptr<mesh_type> mesh_ptrtype;
        // basis
        static const uint16_type nOrderTemperature = BasisTemperatureType::nOrder;
        static const uint16_type nOrderPoly = nOrderTemperature;
        typedef BasisTemperatureType basis_temperature_type;
        typedef Lagrange<nOrderPoly, Vectorial,Continuous,PointSetFekete> basis_velocityconvection_type;
        // function space temperature
        typedef FunctionSpace<mesh_type, bases<basis_temperature_type> > space_temperature_type;
        typedef std::shared_ptr<space_temperature_type> space_temperature_ptrtype;
        typedef typename space_temperature_type::element_type element_temperature_type;
        typedef std::shared_ptr<element_temperature_type> element_temperature_ptrtype;
        typedef typename space_temperature_type::element_external_storage_type element_temperature_external_storage_type;
        // function space velocity convection
        typedef FunctionSpace<mesh_type, bases<basis_velocityconvection_type> > space_velocityconvection_type;
        typedef std::shared_ptr<space_velocityconvection_type> space_velocityconvection_ptrtype;
        typedef typename space_velocityconvection_type::element_type element_velocityconvection_type;
        typedef std::shared_ptr<element_velocityconvection_type> element_velocityconvection_ptrtype;
        // mechanical properties desc
        typedef bases<Lagrange<0, Scalar,Discontinuous> > basis_scalar_P0_type;
        typedef FunctionSpace<mesh_type, basis_scalar_P0_type> space_scalar_P0_type;
        typedef std::shared_ptr<space_scalar_P0_type> space_scalar_P0_ptrtype;
        typedef ThermalPropertiesDescription<space_scalar_P0_type> thermalproperties_type;
        typedef std::shared_ptr<thermalproperties_type> thermalproperties_ptrtype;
        // time scheme
        typedef Bdf<space_temperature_type>  bdf_temperature_type;
        typedef std::shared_ptr<bdf_temperature_type> bdf_temperature_ptrtype;
        // stabilization
        typedef StabilizationGLSParameterBase<mesh_type> stab_gls_parameter_type;
        typedef std::shared_ptr<stab_gls_parameter_type> stab_gls_parameter_ptrtype;
        // exporter
        typedef Exporter<mesh_type,nOrderGeo> export_type;
        typedef std::shared_ptr<export_type> export_ptrtype;

        // algebraic solver
        typedef ModelAlgebraicFactory model_algebraic_factory_type;
        typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

        // context for evaluation
        typedef typename space_temperature_type::Context context_temperature_type;
        typedef std::shared_ptr<context_temperature_type> context_temperature_ptrtype;


        Heat( std::string const& prefix,
                      bool buildMesh = true,
                      worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                      std::string const& subPrefix  = "",
                      ModelBaseRepository const& modelRep = ModelBaseRepository() );

        std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"HeatMesh.path"); }
        //___________________________________________________________________________________//
        // mesh, space, element temperature
        mesh_ptrtype const& mesh() const { return M_mesh; }
        elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

        space_temperature_ptrtype const& spaceTemperature() const { return M_Xh; }
        element_temperature_ptrtype const& fieldTemperaturePtr() const { return M_fieldTemperature; }
        element_temperature_type const& fieldTemperature() const { return *M_fieldTemperature; }
        element_velocityconvection_ptrtype const& fieldVelocityConvectionPtr() const { return M_fieldVelocityConvection; }
        element_velocityconvection_ptrtype & fieldVelocityConvectionPtr() { return M_fieldVelocityConvection; }
        element_velocityconvection_type const& fieldVelocityConvection() const { return *M_fieldVelocityConvection; }
        bool fieldVelocityConvectionIsUsed() const { return M_fieldVelocityConvectionIsUsed; }
        bool fieldVelocityConvectionIsIncompressible() const { return M_fieldVelocityConvectionIsIncompressible; }
        void setFieldVelocityConvectionIsUsed(bool b) { M_fieldVelocityConvectionIsUsed=b; }
        bool fieldVelocityConvectionIsOperational() const { return (M_fieldVelocityConvection.use_count() > 0); }
        bool fieldVelocityConvectionIsUsedAndOperational() const { return this->fieldVelocityConvectionIsUsed() && this->fieldVelocityConvectionIsOperational(); }
        void setFieldVelocityConvectionIsIncompressible(bool b) { M_fieldVelocityConvectionIsIncompressible=b; }
        // stabilization
        bool stabilizationGLS() const { return M_stabilizationGLS; }
        std::string const& stabilizationGLSType() const { return M_stabilizationGLSType; }
        stab_gls_parameter_ptrtype const& stabilizationGLSParameter() const { return M_stabilizationGLSParameter; }
        //___________________________________________________________________________________//
        // physical parameters
        thermalproperties_ptrtype const& thermalProperties() const { return M_thermalProperties; }
        thermalproperties_ptrtype & thermalProperties() { return M_thermalProperties; }
        // boundary condition + body forces
        map_scalar_field<2> const& bcDirichlet() const { return M_bcDirichlet; }
        map_scalar_field<2> const& bcNeumann() const { return M_bcNeumann; }
        map_scalar_fields<2> const& bcRobin() const { return M_bcRobin; }
        map_scalar_field<2> const& bodyForces() const { return M_volumicForcesProperties; }
        //___________________________________________________________________________________//
        // algebraic data and solver
        backend_ptrtype const& backend() const { return  M_backend; }
        BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
        BlocksBaseVector<double> & blockVectorSolution() { return M_blockVectorSolution; }
        size_type nLocalDof() const;
        model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
        model_algebraic_factory_ptrtype & algebraicFactory() { return M_algebraicFactory; }
        //___________________________________________________________________________________//
        // time step scheme
        bdf_temperature_ptrtype const& timeStepBdfTemperature() const { return M_bdfTemperature; }
        std::shared_ptr<TSBase> timeStepBase() { return this->timeStepBdfTemperature(); }
        std::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBdfTemperature(); }
        void updateBdf();
        void updateTimeStep() { this->updateBdf(); }
        //___________________________________________________________________________________//

        std::shared_ptr<std::ostringstream> getInfo() const override;
        void updateInformationObject( pt::ptree & p ) override;
    private :
        void loadParameterFromOptionsVm();
        void initMesh();
        void initMaterialProperties();
        void initFunctionSpaces();
        void initBoundaryConditions();
        void initTimeStep();
        void initPostProcess();

        constexpr auto symbolsExpr( hana::int_<2> /**/ ) const
            {
                return Feel::vf::symbolsExpr( symbolExpr("heat_T",idv(this->fieldTemperature()) ),
                                              symbolExpr("heat_dxT",dxv(this->fieldTemperature()) ),
                                              symbolExpr("heat_dyT",dyv(this->fieldTemperature()) )
                                              );
            }
        constexpr auto symbolsExpr( hana::int_<3> /**/ ) const
            {
                return Feel::vf::symbolsExpr( symbolExpr("heat_T",idv(this->fieldTemperature()) ),
                                              symbolExpr("heat_dxT",dxv(this->fieldTemperature()) ),
                                              symbolExpr("heat_dyT",dyv(this->fieldTemperature()) ),
                                              symbolExpr("heat_dzT",dzv(this->fieldTemperature()) )
                                              );
            }
    public :
        void initAlgebraicFactory();

        void setMesh( mesh_ptrtype const& mesh ) { M_mesh = mesh; }

        BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
        int nBlockMatrixGraph() const { return 1; }
        void init( bool buildModelAlgebraicFactory=true );
        void updateForUseFunctionSpacesVelocityConvection();

        std::set<std::string> postProcessFieldExported( std::set<std::string> const& ifields, std::string const& prefix = "" ) const;
        bool hasPostProcessFieldExported( std::string const& fieldName ) const { return M_postProcessFieldExported.find( fieldName ) != M_postProcessFieldExported.end(); }

        void exportResults() { this->exportResults( this->currentTime() ); }
        void exportResults( double time );
        void exportFields( double time );
        bool updateExportedFields( export_ptrtype exporter, std::set<std::string> const& fields, double time );
        void exportMeasures( double time );
        void setDoExportResults( bool b ) { if (M_exporter) M_exporter->setDoExport( b ); }

        void updateParameterValues();
        constexpr auto symbolsExpr() const { return this->symbolsExpr( hana::int_<nDim>() ); }
        //___________________________________________________________________________________//
        //___________________________________________________________________________________//
        // apply assembly and solver
        /*virtual*/ void solve();

        void updateLinearPDE( DataUpdateLinear & data ) const override;
        //void updateLinearPDEStabilizationGLS( DataUpdateLinear & data ) const;
        void updateLinearPDEWeakBC( sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const;
        void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;
        void updateLinearPDESourceTerm( vector_ptrtype& F, bool buildCstPart) const;
        template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType>
        void updateLinearPDEStabilizationGLS( Expr<RhoCpExprType> const& rhocp, Expr<ConductivityExprType> const& kappa,
                                              Expr<ConvectionExprType> const& uconv, RangeType const& range, DataUpdateLinear & data ) const;

        // non linear (newton)
        void updateNewtonInitialGuess( vector_ptrtype& U ) const override;
        void updateJacobian( DataUpdateJacobian & data ) const override;
        void updateJacobianRobinBC( sparse_matrix_ptrtype& J, bool buildCstPart ) const;
        void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
        template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType>
        void updateJacobianStabilizationGLS( Expr<RhoCpExprType> const& rhocp, Expr<ConductivityExprType> const& kappa,
                                             Expr<ConvectionExprType> const& uconv, RangeType const& range, DataUpdateJacobian & data ) const;
        void updateResidual( DataUpdateResidual & data ) const override;
        void updateResidualSourceTerm( vector_ptrtype& R, bool buildCstPart ) const;
        void updateResidualNeumannBC( vector_ptrtype& R, bool buildCstPart ) const;
        void updateResidualRobinBC( element_temperature_external_storage_type const& u, vector_ptrtype& R, bool buildCstPart ) const;
        void updateResidualDofElimination( DataUpdateResidual & data ) const override;
        template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType>
        void updateResidualStabilizationGLS( Expr<RhoCpExprType> const& rhocp, Expr<ConductivityExprType> const& kappa,
                                             Expr<ConvectionExprType> const& uconv, RangeType const& range, DataUpdateResidual & data ) const;

        //___________________________________________________________________________________//
        //___________________________________________________________________________________//
        // update field from expr
        void updateFieldVelocityConvection( bool onlyExprWithTimeSymbol = false );
        template < typename ExprT >
        void updateFieldVelocityConvection( vf::Expr<ExprT> const& expr )
        {
            this->updateFieldVelocityConvection( elements(this->mesh()), expr );
        }
        template < typename ExprT >
        void updateFieldVelocityConvection( elements_reference_wrapper_t<mesh_type> const& range, vf::Expr<ExprT> const& expr )
        {
            if ( !M_fieldVelocityConvection )
                this->updateForUseFunctionSpacesVelocityConvection();
            M_exprVelocityConvection.reset();// symbolic expression is remove
            M_fieldVelocityConvection->on(_range=range, _expr=expr );
        }

    protected :

        bool M_hasBuildFromMesh, M_isUpdatedForUse;

        mesh_ptrtype M_mesh;
        elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

        space_temperature_ptrtype M_Xh;
        element_temperature_ptrtype M_fieldTemperature;
        bool M_fieldVelocityConvectionIsUsed, M_fieldVelocityConvectionIsIncompressible;
        space_velocityconvection_ptrtype M_XhVelocityConvection;
        element_velocityconvection_ptrtype M_fieldVelocityConvection; // only define with convection effect
        boost::optional<vector_field_expression<nDim,1,2> > M_exprVelocityConvection;

        bdf_temperature_ptrtype M_bdfTemperature;

        // physical parameter
        space_scalar_P0_ptrtype M_XhScalarP0;
        thermalproperties_ptrtype M_thermalProperties;

        // boundary conditions
        map_scalar_field<2> M_bcDirichlet;
        map_scalar_field<2> M_bcNeumann;
        map_scalar_fields<2> M_bcRobin;
        map_scalar_field<2> M_volumicForcesProperties;

        // stabilization
        bool M_stabilizationGLS;
        std::string M_stabilizationGLSType;
        stab_gls_parameter_ptrtype M_stabilizationGLSParameter;

        // algebraic data/tools
        backend_ptrtype M_backend;
        model_algebraic_factory_ptrtype M_algebraicFactory;
        BlocksBaseVector<double> M_blockVectorSolution;
        std::map<std::string,std::set<size_type> > M_dofsWithValueImposed;

        // post-process
        std::set<std::string> M_postProcessFieldExported;
        export_ptrtype M_exporter;
        bool M_doExportAll, M_doExportVelocityConvection;
        std::vector< ModelMeasuresForces > M_postProcessMeasuresForces;
        context_temperature_ptrtype M_postProcessMeasuresContextTemperature;


        typedef boost::function<void ( vector_ptrtype& F, bool buildCstPart )> updateSourceTermLinearPDE_function_type;
        updateSourceTermLinearPDE_function_type M_overwritemethod_updateSourceTermLinearPDE;

    };

} // namespace FeelModels
} // namespace Feel

#include <feel/feelmodels/heat/heatupdatestabilizationgls.hpp>

#endif /* FEELPP_TOOLBOXES_HEAT_HPP */
