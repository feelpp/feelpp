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
#include <feel/feelmodels/modelcore/modelphysics.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>

//#include <feel/feelmodels/heat/thermalpropertiesdescription.hpp>
#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>

#include <feel/feelmodels/modelcore/stabilizationglsparameterbase.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisTemperatureType>
class Heat : public ModelNumerical,
             public ModelPhysics<ConvexType::nDim>,
             public std::enable_shared_from_this< Heat<ConvexType,BasisTemperatureType> >,
             public MarkerManagementDirichletBC,
             public MarkerManagementNeumannBC,
             public MarkerManagementRobinBC
    {
    public:
        typedef ModelNumerical super_type;
        using size_type = typename super_type::size_type;
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
        typedef MaterialsProperties<mesh_type> materialsproperties_type;
        typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;
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

        // measure tools for points evaluation
        typedef MeasurePointsEvaluation<space_temperature_type> measure_points_evaluation_type;
        typedef std::shared_ptr<measure_points_evaluation_type> measure_points_evaluation_ptrtype;

        Heat( std::string const& prefix,
              std::string const& keyword = "heat",
              worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
              std::string const& subPrefix  = "",
              ModelBaseRepository const& modelRep = ModelBaseRepository() );

        std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"HeatMesh.path"); }
        //std::string physic() const { return "heat"; }
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
        materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
        materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
        void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }
        materialsproperties_ptrtype const& thermalProperties() const { return M_materialsProperties; } // DEPRECATED
        //thermalproperties_ptrtype const& thermalProperties() const { return M_thermalProperties; }
        //thermalproperties_ptrtype & thermalProperties() { return M_thermalProperties; }
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
        std::string const& timeStepping() const { return M_timeStepping; }
        bdf_temperature_ptrtype const& timeStepBdfTemperature() const { return M_bdfTemperature; }
        std::shared_ptr<TSBase> timeStepBase() { return this->timeStepBdfTemperature(); }
        std::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBdfTemperature(); }
        void startTimeStep();
        void updateTimeStep();
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
        void initInitialConditions();
        void initPostProcess() override;

        template <typename SymbExprType>
        auto symbolsExprFit( SymbExprType const& se ) const { return super_type::symbolsExprFit( se ); }

    public :
        void initAlgebraicFactory();

        void setMesh( mesh_ptrtype const& mesh ) { M_mesh = mesh; }

        BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
        int nBlockMatrixGraph() const { return 1; }
        void init( bool buildModelAlgebraicFactory=true );
        void updateForUseFunctionSpacesVelocityConvection();

        void exportResults() { this->exportResults( this->currentTime() ); }
        void exportResults( double time );
        template <typename SymbolsExpr>
        void exportResults( double time, SymbolsExpr const& symbolsExpr );

        void executePostProcessMeasures( double time );
        template <typename TupleFieldsType,typename SymbolsExpr>
        void executePostProcessMeasures( double time, TupleFieldsType const& tupleFields, SymbolsExpr const& symbolsExpr );
        FEELPP_DEPRECATED void exportMeasures( double time ) { this->executePostProcessMeasures( time ); }
        void setDoExportResults( bool b ) { if (M_exporter) M_exporter->setDoExport( b ); }

        void updateParameterValues();

        template <typename SymbolsExpr>
        void updateFields( SymbolsExpr const& symbolsExpr )
            {
                //this->thermalProperties()->updateFields( symbolsExpr );
            }

        auto allFields( std::string const& prefix = "" ) const
            {
                return hana::make_tuple( std::make_pair( prefixvm( prefix,"temperature" ),this->fieldTemperaturePtr() ),
                                         std::make_pair( prefixvm( prefix,"velocity-convection" ), this->fieldVelocityConvectionIsOperational()?this->fieldVelocityConvectionPtr() : element_velocityconvection_ptrtype() )
                                         );
            }

        template <typename SymbExprType>
        auto exprPostProcessExports( SymbExprType const& se, std::string const& prefix = "" ) const
            {
                //std::map<std::string,std::vector<std::tuple<ModelExpression, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExpr;
                typedef decltype(expr(typename ModelExpression::expr_scalar_type{},se) ) _expr_scalar_type;
                std::map<std::string,std::vector<std::tuple<_expr_scalar_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprScalar;

                typedef decltype(expr(ModelExpression{}.template expr<nDim,nDim>(),se)) _expr_tensor2_type;
                std::map<std::string,std::vector<std::tuple<_expr_tensor2_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprTensor2;

                auto Id = eye<nDim,nDim>();
                typedef decltype(expr(typename ModelExpression::expr_scalar_type{},se)*Id) _expr_tensor2_from_scalar_type;
                std::map<std::string,std::vector<std::tuple<_expr_tensor2_from_scalar_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprTensor2FromScalar;

                bool onlyScalarThermalConductivity = this->materialsProperties()->allThermalConductivitiesAreScalar();

                for ( auto const& rangeData : this->materialsProperties()->rangeMeshElementsByMaterial() )
                {
                    std::string const& _matName = rangeData.first;
                    auto const& range = rangeData.second;
                    if ( this->materialsProperties()->hasThermalConductivity( _matName ) )
                    {
                        auto const& thermalConductivity = this->materialsProperties()->thermalConductivity( _matName );
                        if ( thermalConductivity.isMatrix() )
                        {
                            auto thermalConductivityExpr = expr( thermalConductivity.template expr<nDim,nDim>(), se);
                            mapExprTensor2[prefixvm(prefix,"thermal-conductivity")].push_back( std::make_tuple( thermalConductivityExpr, range, "element" ) );
                        }
                        else
                        {
                            auto thermalConductivityExpr = expr( thermalConductivity.exprScalar(), se );
                            if ( onlyScalarThermalConductivity )
                            {
                                mapExprScalar[prefixvm(prefix,"thermal-conductivity")].push_back( std::make_tuple( thermalConductivityExpr, range, "element" ) );
                            }
                            else
                            {
                                mapExprTensor2FromScalar[prefixvm(prefix,"thermal-conductivity")].push_back( std::make_tuple( thermalConductivityExpr*Id, range, "element" ) );
                            }
                        }
                    }
                    //mapExpr[prefixvm(prefix,"thermal-conductivity")].push_back( std::make_tuple( thermalConductivity, range, "element" ) );

                    if ( this->materialsProperties()->hasRho( _matName ) )
                    {
                        auto densityExpr = expr( this->materialsProperties()->rho( _matName ).expr(), se );
                        mapExprScalar[prefixvm(prefix,"density")].push_back( std::make_tuple( densityExpr, range, "element" ) );
                    }
                }
                return hana::make_tuple( mapExprScalar,mapExprTensor2,mapExprTensor2FromScalar );
            }

        auto symbolsExpr( std::string const& prefix_symbol = "heat_" ) const { return this->symbolsExpr( this->fieldTemperature(), prefix_symbol ); }
        template <typename FieldTemperatureType>
        auto symbolsExpr( FieldTemperatureType const& t, std::string const& prefix_symbol = "heat_" ) const
        {
            auto seField = this->symbolsExprField( t, prefix_symbol );
            auto seFit = this->symbolsExprFit( seField );
            auto seMat = this->symbolsExprMaterial( Feel::vf::symbolsExpr( seField, seFit ), prefix_symbol );
            return Feel::vf::symbolsExpr( seField, seFit, seMat );
        }

        auto symbolsExprField( std::string const& prefix_symbol = "heat_" ) const { return this->symbolsExprField( this->fieldTemperature(), prefix_symbol ); }
        template <typename FieldTemperatureType>
        auto symbolsExprField( FieldTemperatureType const& t, std::string const& prefix_symbol = "heat_" ) const
            {
                // generate symbols heat_T, heat_grad_T(_x,_y,_z), heat_dnT
                return Feel::vf::symbolsExpr( symbolExpr( (boost::format("%1%T")%prefix_symbol).str(),idv(t) ),
                                              symbolExpr( (boost::format("%1%grad_T")%prefix_symbol).str(),gradv(t), SymbolExprComponentSuffix( 1,nDim,true ) ),
                                              symbolExpr( (boost::format("%1%dn_T")%prefix_symbol).str(),dnv(t) )
                                              );
            }


        template <typename SymbExprType>
        auto symbolsExprMaterial( SymbExprType const& se, std::string const& prefix_symbol = "heat_" ) const
        {
            return this->materialsProperties()->symbolsExpr( se, prefix_symbol );
        }

        //___________________________________________________________________________________//
        //___________________________________________________________________________________//
        // apply assembly and solver
        /*virtual*/ void solve();

        void updateLinearPDE( DataUpdateLinear & data ) const override;
        template <typename SymbolsExpr>
        void updateLinearPDE( DataUpdateLinear & data, SymbolsExpr const& symbolsExpr ) const;
        template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType>
        void updateLinearPDEStabilizationGLS( Expr<RhoCpExprType> const& rhocp, Expr<ConductivityExprType> const& kappa,
                                              Expr<ConvectionExprType> const& uconv, RangeType const& range, DataUpdateLinear & data ) const;
        void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;

        // non linear (newton)
        void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;

        void updateJacobian( DataUpdateJacobian & data ) const override;
        template <typename SymbolsExpr>
        void updateJacobian( DataUpdateJacobian & data, SymbolsExpr const& symbolsExpr ) const;
        template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType>
        void updateJacobianStabilizationGLS( Expr<RhoCpExprType> const& rhocp, Expr<ConductivityExprType> const& kappa,
                                             Expr<ConvectionExprType> const& uconv, RangeType const& range, DataUpdateJacobian & data ) const;
        void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;

        void updateResidual( DataUpdateResidual & data ) const override;
        template <typename SymbolsExpr>
        void updateResidual( DataUpdateResidual & data, SymbolsExpr const& symbolsExpr ) const;
        template <typename RhoCpExprType,typename ConductivityExprType,typename ConvectionExprType,typename RangeType,typename... ExprT>
        void updateResidualStabilizationGLS( Expr<RhoCpExprType> const& rhocp, Expr<ConductivityExprType> const& kappa,
                                             Expr<ConvectionExprType> const& uconv, RangeType const& range, DataUpdateResidual & data,
                                             const ExprT&... exprs ) const;
        void updateResidualDofElimination( DataUpdateResidual & data ) const override;

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

    private :
        void updateTimeStepCurrentResidual();

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

        // time discretisation
        std::string M_timeStepping;
        bdf_temperature_ptrtype M_bdfTemperature;
        double M_timeStepThetaValue;
        vector_ptrtype M_timeStepThetaSchemePreviousContrib;

        // physical parameter
        space_scalar_P0_ptrtype M_XhScalarP0;
        //thermalproperties_ptrtype M_thermalProperties;
        materialsproperties_ptrtype M_materialsProperties;

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

        // post-process
        export_ptrtype M_exporter;
        std::map<std::string,ModelMeasuresNormalFluxGeneric> M_postProcessMeasuresNormalHeatFlux;
        measure_points_evaluation_ptrtype M_measurePointsEvaluation;
    };


template< typename ConvexType, typename BasisTemperatureType>
template <typename SymbolsExpr>
void
Heat<ConvexType,BasisTemperatureType>::exportResults( double time, SymbolsExpr const& symbolsExpr )
{
    this->log("Heat","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    auto fields = this->allFields();
    this->executePostProcessExports( M_exporter, time, fields, symbolsExpr, this->exprPostProcessExports( symbolsExpr ) );
    this->executePostProcessMeasures( time, fields, symbolsExpr );
    this->executePostProcessSave( (this->isStationary())? invalid_uint32_type_value : M_bdfTemperature->iteration(), fields );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("Heat","exportResults", "finish");
}

template< typename ConvexType, typename BasisTemperatureType>
template <typename TupleFieldsType, typename SymbolsExpr>
void
Heat<ConvexType,BasisTemperatureType>::executePostProcessMeasures( double time, TupleFieldsType const& tupleFields, SymbolsExpr const& symbolsExpr )
{
    bool hasMeasure = false;

    // compute measures
    for ( auto const& [ppName,ppFlux] : M_postProcessMeasuresNormalHeatFlux )
    {
        auto const& u = this->fieldTemperature();
        double signFlux = (ppFlux.isOutward())? -1.0 : 1.0;
        auto kappa = this->materialsProperties()->thermalConductivityExpr( symbolsExpr );
        double heatFlux = 0;
        if constexpr ( std::decay_t<decltype(kappa)>::template evaluator_t<typename mesh_type::element_type>::shape::is_scalar )
                         heatFlux = integrate(_range=markedfaces(this->mesh(),ppFlux.markers() ),
                                              _expr=signFlux*kappa*gradv(u)*N() ).evaluate()(0,0);
        else
            heatFlux = integrate(_range=markedfaces(this->mesh(),ppFlux.markers() ),
                                 _expr=signFlux*inner(kappa*trans(gradv(u)),N()) ).evaluate()(0,0);
        this->postProcessMeasuresIO().setMeasure("Normal_Heat_Flux_"+ppName,heatFlux);
        hasMeasure = true;
    }

    bool hasMeasureNorm = this->updatePostProcessMeasuresNorm( this->mesh(), M_rangeMeshElements, symbolsExpr, tupleFields );
    bool hasMeasureStatistics = this->updatePostProcessMeasuresStatistics( this->mesh(), M_rangeMeshElements, symbolsExpr, tupleFields );
    bool hasMeasurePoint = this->updatePostProcessMeasuresPoint( M_measurePointsEvaluation, tupleFields );
    if ( hasMeasureNorm || hasMeasureStatistics || hasMeasurePoint )
        hasMeasure = true;

    if ( hasMeasure )
    {
        if ( !this->isStationary() )
            this->postProcessMeasuresIO().setMeasure( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
        this->upload( this->postProcessMeasuresIO().pathFile() );
    }
}

} // namespace FeelModels
} // namespace Feel

#include <feel/feelmodels/heat/heatassembly.hpp>
#include <feel/feelmodels/heat/heatupdatestabilizationgls.hpp>


#endif /* FEELPP_TOOLBOXES_HEAT_HPP */
