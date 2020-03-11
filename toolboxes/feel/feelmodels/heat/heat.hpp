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
        using size_type = typename super_type::size_type;
        typedef Heat<ConvexType,BasisTemperatureType> self_type;
        typedef std::shared_ptr<self_type> self_ptrtype;
        //___________________________________________________________________________________//
        // mesh
        typedef ConvexType convex_type;
        static const uint16_type nDim = convex_type::nDim;
        static const uint16_type nOrderGeo = convex_type::nOrder;
        typedef Mesh<convex_type,double,0,uint32_type> mesh_type;
        typedef std::shared_ptr<mesh_type> mesh_ptrtype;
        // basis
        static const uint16_type nOrderTemperature = BasisTemperatureType::nOrder;
        static const uint16_type nOrderPoly = nOrderTemperature;
        typedef BasisTemperatureType basis_temperature_type;
        typedef Lagrange<nOrderPoly, Vectorial,Continuous,PointSetFekete,0> basis_velocityconvection_type;
        // function space temperature
        typedef FunctionSpace<mesh_type, bases<basis_temperature_type>, double, Periodicity <NoPeriodicity>,mortars<NoMortar> > space_temperature_type;
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

        // measure tools for points evaluation
        typedef MeasurePointsEvaluation<space_temperature_type> measure_points_evaluation_type;
        typedef std::shared_ptr<measure_points_evaluation_type> measure_points_evaluation_ptrtype;

        Heat( std::string const& prefix,
              std::string const& keyword = "heat",
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

        template <typename FieldTemperatureType>
        constexpr auto symbolsExprField( FieldTemperatureType const& t, hana::int_<2> /**/ ) const
            {
                return Feel::vf::symbolsExpr( symbolExpr("heat_T",idv(t) ),
                                              symbolExpr("heat_dxT",dxv(t) ),
                                              symbolExpr("heat_dyT",dyv(t) ),
                                              symbolExpr("heat_dnT",dnv(t) )
                                              );
            }
        template <typename FieldTemperatureType>
        constexpr auto symbolsExprField( FieldTemperatureType const& t, hana::int_<3> /**/ ) const
            {
                return Feel::vf::symbolsExpr( symbolExpr("heat_T",idv(t) ),
                                              symbolExpr("heat_dxT",dxv(t) ),
                                              symbolExpr("heat_dyT",dyv(t) ),
                                              symbolExpr("heat_dzT",dzv(t) ),
                                              symbolExpr("heat_dnT",dnv(t) )
                                              );
            }
        //auto symbolsExprFit() const { return super_type::symbolsExprFit( this->symbolsExprField() ); }

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
                this->thermalProperties()->updateFields( symbolsExpr );
            }

        auto allFields() const
            {
                return hana::make_tuple( std::make_pair( "temperature",this->fieldTemperaturePtr() ),
                                         std::make_pair( "velocity-convection", this->fieldVelocityConvectionIsOperational()?this->fieldVelocityConvectionPtr() : element_velocityconvection_ptrtype() ),
                                         std::make_pair( "thermal-conductivity",this->thermalProperties()->fieldThermalConductivityPtr() ),
                                         std::make_pair( "density", this->thermalProperties()->fieldRhoPtr() )
                                         );
            }

        template <typename FieldTemperatureType>
        /*constexpr*/auto symbolsExpr( FieldTemperatureType const& t ) const
        {
            auto seField = this->symbolsExprField( t );
            auto seFit = this->symbolsExprFit( seField );
            auto seMat = this->symbolsExprMaterial( Feel::vf::symbolsExpr( seField, seFit ) );
            return Feel::vf::symbolsExpr( seField, seFit, seMat );
        }
        auto symbolsExpr() const { return this->symbolsExpr( this->fieldTemperature() ); }

        constexpr auto symbolsExprField() const { return this->symbolsExprField( this->fieldTemperature() ); }
        template <typename FieldTemperatureType>
        constexpr auto symbolsExprField( FieldTemperatureType const& t ) const { return this->symbolsExprField( t, hana::int_<nDim>() ); }

        template <typename SymbExprType>
        auto symbolsExprMaterial( SymbExprType const& se ) const
        {
            typedef decltype(expr(scalar_field_expression<2>{},se)) _expr_scalar_type;
            std::vector<std::pair<std::string,_expr_scalar_type>> matPropSymbsScalar;
            typedef decltype(expr(matrix_field_expression<nDim,nDim,2>{},se)(0,0)) _expr_matrix_comp_type;
            std::vector<std::pair<std::string,_expr_matrix_comp_type>> matPropSymbsMatrixComp;
            for ( auto const& rangeData : this->thermalProperties()->rangeMeshElementsByMaterial() )
            {
                std::string const& _matName = rangeData.first;
                auto const& thermalConductivity = this->thermalProperties()->thermalConductivity( _matName );
                if ( thermalConductivity.isMatrix() )
                {
#if 0
                    // generate compilation error because need to fix/improve the return type of expr.evaluate(bool, worldcomm_ptr_t)
                    for ( int i=0;i<nDim;++i )
                        for ( int j=0;j<nDim;++j )
                            matPropSymbsMatrixComp.push_back( std::make_pair( (boost::format("heat_%1%_%2%%3%")%_matName%i%j).str(), expr( thermalConductivity.template exprMatrix<nDim,nDim>(), se )(i,j) ) );
#endif
                }
                else
                    matPropSymbsScalar.push_back( std::make_pair( (boost::format("heat_%1%_k")% _matName).str(), expr( thermalConductivity.exprScalar(), se ) ) );
            }
            return Feel::vf::symbolsExpr( symbolExpr( matPropSymbsScalar )/*, symbolExpr( matPropSymbsMatrixComp )*/ );
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

        void assembleLinear() { this->algebraicFactory()->assembleLinear(this->blockVectorSolution().vectorMonolithic()); }
        vector_ptrtype rhs() { return this->algebraicFactory()->rhs(); }
        sparse_matrix_ptrtype matrix() { return this->algebraicFactory()->matrix(); }

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

        // post-process
        export_ptrtype M_exporter;
        bool M_doExportAll, M_doExportVelocityConvection;
        std::vector< ModelMeasuresForces > M_postProcessMeasuresForces;
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
    this->executePostProcessExports( M_exporter, time, fields );
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
    for ( auto const& ppForces : M_postProcessMeasuresForces )
    {
        CHECK( ppForces.meshMarkers().size() == 1 ) << "TODO";
        auto const& u = this->fieldTemperature();
        auto kappa = idv(this->thermalProperties()->fieldThermalConductivity());
        double heatFlux = integrate(_range=markedfaces(this->mesh(),ppForces.meshMarkers() ),
                                    _expr=kappa*gradv(u)*N() ).evaluate()(0,0);
        std::string name = ppForces.name();
        this->postProcessMeasuresIO().setMeasure("NormalHeatFlux_"+name,heatFlux);
        hasMeasure = true;
    }

    bool hasMeasureNorm = this->executePostProcessMeasuresNorm( this->mesh(), M_rangeMeshElements, tupleFields, symbolsExpr );
    bool hasMeasureStatistics = this->executePostProcessMeasuresStatistics( this->mesh(), M_rangeMeshElements, tupleFields, symbolsExpr );
    bool hasMeasurePoint = this->executePostProcessMeasuresPoint( M_measurePointsEvaluation, tupleFields );
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
