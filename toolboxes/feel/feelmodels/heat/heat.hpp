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

#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>

#include <feel/feelmodels/modelcore/stabilizationglsparameterbase.hpp>

#include <feel/feelmodels/heat/heatboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisTemperatureType>
class Heat : public ModelNumerical,
             public ModelPhysics<ConvexType::nDim>,
             public std::enable_shared_from_this< Heat<ConvexType,BasisTemperatureType> >
    {
        typedef ModelPhysics<ConvexType::nDim> super_physics_type;
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
        static const uint16_type nRealDim = convex_type::nRealDim;
        typedef Mesh<convex_type> mesh_type;
        typedef std::shared_ptr<mesh_type> mesh_ptrtype;
        // basis
        static const uint16_type nOrderTemperature = BasisTemperatureType::nOrder;
        static const uint16_type nOrderPoly = nOrderTemperature;
        typedef BasisTemperatureType basis_temperature_type;
        // function space temperature
        typedef FunctionSpace<mesh_type, bases<basis_temperature_type> > space_temperature_type;
        typedef std::shared_ptr<space_temperature_type> space_temperature_ptrtype;
        typedef typename space_temperature_type::element_type element_temperature_type;
        typedef std::shared_ptr<element_temperature_type> element_temperature_ptrtype;
        typedef typename space_temperature_type::element_external_storage_type element_temperature_external_storage_type;
        // materials properties
        typedef MaterialsProperties<nRealDim> materialsproperties_type;
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

        struct FieldTag
        {
            static auto temperature( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
        };

        template <typename ... Ts>
        static self_ptrtype New( Ts && ... v )
            {
                auto args = NA::make_arguments( std::forward<Ts>(v)... );
                std::string const& prefix = args.get(_prefix);
                std::string const& keyword = args.get_else(_keyword,"heat");
                worldcomm_ptr_t worldcomm = args.get_else(_worldcomm,Environment::worldCommPtr());
                auto && repository = args.get_else(_repository,ModelBaseRepository{});
                return std::make_shared<self_type>( prefix, keyword, worldcomm, "", repository );
            }


        Heat( std::string const& prefix,
              std::string const& keyword = "heat",
              worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
              std::string const& subPrefix  = "",
              ModelBaseRepository const& modelRep = ModelBaseRepository() );

        //___________________________________________________________________________________//
        // mesh, space, element temperature
        mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
        void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }
        elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

        void applyRemesh( mesh_ptrtype const& newMesh );

        space_temperature_ptrtype const& spaceTemperature() const { return M_Xh; }
        element_temperature_ptrtype const& fieldTemperaturePtr() const { return M_fieldTemperature; }
        element_temperature_type const& fieldTemperature() const { return *M_fieldTemperature; }

        // stabilization
        bool stabilizationGLS() const { return M_stabilizationGLS; }
        std::string const& stabilizationGLSType() const { return M_stabilizationGLSType; }
        stab_gls_parameter_ptrtype const& stabilizationGLSParameter() const { return M_stabilizationGLSParameter; }
        bool stabilizationGLS_checkConductivityDependencyOnCoordinates() const { return M_stabilizationGLS_checkConductivityDependencyOnCoordinates; }

        //___________________________________________________________________________________//
        // physical parameters
        materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
        materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
        void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

        //___________________________________________________________________________________//
        // time step scheme
        std::string const& timeStepping() const { return M_timeStepping; }
        bdf_temperature_ptrtype const& timeStepBdfTemperature() const { return M_bdfTemperature; }
        std::shared_ptr<TSBase> timeStepBase() { return this->timeStepBdfTemperature(); }
        std::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBdfTemperature(); }
        void startTimeStep();
        void updateTimeStep();

        template <typename SymbolsExprType>
        void updateInitialConditions( SymbolsExprType const& se );
        //___________________________________________________________________________________//

        std::shared_ptr<std::ostringstream> getInfo() const override;
        void updateInformationObject( nl::json & p ) const override;
        tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

    private :
        void loadParameterFromOptionsVm();
        void initMesh();
        void initMaterialProperties();
        void initFunctionSpaces();
        void initBoundaryConditions();
        void initTimeStep();
        void initPostProcess() override;

        void initAlgebraicModel();
        void updateAlgebraicDofEliminationIds();

    public :
        void initAlgebraicFactory();

        BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
        int nBlockMatrixGraph() const { return 1; }
        void init( bool buildModelAlgebraicFactory=true );

        void updateParameterValues();
        void setParameterValues( std::map<std::string,double> const& paramValues );

        //___________________________________________________________________________________//
        // execute post-processing
        //___________________________________________________________________________________//

        void exportResults() { this->exportResults( this->currentTime() ); }
        void exportResults( double time );

        template <typename ModelFieldsType,typename SymbolsExpr,typename ExportsExprType>
        void exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr );

        template <typename SymbolsExpr>
        void exportResults( double time, SymbolsExpr const& symbolsExpr )
            {
                return this->exportResults( time, this->modelFields(), symbolsExpr, this->exprPostProcessExports( symbolsExpr ) );
            }

        template <typename ModelFieldsType,typename SymbolsExpr, typename ModelMeasuresQuantitiesType>
        void executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ModelMeasuresQuantitiesType const& mquantities );

        bool checkResults() const override;
        //___________________________________________________________________________________//
        // export expressions
        //___________________________________________________________________________________//

        template <typename SymbExprType>
        auto exprPostProcessExportsToolbox( SymbExprType const& se, std::string const& prefix ) const
            {
                using _expr_velocity_convection_type = std::decay_t<decltype( std::declval<ModelPhysicHeat<nDim>>().convection().expr( se ) )>;
                std::map<std::string,std::vector<std::tuple< _expr_velocity_convection_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprVelocityConvection;

                for ( auto const& [physicId,physicData] : this->physicsFromCurrentType() )
                {
                    auto physicHeatData = std::static_pointer_cast<ModelPhysicHeat<nDim>>(physicData);
                    for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicId ) )
                    {
                        auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                        if ( physicHeatData->hasConvectionEnabled() )
                        {
                             auto velocityConvectionExpr = physicHeatData->convection().expr( se );
                             mapExprVelocityConvection[prefixvm(prefix,"velocity-convection")].push_back( std::make_tuple( velocityConvectionExpr, range, "nodal" ) );
                        }
                    }
                }
                return hana::make_tuple( mapExprVelocityConvection );
            }
        template <typename SymbExprType>
        auto exprPostProcessExports( SymbExprType const& se, std::string const& prefix = "" ) const
            {
                return hana::concat( this->materialsProperties()->exprPostProcessExports( this->mesh(),this->physicsAvailable(),se ),
                                     this->exprPostProcessExportsToolbox( se, prefix ) );
            }
        //___________________________________________________________________________________//
        // toolbox fields
        //___________________________________________________________________________________//

        auto modelFields( std::string const& prefix = "" ) const
            {
                return this->modelFields( this->fieldTemperaturePtr(), prefix );
            }
        auto modelFields( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
            {
                auto field_t = this->spaceTemperature()->elementPtr( *sol, rowStartInVector + this->startSubBlockSpaceIndex( "temperature" ) );
                return this->modelFields( field_t, prefix );
            }
        auto modelFields( std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorData, std::string const& prefix = "" ) const
            {
                auto itFindSolution = vectorData.find( "solution" );
                CHECK( itFindSolution != vectorData.end() ) << "require solution data";
                vector_ptrtype sol = std::get<0>( itFindSolution->second );
                size_type rowStartInVector =  std::get<1>( itFindSolution->second );
                auto field_t = this->spaceTemperature()->elementPtr( *sol, rowStartInVector + this->startSubBlockSpaceIndex( "temperature" ) );
                return this->modelFields( field_t, prefix );
            }
        template <typename TemperatureFieldType>
        auto modelFields( TemperatureFieldType const& field_t, std::string const& prefix = "" ) const
            {
                return Feel::FeelModels::modelFields( modelField<FieldCtx::FULL>( FieldTag::temperature(this), prefix, "temperature", field_t, "T", this->keyword() ) );
            }

        auto trialSelectorModelFields( size_type startBlockSpaceIndex = 0 ) const
            {
                return Feel::FeelModels::selectorModelFields( selectorModelField( FieldTag::temperature(this), "temperature", startBlockSpaceIndex ) );
            }

        //___________________________________________________________________________________//
        // symbols expressions
        //___________________________________________________________________________________//

        template <typename ModelFieldsType>
        auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            auto seToolbox = this->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            auto seMeshes = this->template symbolsExprMeshes<mesh_type>();
            auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr(); // generate symbols heat_T, heat_grad_T(_x,_y,_z), heat_dn_T
            return Feel::vf::symbolsExpr( seToolbox, seParam, seMeshes, seMat, seFields );
        }
        auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

        template <typename ModelFieldsType>
        auto symbolsExprToolbox( ModelFieldsType const& mfields ) const
            {
                auto const& t = mfields.field( FieldTag::temperature(this), "temperature" );
                // generate symbol heat_nflux
                typedef decltype( this->normalHeatFluxExpr(t) ) _expr_nflux_type;
                symbol_expression_t<_expr_nflux_type> se_nflux;
                std::string symbolNormalHeatFluxStr = prefixvm( this->keyword(), "nflux", "_");
                auto _normalHeatFluxExpr = this->normalHeatFluxExpr( t );
                se_nflux.add( symbolNormalHeatFluxStr, _normalHeatFluxExpr, SymbolExprComponentSuffix( 1,1 ) );

#if 0
                // velocity convection : on each material
                symbol_expression_t<velocity_convection_expr_type> se_vconv_bymat;
                for ( auto const& [matName,uExpr] : M_exprVelocityConvection )
                {
                    std::string symbolstr_vconv_bymat = (boost::format("%1%_%2%_vconv")%this->keyword() %matName).str();
                    se_vconv_bymat.add( symbolstr_vconv_bymat, uExpr, SymbolExprComponentSuffix( nDim,1 ) );
                }
                // velocity convection : for all materials
                typedef decltype( this->velocityConvectionExpr() ) _expr_vconv_type;
                symbol_expression_t<_expr_vconv_type> se_vconv;
                std::string symbolstr_vconv = prefixvm( this->keyword(), "vconv", "_");
                se_vconv.add( symbolstr_vconv, this->velocityConvectionExpr(), SymbolExprComponentSuffix( nDim,1 ) );
#endif
                return Feel::vf::symbolsExpr( se_nflux/*,se_vconv,se_vconv_bymat*/ );
            }

        template <typename ModelFieldsType, typename TrialSelectorModelFieldsType>
        auto trialSymbolsExpr( ModelFieldsType const& mfields, TrialSelectorModelFieldsType const& tsmf ) const
            {
                return mfields.trialSymbolsExpr( tsmf );
            }

        //___________________________________________________________________________________//
        // model context helper
        //___________________________________________________________________________________//

        template <typename ModelFieldsType>
        auto modelContext( ModelFieldsType const& mfields, std::string const& prefix = "" ) const
            {
                auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
                return Feel::FeelModels::modelContext( mfields, std::move( se ) );
            }
        auto modelContext( std::string const& prefix = "" ) const
            {
                auto mfields = this->modelFields( prefix );
                auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
                return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
            }
        auto modelContext( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
            {
                auto mfields = this->modelFields( sol, rowStartInVector, prefix );
                auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
                auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( rowStartInVector ) );
                return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ), std::move( tse ) );
            }
        auto modelContextNoTrialSymbolsExpr( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
            {
                auto mfields = this->modelFields( sol, rowStartInVector, prefix );
                auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
                return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
            }

        //___________________________________________________________________________________//
        // toolbox expressions
        //___________________________________________________________________________________//

        template <typename FieldTemperatureType, typename SymbolsExpr = symbols_expression_empty_t>
        auto normalHeatFluxExpr( FieldTemperatureType const& t, bool isOutward = true, SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
            {
                double signFlux = isOutward? -1.0 : 1.0;
                auto kappa = this->materialsProperties()->thermalConductivityExpr( symbolsExpr );
                if constexpr ( std::decay_t<decltype(kappa)>::template evaluator_t<typename mesh_type::element_type>::shape::is_scalar )
                                 return signFlux*kappa*gradv(t)*N();
                else
                    return signFlux*inner(kappa*trans(gradv(t)),N());
            }

#if 0
        velocity_convection_expr_type const& velocityConvectionExpr( std::string const& matName ) const
            {
                auto itFindVel = M_exprVelocityConvection.find( matName );
                CHECK( itFindVel != M_exprVelocityConvection.end() ) << "no velocity convection with material " << matName << std::endl;
                return itFindVel->second;
            }

        template <typename SymbolsExprType>
        auto velocityConvectionExpr( std::string const& matName, SymbolsExprType const& se ) const
            {
                return expr( this->velocityConvectionExpr( matName ), se );
            }

        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto velocityConvectionExpr( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
            {
                std::vector<std::pair<std::string,velocity_convection_expr_type>> theExprs;
                for ( auto const& [matName,uExpr] : M_exprVelocityConvection )
                    theExprs.push_back( std::make_pair( matName, uExpr ) );

                if constexpr ( std::is_same_v<SymbolsExprType,symbols_expression_empty_t> )
                    return expr<typename mesh_type::index_type>( this->materialsProperties()->exprSelectorByMeshElementMapping(), theExprs );
                else
                    return expr<typename mesh_type::index_type>( this->materialsProperties()->exprSelectorByMeshElementMapping(), theExprs ).applySymbolsExpr( se );
            };
#endif
        //___________________________________________________________________________________//
        // apply assembly and solver
        //___________________________________________________________________________________//

        void solve();

        void updateLinearPDE( DataUpdateLinear & data ) const override;
        template <typename ModelContextType>
        void updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mfields ) const;
        template <typename ModelContextType,typename RangeType>
        void updateLinearPDEStabilizationGLS(  DataUpdateLinear & data, ModelContextType const& mctx,
                                               ModelPhysicHeat<nDim> const& physicHeatData,
                                               MaterialProperties const& matProps, RangeType const& range ) const;
        void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;
        template <typename ModelContextType>
        void updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mfields ) const;

        // non linear (newton)
        void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
        template <typename ModelContextType>
        void updateNewtonInitialGuess( DataNewtonInitialGuess & data, ModelContextType const& mfields ) const;

        void updateJacobian( DataUpdateJacobian & data ) const override;
        template <typename ModelContextType>
        void updateJacobian( DataUpdateJacobian & data, ModelContextType const& mfields ) const;
        template <typename ModelContextType,typename RangeType>
        void updateJacobianStabilizationGLS( DataUpdateJacobian & data, ModelContextType const& mctx,
                                             ModelPhysicHeat<nDim> const& physicHeatData,
                                             MaterialProperties const& matProps, RangeType const& range ) const;
        void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;

        void updateResidual( DataUpdateResidual & data ) const override;
        template <typename ModelContextType>
        void updateResidual( DataUpdateResidual & data, ModelContextType const& mfields ) const;
        template <typename ModelContextType,typename RangeType,typename... ExprAddedType>
        void updateResidualStabilizationGLS( DataUpdateResidual & data, ModelContextType const& mctx,
                                             ModelPhysicHeat<nDim> const& physicHeatData,
                                             MaterialProperties const& matProps, RangeType const& range,
                                             const ExprAddedType&... exprsAddedInResidual ) const;
        void updateResidualDofElimination( DataUpdateResidual & data ) const override;

        //___________________________________________________________________________________//
        //___________________________________________________________________________________//
    private :
        void updateTimeStepCurrentResidual();

        auto modelMeasuresQuantities( std::string const& prefix = "" ) const
            {
#if 0
                return Feel::FeelModels::modelMeasuresQuantities( modelMeasuresQuantity( prefix, "q1", 3.1 ),
                                                                  modelMeasuresQuantity( prefix, "q2", std::vector<double>({3.8,3.9,3.10}) ),
                                                                  modelMeasuresQuantity( prefix, "q3", eigen_matrix_type<nDim, nDim>::Constant(5.0) )
                                                                  //modelMeasuresQuantity( prefix, "q4", std::bind( &self_type::volume, this ) )
                                                                  );
#endif
                return model_measures_quantities_empty_t{};
            }

    protected :

        elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

        space_temperature_ptrtype M_Xh;
        bool M_useExtendedDoftable = false;
        element_temperature_ptrtype M_fieldTemperature;

        // time discretisation
        std::string M_timeStepping;
        bdf_temperature_ptrtype M_bdfTemperature;
        double M_timeStepThetaValue;
        vector_ptrtype M_timeStepThetaSchemePreviousContrib;

        std::map<std::string,double> M_currentParameterValues;

        // physical parameter
        materialsproperties_ptrtype M_materialsProperties;

        // boundary conditions
        using boundary_conditions_type = HeatBoundaryConditions;
        std::shared_ptr<boundary_conditions_type> M_boundaryConditions;

        // stabilization
        bool M_stabilizationGLS;
        std::string M_stabilizationGLSType;
        stab_gls_parameter_ptrtype M_stabilizationGLSParameter;
        bool M_stabilizationGLS_checkConductivityDependencyOnCoordinates = true;

        std::string M_solverName;

        // post-process
        export_ptrtype M_exporter;
        std::map<std::string,ModelMeasuresNormalFluxGeneric> M_postProcessMeasuresNormalHeatFlux;
    };


template< typename ConvexType, typename BasisTemperatureType>
template <typename SymbolsExprType>
void
Heat<ConvexType,BasisTemperatureType>::updateInitialConditions( SymbolsExprType const& se )
{
    if ( !this->doRestart() )
    {
        std::vector<element_temperature_ptrtype> icTemperatureFields;
        std::map<int, double> icPriorTimes;
        if ( this->isStationary() )
        {
            icTemperatureFields = { this->fieldTemperaturePtr() };
            icPriorTimes = {{0,0}};
        }
        else
        {
            icTemperatureFields = M_bdfTemperature->unknowns();
            icPriorTimes = M_bdfTemperature->priorTimes();
        }

        super_type::updateInitialConditions( "temperature", M_rangeMeshElements, se, icTemperatureFields, icPriorTimes );

        if ( Environment::vm().count( prefixvm(this->prefix(),"initial-solution.temperature").c_str() ) )
        {
            auto myexpr = expr( soption(_prefix=this->prefix(),_name="initial-solution.temperature"),
                                "",this->worldComm(),this->repository().expr() );
            icTemperatureFields[0]->on(_range=M_rangeMeshElements,_expr=myexpr);
            for ( int k=1;k<icTemperatureFields.size();++k )
                *icTemperatureFields[k] = *icTemperatureFields[0];
        }

        if ( !this->isStationary() )
            *this->fieldTemperaturePtr() = M_bdfTemperature->unknown(0);
    }

}

template< typename ConvexType, typename BasisTemperatureType>
template <typename ModelFieldsType, typename SymbolsExpr, typename ExportsExprType>
void
Heat<ConvexType,BasisTemperatureType>::exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr )
{
    this->log("Heat","exportResults", "start");
    this->timerTool("PostProcessing").start();

    if ( M_exporter && M_exporter->exporterGeometry() == EXPORTER_GEOMETRY_CHANGE ) // TODO mv this code
        M_exporter->defaultTimeSet()->setMesh( this->mesh() );
    this->executePostProcessExports( M_exporter, time, mfields, symbolsExpr, exportsExpr );
    this->executePostProcessMeasures( time, mfields, symbolsExpr, this->modelMeasuresQuantities() );
    this->executePostProcessSave( (this->isStationary())? invalid_uint32_type_value : M_bdfTemperature->iteration(), mfields );

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
template <typename ModelFieldsType, typename SymbolsExpr, typename ModelMeasuresQuantitiesType>
void
Heat<ConvexType,BasisTemperatureType>::executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ModelMeasuresQuantitiesType const& mquantities )
{
    auto const& t = mfields.field( FieldTag::temperature(this), "temperature" );

    // compute measures
    for ( auto const& [ppName,ppFlux] : M_postProcessMeasuresNormalHeatFlux )
    {
        auto range = markedfaces(this->mesh(),ppFlux.markers() );
        auto heatFluxExpr = this->normalHeatFluxExpr( t, ppFlux.isOutward(), symbolsExpr );
        auto heatFluxExprUsed = evalOnFaces( std::move(heatFluxExpr),ppFlux.requiresMarkersConnection(),ppFlux.internalFacesEvalutationType() );
        double heatFlux = integrate(_range=range,
                                    _expr=heatFluxExprUsed ).evaluate()(0,0);
        this->postProcessMeasures().setValue("Normal_Heat_Flux_"+ppName,heatFlux);
    }

    // execute common post process and save measures
    super_type::executePostProcessMeasures( time, this->mesh(), M_rangeMeshElements, symbolsExpr, mfields, mquantities );
}

} // namespace FeelModels
} // namespace Feel

#include <feel/feelmodels/heat/heatassembly.hpp>
#include <feel/feelmodels/heat/heatupdatestabilizationgls.hpp>


#endif /* FEELPP_TOOLBOXES_HEAT_HPP */
