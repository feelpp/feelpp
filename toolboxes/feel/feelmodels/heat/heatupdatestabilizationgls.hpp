/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_HEAT_UPDATESTABILIZATIONGLS_HPP
#define FEELPP_TOOLBOXES_HEAT_UPDATESTABILIZATIONGLS_HPP 1

#include <feel/feelmodels/modelvf/stabilizationglsparameter.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>
//#include <feel/feeldiscr/pchv.hpp>
//#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelmesh/intersect.hpp>

#include <feel/feelmodels/modelvf/exproperations.hpp>

namespace Feel
{
namespace FeelModels
{

namespace Heat_detail
{

template<typename HeatType,typename ModelContextType, typename ... ExprAddedType>
auto
exprResidualImpl( HeatType const& heat, ModelPhysicHeat<HeatType::nDim> const& physicHeatData, MaterialProperties const& matProps,
                  ModelContextType const& mctx, double timeSteppingScaling, bool timeSteppingEvaluateResidualWithoutTimeDerivative = false )
{
    auto const& u = mctx.field( HeatType::FieldTag::temperature(&heat), "temperature" );
    auto const& se = mctx.symbolsExpr();
    std::string const& matName = matProps.materialName();

    using expr_density_type = std::decay_t<decltype( expr( matProps.property("density").template expr<1,1>(), se ) )>;
    using expr_heat_capacity_type = std::decay_t<decltype( expr( matProps.property("specific-heat-capacity").template expr<1,1>(), se ) )>;
    using expr_velocity_convection_type = std::decay_t<decltype( std::declval<ModelPhysicHeat<HeatType::nDim>>().convection().expr( se ) )>;
    using expr_thermal_conductivity_scalar_type = std::decay_t<decltype( expr( matProps.property("thermal-conductivity").template expr<1,1>(), se ) )>;
    using expr_thermal_conductivity_matrix_type = std::decay_t<decltype( expr( matProps.property("thermal-conductivity").template expr<HeatType::nDim,HeatType::nDim>(), se ) )>;
    using expr_grad_thermal_conductivity_scalar_type = std::decay_t<decltype( grad<HeatType::nDim>( expr_thermal_conductivity_scalar_type{} ) )>;
    using expr_heatsource_type = std::decay_t<decltype( std::declval<ModelPhysicHeat<HeatType::nDim>>().heatSources().front().expr( se ) )>;

    using expr_residual_time_derivative_type = std::decay_t<decltype( expr_density_type{}*expr_heat_capacity_type{}*(heat.timeStepBdfTemperature()->polyDerivCoefficient(0)*idv(u)-idv(heat.timeStepBdfTemperature()->polyDeriv())) )>;
    using expr_residual_convection_type = std::decay_t<decltype( timeSteppingScaling*expr_density_type{}*expr_heat_capacity_type{}*gradv(u)*expr_velocity_convection_type{} )>;
    using expr_residual_diffusion_scalar_type = std::decay_t<decltype( -timeSteppingScaling*expr_thermal_conductivity_scalar_type{}*laplacianv(u) )>;
    using expr_residual_diffusion_matrix_type = std::decay_t<decltype( -timeSteppingScaling*inner(expr_thermal_conductivity_matrix_type{},hessv(u)) )>;
    using expr_residual_diffusion_variable_conductivity_scalar_type = std::decay_t<decltype( -timeSteppingScaling*inner(expr_grad_thermal_conductivity_scalar_type{},gradv(u)) )>;
    using expr_residual_heatsource_type = std::decay_t<decltype( -timeSteppingScaling*expr_heatsource_type{} )>;

    auto residual_full = exprOptionalConcat<expr_residual_time_derivative_type,expr_residual_convection_type,
                                            expr_residual_diffusion_scalar_type,expr_residual_diffusion_matrix_type,
                                            expr_residual_diffusion_variable_conductivity_scalar_type,
                                            expr_residual_heatsource_type,
                                            ExprAddedType...>();


    auto thermalConductivityExprBase = matProps.property("thermal-conductivity");
    bool thermalConductivityIsMatrix = thermalConductivityExprBase.template hasExpr<HeatType::nDim,HeatType::nDim>();
    bool thermalConductivityDependsOnCoordinatesInSpace = thermalConductivityExprBase.template hasSymbolDependencyOnCoordinatesInSpace<HeatType::nDim>( se ) && heat.stabilizationGLS_checkConductivityDependencyOnCoordinates();
    if (  HeatType::nOrderTemperature > 1 || thermalConductivityDependsOnCoordinatesInSpace )
    {
        if ( thermalConductivityIsMatrix )
        {
            auto thermalConductivityExpr = expr( thermalConductivityExprBase.template expr<HeatType::nDim,HeatType::nDim>(), se );
            if constexpr ( HeatType::nOrderTemperature > 1 )
            {
                residual_full.expression().add( -timeSteppingScaling*inner(thermalConductivityExpr,hessv(u)) );
            }
            if ( thermalConductivityDependsOnCoordinatesInSpace )
            {
                CHECK( false ) << "thermalConductivity depends on CoordinatesInSpace not implemented with matrix case";
            }
        }
        else
        {
            auto thermalConductivityExpr = expr( thermalConductivityExprBase.template expr<1,1>(), se );
            if constexpr ( HeatType::nOrderTemperature > 1 )
            {
                residual_full.expression().add( -timeSteppingScaling*thermalConductivityExpr*laplacianv(u) );
            }
            if ( thermalConductivityDependsOnCoordinatesInSpace )
            {
                auto gradThermalConductivityExpr = grad<HeatType::nDim>( thermalConductivityExpr, "", heat.worldComm(),heat.repository().expr() );
                residual_full.expression().add( -timeSteppingScaling*inner(gradThermalConductivityExpr,gradv(u)) );
            }
        }
    }

    if ( !heat.isStationary() && !timeSteppingEvaluateResidualWithoutTimeDerivative )
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        auto heatCapacityExpr = expr( matProps.property("specific-heat-capacity").template expr<1,1>(), se );
        residual_full.expression().add( densityExpr*heatCapacityExpr*(heat.timeStepBdfTemperature()->polyDerivCoefficient(0)*idv(u) - idv(heat.timeStepBdfTemperature()->polyDeriv()) ) );
    }

    if ( physicHeatData.hasConvectionEnabled() )
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        auto heatCapacityExpr = expr( matProps.property("specific-heat-capacity").template expr<1,1>(), se );
        auto velConvExpr = physicHeatData.convection().expr( se );
        residual_full.expression().add( timeSteppingScaling*densityExpr*heatCapacityExpr*gradv(u)*velConvExpr );
    }

    for ( auto const& heatSource : physicHeatData.heatSources() )
    {
        auto hsExpr = heatSource.expr( se );
        residual_full.expression().add( -timeSteppingScaling*hsExpr );
    }


    return residual_full;
}

template<typename HeatType,typename ModelContextType>
auto
exprResidual( HeatType const& heat, ModelPhysicHeat<HeatType::nDim> const& physicHeatData, MaterialProperties const& matProps,
              ModelContextType const& mctx, double timeSteppingScaling, bool timeSteppingEvaluateResidualWithoutTimeDerivative = false )
{
    return exprResidualImpl( heat,physicHeatData,matProps,mctx,timeSteppingScaling,timeSteppingEvaluateResidualWithoutTimeDerivative );
}
template<typename ExprAddedType,typename HeatType,typename ModelContextType>
auto
exprResidual( HeatType const& heat, ModelPhysicHeat<HeatType::nDim> const& physicHeatData, MaterialProperties const& matProps,
              ModelContextType const& mctx, double timeSteppingScaling, bool timeSteppingEvaluateResidualWithoutTimeDerivative = false )
{
    return exprResidualImpl<HeatType,ModelContextType,ExprAddedType>( heat,physicHeatData,matProps,mctx,timeSteppingScaling,timeSteppingEvaluateResidualWithoutTimeDerivative );
}

} //  namespace Heat_detail




template< typename ConvexType, typename BasisTemperatureType >
template <typename ModelContextType,typename RangeType>
void
Heat<ConvexType,BasisTemperatureType>::updateJacobianStabilizationGLS( DataUpdateJacobian & data, ModelContextType const& mctx,
                                                                       ModelPhysicHeat<nDim> const& physicHeatData,
                                                                       MaterialProperties const& matProps, RangeType const& range ) const
{
    sparse_matrix_ptrtype& J = data.jacobian();

    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Heat","updateJacobianStabilizationGLS", "start"+sc);

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    std::string const& matName = matProps.materialName();
    auto mesh = this->mesh();
    auto XhT = this->spaceTemperature();
    auto const& u = mctx.field( FieldTag::temperature(this), "temperature" );
    auto const& v = u;
    auto const& se = mctx.symbolsExpr();

    auto bilinearFormTT = form2( _test=XhT,_trial=XhT,_matrix=J,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=this->rowStartInMatrix(),
                                 _colstart=this->colStartInMatrix() );

    int coeffNatureStabilization = ( this->stabilizationGLSType() == "supg" )? 0 : (this->stabilizationGLSType() == "gls")? 1 : -1;

    using expr_density_type = std::decay_t<decltype( expr( matProps.property("density").template expr<1,1>(), se ) )>;
    using expr_heat_capacity_type = std::decay_t<decltype( expr( matProps.property("specific-heat-capacity").template expr<1,1>(), se ) )>;
    using expr_velocity_convection_type = std::decay_t<decltype( std::declval<ModelPhysicHeat<nDim>>().convection().expr( se ) )>;

    using expr_coeff_null_type = std::decay_t<decltype( cst(0.) )>;
    using expr_thermal_conductivity_scalar_type = std::decay_t<decltype( expr( matProps.property("thermal-conductivity").template expr<1,1>(), se ) )>;
    using expr_thermal_conductivity_matrix_type = std::decay_t<decltype( expr( matProps.property("thermal-conductivity").template expr<nDim,nDim>(), se ) )>;
    using expr_grad_thermal_conductivity_scalar_type = std::decay_t<decltype( grad<nDim>( expr_thermal_conductivity_scalar_type{} ) )>;

    using expr_coeff_convection_part1_type = std::decay_t<decltype( expr_density_type{}*expr_heat_capacity_type{}*expr_velocity_convection_type{} )>;
    using expr_coeff_convection_part2_type = std::decay_t<decltype( -trans(expr_grad_thermal_conductivity_scalar_type{}) )>;
    auto exprCoeffConvection = exprOptionalConcat<expr_coeff_convection_part1_type,expr_coeff_convection_part2_type>();
    using expr_coeff_convection_type = std::decay_t<decltype( exprCoeffConvection )>;

    auto tauExprScalarDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_thermal_conductivity_scalar_type,expr_coeff_null_type>( *this->stabilizationGLSParameter() );
    auto tauExprMatrixDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_thermal_conductivity_matrix_type,expr_coeff_null_type>( *this->stabilizationGLSParameter() );

    using expr_residual_lhs_time_derivative_type = std::decay_t<decltype( expr_density_type{}*expr_heat_capacity_type{}*this->timeStepBdfTemperature()->polyDerivCoefficient(0)*idt(u) )>;
    using expr_residual_lhs_convection_type = std::decay_t<decltype( timeSteppingScaling*expr_density_type{}*expr_heat_capacity_type{}*gradt(u)*expr_velocity_convection_type{} )>;
    using expr_residual_lhs_diffusion_scalar_type = std::decay_t<decltype( -timeSteppingScaling*expr_thermal_conductivity_scalar_type{}*laplaciant(u) )>;
    using expr_residual_lhs_diffusion_matrix_type = std::decay_t<decltype( -timeSteppingScaling*inner(expr_thermal_conductivity_matrix_type{},hesst(u)) )>;
    using expr_residual_lhs_diffusion_variable_conductivity_scalar_type = std::decay_t<decltype( -timeSteppingScaling*inner(expr_grad_thermal_conductivity_scalar_type{},gradt(u)) )>;
    auto residual_lhs = exprOptionalConcat<expr_residual_lhs_time_derivative_type,expr_residual_lhs_convection_type,expr_residual_lhs_diffusion_scalar_type,expr_residual_lhs_diffusion_matrix_type,
                                           expr_residual_lhs_diffusion_variable_conductivity_scalar_type>();

    // test functions
    using expr_test_convection_type = std::decay_t<decltype( expr_density_type{}*expr_heat_capacity_type{}*grad(u)*expr_velocity_convection_type{} )>;
    using expr_test_diffusion_scalar_type = std::decay_t<decltype( -coeffNatureStabilization*expr_thermal_conductivity_scalar_type{}*laplacian(u) )>;
    using expr_test_diffusion_matrix_type = std::decay_t<decltype( -coeffNatureStabilization*inner(expr_thermal_conductivity_matrix_type{},hess(u)) )>;
    using expr_test_diffusion_variable_conductivity_scalar_type = std::decay_t<decltype( -/*coeffNatureStabilization**/inner(expr_grad_thermal_conductivity_scalar_type{},grad(u)) )>;
    auto stab_test = exprOptionalConcat<expr_test_convection_type,expr_test_diffusion_scalar_type,expr_test_diffusion_matrix_type,
                                        expr_test_diffusion_variable_conductivity_scalar_type>();


    auto thermalConductivityExprBase = matProps.property("thermal-conductivity");
    bool thermalConductivityIsMatrix = thermalConductivityExprBase.template hasExpr<nDim,nDim>();
    bool thermalConductivityDependsOnCoordinatesInSpace = thermalConductivityExprBase.template hasSymbolDependencyOnCoordinatesInSpace<nDim>( se ) && this->stabilizationGLS_checkConductivityDependencyOnCoordinates();
    if ( thermalConductivityIsMatrix )
    {
        auto thermalConductivityExpr = expr( thermalConductivityExprBase.template expr<nDim,nDim>(), se );
        if constexpr ( nOrderTemperature > 1 )
        {
            residual_lhs.expression().add( -timeSteppingScaling*inner(thermalConductivityExpr,hesst(u)) );
            if ( coeffNatureStabilization != 0 )
                stab_test.expression().add( -coeffNatureStabilization*inner(thermalConductivityExpr,hess(u)) );
        }
        if ( thermalConductivityDependsOnCoordinatesInSpace )
        {
            CHECK( false ) << "thermalConductivity depends on CoordinatesInSpace not implemented with matrix case";
        }
        tauExprMatrixDiffusion.expression().setDiffusion( thermalConductivityExpr );
    }
    else
    {
        auto thermalConductivityExpr = expr( thermalConductivityExprBase.template expr<1,1>(), se );
        if constexpr ( nOrderTemperature > 1 )
        {
            residual_lhs.expression().add( -timeSteppingScaling*thermalConductivityExpr*laplaciant(u) );
            if ( coeffNatureStabilization != 0 )
                stab_test.expression().add( -coeffNatureStabilization*thermalConductivityExpr*laplacian(u) );
        }
        if ( thermalConductivityDependsOnCoordinatesInSpace )
        {
            auto gradThermalConductivityExpr = grad<nDim>( thermalConductivityExpr, "", this->worldComm(),this->repository().expr() );
            residual_lhs.expression().add( -timeSteppingScaling*inner(gradThermalConductivityExpr,gradt(u)) );
            //if ( coeffNatureStabilization != 0 )
            stab_test.expression().add( -/*coeffNatureStabilization**/inner(gradThermalConductivityExpr,grad(u)) );
            exprCoeffConvection.expression().add( -trans(gradThermalConductivityExpr) );
        }
        tauExprScalarDiffusion.expression().setDiffusion( thermalConductivityExpr );
    }

    if ( !this->isStationary() )
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        auto heatCapacityExpr = expr( matProps.property("specific-heat-capacity").template expr<1,1>(), se );
        residual_lhs.expression().add( densityExpr*heatCapacityExpr*this->timeStepBdfTemperature()->polyDerivCoefficient(0)*idt(u) );
    }

    if ( physicHeatData.hasConvectionEnabled() )
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        auto heatCapacityExpr = expr( matProps.property("specific-heat-capacity").template expr<1,1>(), se );
        auto velConvExpr = physicHeatData.convection().expr( se );
        residual_lhs.expression().add( timeSteppingScaling*densityExpr*heatCapacityExpr*gradt(u)*velConvExpr );
        stab_test.expression().add( densityExpr*heatCapacityExpr*grad(u)*velConvExpr );
        exprCoeffConvection.expression().add( densityExpr*heatCapacityExpr*velConvExpr );
    }

    if ( exprCoeffConvection.expression().hasExpr() )
    {
        if ( thermalConductivityIsMatrix )
            tauExprMatrixDiffusion.expression().setConvection( exprCoeffConvection );
        else
            tauExprScalarDiffusion.expression().setConvection( exprCoeffConvection );
    }

    auto tauFieldPtr = this->stabilizationGLSParameter()->fieldTauPtr();
#if 0 // already done with updateResidual
    if ( thermalConductivityIsMatrix )
        tauFieldPtr->on(_range=range,_expr=tauExprMatrixDiffusion);
    else
        tauFieldPtr->on(_range=range,_expr=tauExprScalarDiffusion);
#endif
    auto tau = idv(tauFieldPtr);

    bilinearFormTT +=
        integrate( _range=range,
                   _expr=tau*residual_lhs*stab_test,
                   _geomap=this->geomap() );

    this->log("Heat","updateJacobianStabilizationGLS", "finish"+sc);
}

template< typename ConvexType, typename BasisTemperatureType >
template <typename ModelContextType,typename RangeType,typename... ExprAddedType>
void
Heat<ConvexType,BasisTemperatureType>::updateResidualStabilizationGLS( DataUpdateResidual & data, ModelContextType const& mctx,
                                                                       ModelPhysicHeat<nDim> const& physicHeatData,
                                                                       MaterialProperties const& matProps, RangeType const& range,
                                                                       const ExprAddedType&... exprsAddedInResidual ) const
{
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;

    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false;
    if ( !this->isStationary() )
    {
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm(this->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );
    }
    if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
        return;

    std::string const& matName = matProps.materialName();
    auto XhT = this->spaceTemperature();
    auto const& u = mctx.field( FieldTag::temperature(this), "temperature" );
    auto const& v = u;
    auto const& se = mctx.symbolsExpr();

    auto myLinearFormT = form1( _test=XhT, _vector=R,
                                _rowstart=this->rowStartInVector() );


    int coeffNatureStabilization = ( this->stabilizationGLSType() == "supg" )? 0 : (this->stabilizationGLSType() == "gls")? 1 : -1;

    using expr_density_type = std::decay_t<decltype( expr( matProps.property("density").template expr<1,1>(), se ) )>;
    using expr_heat_capacity_type = std::decay_t<decltype( expr( matProps.property("specific-heat-capacity").template expr<1,1>(), se ) )>;
    using expr_velocity_convection_type = std::decay_t<decltype( std::declval<ModelPhysicHeat<nDim>>().convection().expr( se ) )>;

    using expr_coeff_null_type = std::decay_t<decltype( cst(0.) )>;
    using expr_thermal_conductivity_scalar_type = std::decay_t<decltype( expr( matProps.property("thermal-conductivity").template expr<1,1>(), se ) )>;
    using expr_thermal_conductivity_matrix_type = std::decay_t<decltype( expr( matProps.property("thermal-conductivity").template expr<nDim,nDim>(), se ) )>;
    using expr_grad_thermal_conductivity_scalar_type = std::decay_t<decltype( grad<nDim>( expr_thermal_conductivity_scalar_type{} ) )>;
    using expr_coeff_convection_part1_type = std::decay_t<decltype( expr_density_type{}*expr_heat_capacity_type{}*expr_velocity_convection_type{} )>;
    using expr_coeff_convection_part2_type = std::decay_t<decltype( -trans(expr_grad_thermal_conductivity_scalar_type{}) )>;
    auto exprCoeffConvection = exprOptionalConcat<expr_coeff_convection_part1_type,expr_coeff_convection_part2_type>();
    using expr_coeff_convection_type = std::decay_t<decltype( exprCoeffConvection )>;

    auto tauExprScalarDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_thermal_conductivity_scalar_type,expr_coeff_null_type>( *this->stabilizationGLSParameter() );
    auto tauExprMatrixDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_thermal_conductivity_matrix_type,expr_coeff_null_type>( *this->stabilizationGLSParameter() );

    // test functions
    using expr_test_convection_type = std::decay_t<decltype( expr_density_type{}*expr_heat_capacity_type{}*grad(u)*expr_velocity_convection_type{} )>;
    using expr_test_diffusion_scalar_type = std::decay_t<decltype( -coeffNatureStabilization*expr_thermal_conductivity_scalar_type{}*laplacian(u) )>;
    using expr_test_diffusion_matrix_type = std::decay_t<decltype( -coeffNatureStabilization*inner(expr_thermal_conductivity_matrix_type{},hess(u)) )>;
    using expr_test_diffusion_variable_conductivity_scalar_type = std::decay_t<decltype( -/*coeffNatureStabilization**/inner(expr_grad_thermal_conductivity_scalar_type{},grad(u)) )>;
    auto stab_test = exprOptionalConcat<expr_test_convection_type,expr_test_diffusion_scalar_type,expr_test_diffusion_matrix_type,
                                        expr_test_diffusion_variable_conductivity_scalar_type>();


    using expr_previous_residual_type = std::decay_t<decltype( Heat_detail::exprResidual( *this, physicHeatData, matProps, mctx, 1.0-timeSteppingScaling, true ) )>;
    auto residual_base = Heat_detail::exprResidual<expr_previous_residual_type>( *this, physicHeatData, matProps, mctx, timeSteppingScaling, timeSteppingEvaluateResidualWithoutTimeDerivative );
    if ( !this->isStationary() && this->timeStepping() == "Theta" )
    {
        auto const& mctxPrevious = mctx.additionalContext( "time-stepping.previous-model-context" );
        auto previousResidualExpr = Heat_detail::exprResidual( *this, physicHeatData, matProps, mctxPrevious, 1.0-timeSteppingScaling, true );
        if ( data.hasParameterValuesInfo( "time-stepping.previous-parameter-values" ) )
            previousResidualExpr.setParameterValues( data.parameterValuesInfo( "time-stepping.previous-parameter-values" ) );
        residual_base.expression().add( previousResidualExpr );
    }

    auto residual_full = Feel::FeelModels::vfdetail::addExpr( hana::make_tuple( residual_base, exprsAddedInResidual... ) );

    auto thermalConductivityExprBase = matProps.property("thermal-conductivity");
    bool thermalConductivityIsMatrix = thermalConductivityExprBase.template hasExpr<nDim,nDim>();
    bool thermalConductivityDependsOnCoordinatesInSpace = thermalConductivityExprBase.template hasSymbolDependencyOnCoordinatesInSpace<nDim>( se ) && this->stabilizationGLS_checkConductivityDependencyOnCoordinates();
    if ( thermalConductivityIsMatrix )
    {
        auto thermalConductivityExpr = expr( thermalConductivityExprBase.template expr<nDim,nDim>(), se );
        if constexpr ( nOrderTemperature > 1 )
        {
            if ( coeffNatureStabilization != 0 )
                stab_test.expression().add( -coeffNatureStabilization*inner(thermalConductivityExpr,hess(u)) );
        }
        if ( thermalConductivityDependsOnCoordinatesInSpace )
        {
            CHECK( false ) << "thermalConductivity depends on CoordinatesInSpace not implemented with matrix case";
        }
        tauExprMatrixDiffusion.expression().setDiffusion( thermalConductivityExpr );
    }
    else
    {
        auto thermalConductivityExpr = expr( thermalConductivityExprBase.template expr<1,1>(), se );
        if constexpr ( nOrderTemperature > 1 )
        {
            if ( coeffNatureStabilization != 0 )
                stab_test.expression().add( -coeffNatureStabilization*thermalConductivityExpr*laplacian(u) );
        }
        if ( thermalConductivityDependsOnCoordinatesInSpace )
        {
            //if ( coeffNatureStabilization != 0 )
            {
                auto gradThermalConductivityExpr = grad<nDim>( thermalConductivityExpr, "", this->worldComm(),this->repository().expr() );
                stab_test.expression().add( -/*coeffNatureStabilization**/inner(gradThermalConductivityExpr,grad(u)) );
                exprCoeffConvection.expression().add( -trans(gradThermalConductivityExpr) );
            }
        }
        tauExprScalarDiffusion.expression().setDiffusion( thermalConductivityExpr );
    }

    if ( physicHeatData.hasConvectionEnabled() )
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        auto heatCapacityExpr = expr( matProps.property("specific-heat-capacity").template expr<1,1>(), se );
        auto velConvExpr = physicHeatData.convection().expr( se );
        exprCoeffConvection.expression().add( densityExpr*heatCapacityExpr*velConvExpr );
        stab_test.expression().add( densityExpr*heatCapacityExpr*grad(u)*velConvExpr );
    }

    if ( exprCoeffConvection.expression().hasExpr() )
    {
        if ( thermalConductivityIsMatrix )
            tauExprMatrixDiffusion.expression().setConvection( exprCoeffConvection );
        else
            tauExprScalarDiffusion.expression().setConvection( exprCoeffConvection );
    }

    auto tauFieldPtr = this->stabilizationGLSParameter()->fieldTauPtr();
    if ( thermalConductivityIsMatrix )
        tauFieldPtr->on(_range=range,_expr=tauExprMatrixDiffusion);
    else
        tauFieldPtr->on(_range=range,_expr=tauExprScalarDiffusion);
    auto tau = idv(tauFieldPtr);

    myLinearFormT +=
        integrate( _range=range,
                   _expr=tau*residual_full*stab_test,
                   _geomap=this->geomap() );
}

template< typename ConvexType, typename BasisTemperatureType >
template <typename ModelContextType,typename RangeType>
void
Heat<ConvexType,BasisTemperatureType>::updateLinearPDEStabilizationGLS( DataUpdateLinear & data, ModelContextType const& mctx,
                                                                        ModelPhysicHeat<nDim> const& physicHeatData,
                                                                        MaterialProperties const& matProps, RangeType const& range ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    if ( buildCstPart )
        return;
    bool doAssemblyRhs = !data.hasInfo( "ignore-assembly.rhs" );
    bool doAssemblyLhs = !data.hasInfo( "ignore-assembly.lhs" );

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("Heat","updateLinearPDEStabilizationGLS", "start"+sc);

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->prefix(),"time-stepping.scaling") );

    std::string const& matName = matProps.materialName();
    auto mesh = this->mesh();
    auto XhT = this->spaceTemperature();
    auto const& u = this->fieldTemperature();
    auto const& v = this->fieldTemperature();
    auto const& se = mctx.symbolsExpr();

    auto bilinearFormTT = form2( _test=XhT,_trial=XhT,_matrix=A,
                                 _pattern=size_type(Pattern::COUPLED),
                                 _rowstart=this->rowStartInMatrix(),
                                 _colstart=this->colStartInMatrix() );
    auto myLinearFormT = form1( _test=XhT, _vector=F,
                                _rowstart=this->rowStartInVector() );

    int coeffNatureStabilization = ( this->stabilizationGLSType() == "supg" )? 0 : (this->stabilizationGLSType() == "gls")? 1 : -1;

    using expr_density_type = std::decay_t<decltype( expr( matProps.property("density").template expr<1,1>(), se ) )>;
    using expr_heat_capacity_type = std::decay_t<decltype( expr( matProps.property("specific-heat-capacity").template expr<1,1>(), se ) )>;
    using expr_velocity_convection_type = std::decay_t<decltype( std::declval<ModelPhysicHeat<nDim>>().convection().expr( se ) )>;
    using expr_heatsource_type = std::decay_t<decltype( std::declval<ModelPhysicHeat<nDim>>().heatSources().front().expr( se ) )>;

    using expr_coeff_null_type = std::decay_t<decltype( cst(0.) )>;
    using expr_thermal_conductivity_scalar_type = std::decay_t<decltype( expr( matProps.property("thermal-conductivity").template expr<1,1>(), se ) )>;
    using expr_thermal_conductivity_matrix_type = std::decay_t<decltype( expr( matProps.property("thermal-conductivity").template expr<nDim,nDim>(), se ) )>;
    using expr_grad_thermal_conductivity_scalar_type = std::decay_t<decltype( grad<nDim>( expr_thermal_conductivity_scalar_type{} ) )>;
    using expr_coeff_convection_part1_type = std::decay_t<decltype( expr_density_type{}*expr_heat_capacity_type{}*expr_velocity_convection_type{} )>;
    using expr_coeff_convection_part2_type = std::decay_t<decltype( -trans(expr_grad_thermal_conductivity_scalar_type{}) )>;
    auto exprCoeffConvection = exprOptionalConcat<expr_coeff_convection_part1_type,expr_coeff_convection_part2_type>();
    using expr_coeff_convection_type = std::decay_t<decltype( exprCoeffConvection )>;

    auto tauExprScalarDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_thermal_conductivity_scalar_type,expr_coeff_null_type>( *this->stabilizationGLSParameter() );
    auto tauExprMatrixDiffusion = Feel::FeelModels::stabilizationGLSParameterExpr<expr_coeff_convection_type,expr_thermal_conductivity_matrix_type,expr_coeff_null_type>( *this->stabilizationGLSParameter() );

    using expr_residual_lhs_time_derivative_type = std::decay_t<decltype( expr_density_type{}*expr_heat_capacity_type{}*this->timeStepBdfTemperature()->polyDerivCoefficient(0)*idt(u) )>;
    using expr_residual_lhs_convection_type = std::decay_t<decltype( timeSteppingScaling*expr_density_type{}*expr_heat_capacity_type{}*gradt(u)*expr_velocity_convection_type{} )>;
    using expr_residual_lhs_diffusion_scalar_type = std::decay_t<decltype( -timeSteppingScaling*expr_thermal_conductivity_scalar_type{}*laplaciant(u) )>;
    using expr_residual_lhs_diffusion_matrix_type = std::decay_t<decltype( -timeSteppingScaling*inner(expr_thermal_conductivity_matrix_type{},hesst(u)) )>;
    using expr_residual_lhs_diffusion_variable_conductivity_scalar_type = std::decay_t<decltype( -timeSteppingScaling*inner(expr_grad_thermal_conductivity_scalar_type{},gradt(u)) )>;

    using expr_residual_rhs_time_derivative_type = std::decay_t<decltype( expr_density_type{}*expr_heat_capacity_type{}*idv(this->timeStepBdfTemperature()->polyDeriv()) )>;
    using expr_residual_rhs_previous_residual_type = std::decay_t<decltype( -Heat_detail::exprResidual( *this, physicHeatData, matProps, mctx, 1.0-timeSteppingScaling, true ) )>;
    using expr_residual_rhs_heatsource_type = std::decay_t<decltype( timeSteppingScaling*expr_heatsource_type{} )>;

    auto residual_lhs = exprOptionalConcat<expr_residual_lhs_time_derivative_type,expr_residual_lhs_convection_type,expr_residual_lhs_diffusion_scalar_type,expr_residual_lhs_diffusion_matrix_type,
                                           expr_residual_lhs_diffusion_variable_conductivity_scalar_type>();
    auto residual_rhs = exprOptionalConcat<expr_residual_rhs_time_derivative_type,expr_residual_rhs_previous_residual_type,expr_residual_rhs_heatsource_type>();


    // test functions
    using expr_test_convection_type = std::decay_t<decltype( expr_density_type{}*expr_heat_capacity_type{}*grad(u)*expr_velocity_convection_type{} )>;
    using expr_test_diffusion_scalar_type = std::decay_t<decltype( -coeffNatureStabilization*expr_thermal_conductivity_scalar_type{}*laplacian(u) )>;
    using expr_test_diffusion_matrix_type = std::decay_t<decltype( -coeffNatureStabilization*inner(expr_thermal_conductivity_matrix_type{},hess(u)) )>;
    using expr_test_diffusion_variable_conductivity_scalar_type = std::decay_t<decltype( -/*coeffNatureStabilization* */inner(expr_grad_thermal_conductivity_scalar_type{},grad(u)) )>;
    auto stab_test = exprOptionalConcat<expr_test_convection_type,expr_test_diffusion_scalar_type,expr_test_diffusion_matrix_type,
                                        expr_test_diffusion_variable_conductivity_scalar_type>();



    auto thermalConductivityExprBase = matProps.property("thermal-conductivity");
    bool thermalConductivityIsMatrix = thermalConductivityExprBase.template hasExpr<nDim,nDim>();
    bool thermalConductivityDependsOnCoordinatesInSpace = thermalConductivityExprBase.template hasSymbolDependencyOnCoordinatesInSpace<nDim>( se ) && this->stabilizationGLS_checkConductivityDependencyOnCoordinates();
    if ( thermalConductivityIsMatrix )
    {
        auto thermalConductivityExpr = expr( thermalConductivityExprBase.template expr<nDim,nDim>(), se );
        if constexpr ( nOrderTemperature > 1 )
        {
            residual_lhs.expression().add( -timeSteppingScaling*inner(thermalConductivityExpr,hesst(u)) );
            if ( coeffNatureStabilization != 0 )
                stab_test.expression().add( -coeffNatureStabilization*inner(thermalConductivityExpr,hess(u)) );
        }
        if ( thermalConductivityDependsOnCoordinatesInSpace )
        {
            CHECK( false ) << "thermalConductivity depends on CoordinatesInSpace not implemented with matrix case";
        }
        tauExprMatrixDiffusion.expression().setDiffusion( thermalConductivityExpr );
    }
    else
    {
        auto thermalConductivityExpr = expr( thermalConductivityExprBase.template expr<1,1>(), se );
        if constexpr ( nOrderTemperature > 1 )
        {
            residual_lhs.expression().add( -timeSteppingScaling*thermalConductivityExpr*laplaciant(u) );
            if ( coeffNatureStabilization != 0 )
                stab_test.expression().add( -coeffNatureStabilization*thermalConductivityExpr*laplacian(u) );
        }
        if ( thermalConductivityDependsOnCoordinatesInSpace )
        {
            auto gradThermalConductivityExpr = grad<nDim>( thermalConductivityExpr, "", this->worldComm(),this->repository().expr() );
            residual_lhs.expression().add( -timeSteppingScaling*inner(gradThermalConductivityExpr,gradt(u)) );
            //if ( coeffNatureStabilization != 0 )
            stab_test.expression().add( -/*coeffNatureStabilization* */inner(gradThermalConductivityExpr,grad(u)) );
            exprCoeffConvection.expression().add( -trans(gradThermalConductivityExpr) );
        }
        tauExprScalarDiffusion.expression().setDiffusion( thermalConductivityExpr );
    }

    if ( !this->isStationary() )
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        auto heatCapacityExpr = expr( matProps.property("specific-heat-capacity").template expr<1,1>(), se );
        residual_lhs.expression().add( densityExpr*heatCapacityExpr*this->timeStepBdfTemperature()->polyDerivCoefficient(0)*idt(u) );
        residual_rhs.expression().add( densityExpr*heatCapacityExpr*idv(this->timeStepBdfTemperature()->polyDeriv()) );
    }

    if ( physicHeatData.hasConvectionEnabled() )
    {
        auto densityExpr = expr( matProps.property("density").template expr<1,1>(), se );
        auto heatCapacityExpr = expr( matProps.property("specific-heat-capacity").template expr<1,1>(), se );
        auto velConvExpr = physicHeatData.convection().expr( se );
        residual_lhs.expression().add( timeSteppingScaling*densityExpr*heatCapacityExpr*gradt(u)*velConvExpr );
        stab_test.expression().add( densityExpr*heatCapacityExpr*grad(u)*velConvExpr );
        exprCoeffConvection.expression().add( densityExpr*heatCapacityExpr*velConvExpr );
    }

    for ( auto const& heatSource : physicHeatData.heatSources() )
    {
        auto hsExpr = heatSource.expr( se );
        residual_rhs.expression().add( timeSteppingScaling*hsExpr );
    }


    if ( !this->isStationary() && this->timeStepping() == "Theta" )
    {
        auto const& mctxPrevious = mctx.additionalContext( "time-stepping.previous-model-context" );
        auto previousResidualExpr = -Heat_detail::exprResidual( *this, physicHeatData, matProps, mctxPrevious, 1.0-timeSteppingScaling, true );
        if ( data.hasParameterValuesInfo( "time-stepping.previous-parameter-values" ) )
            previousResidualExpr.setParameterValues( data.parameterValuesInfo( "time-stepping.previous-parameter-values" ) );
        residual_rhs.expression().add( previousResidualExpr );
    }

    if ( exprCoeffConvection.expression().hasExpr() )
    {
        if ( thermalConductivityIsMatrix )
            tauExprMatrixDiffusion.expression().setConvection( exprCoeffConvection );
        else
            tauExprScalarDiffusion.expression().setConvection( exprCoeffConvection );
    }

    auto tauFieldPtr = this->stabilizationGLSParameter()->fieldTauPtr();
    if ( thermalConductivityIsMatrix )
        tauFieldPtr->on(_range=range,_expr=tauExprMatrixDiffusion);
    else
        tauFieldPtr->on(_range=range,_expr=tauExprScalarDiffusion);
    auto tau = idv(tauFieldPtr);

    if ( doAssemblyLhs )
    {
        bilinearFormTT +=
            integrate( _range=range,
                       _expr=tau*residual_lhs*stab_test,
                       _geomap=this->geomap() );
    }

    if ( doAssemblyRhs && residual_rhs.expression().hasExpr() )
    {
        myLinearFormT +=
            integrate( _range=range,
                       _expr=tau*residual_rhs*stab_test,
                       _geomap=this->geomap() );
    }

    this->log("Heat","updateLinearPDEStabilizationGLS", "finish"+sc);
}


} // namespace Feel
} // namespace FeelModels

#endif
