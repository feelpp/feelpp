/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/heatfluid/heatfluid.hpp>

namespace Feel
{
namespace FeelModels
{

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("HeatFluid","updateLinearPDE", "start"+sc);

    auto mctx = this->modelContext( vecCurrentPicardSolution, this->heatModel(), this->fluidModel() );
    if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    {
        auto previousSol = data.vectorInfo( "time-stepping.previous-solution");
        //auto mctxPrevious = this->modelContextNoTrialSymbolsExpr( previousSol, this->heatModel()->rowStartInVector(), this->fluidModel()->rowStartInVector() );
        auto mctxPrevious = this->modelContextNoTrialSymbolsExpr(
            { { "solution", std::make_tuple( previousSol, this->heatModel()->startBlockSpaceIndexVector() ) } },
            {
                { "solution", std::make_tuple( previousSol, this->fluidModel()->startBlockSpaceIndexVector()) },
                { "velocity_extrapolated", std::make_tuple( this->fluidModel()->vectorPreviousVelocityExtrapolated(), 0 ) }
            } );
        mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    }


    M_heatModel->updateLinearPDE( data, mctx );
    M_fluidModel->updateLinearPDE( data, mctx );

    if ( buildNonCstPart && M_useNaturalConvection )
    {
        double timeSteppingScaling_fluid = 1.;
        if ( !M_fluidModel->isStationaryModel() )
            timeSteppingScaling_fluid = data.doubleInfo( prefixvm(M_fluidModel->prefix(),"time-stepping.scaling") );
        double timeSteppingScaling_heat = 1.;
        if ( !M_heatModel->isStationary() )
            timeSteppingScaling_heat = data.doubleInfo( prefixvm(M_heatModel->prefix(),"time-stepping.scaling") );

        auto const& se = mctx.symbolsExpr();
        auto const& t = mctx.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );
        auto const& u = mctx.field( fluid_model_type::FieldTag::velocity(this->fluidModel().get()), "velocity" );

        auto XhT = this->heatModel()->spaceTemperature();
        auto XhV = this->fluidModel()->functionSpaceVelocity();
        size_type blockIndexTemperature = this->heatModel()->startBlockSpaceIndexVector() + this->heatModel()->startSubBlockSpaceIndex("temperature");
        size_type blockIndexVelocity = this->fluidModel()->startBlockSpaceIndexVector() + this->fluidModel()->startSubBlockSpaceIndex("velocity");

        auto bfVT = form2( _test=XhV,_trial=XhT,_matrix=A,
                           _rowstart=M_fluidModel->rowStartInMatrix(),
                           _colstart=M_heatModel->colStartInMatrix() );
        auto lfV = form1( _test=XhV, _vector=F,
                          _rowstart=M_fluidModel->rowStartInVector() );

        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        {
            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
            {
                auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                auto const& matProps = this->materialsProperties()->materialProperties( matName );
                auto const& rho = this->materialsProperties()->rho( matName );
                //auto const& rhoHeatCapacity = this->materialsProperties()->rhoHeatCapacity( matName );
                auto const& thermalExpansion = this->materialsProperties()->thermalExpansion( matName );

                //auto rhoHeatCapacityExpr = expr( rhoHeatCapacity.expr(), se );
                auto rhoExpr = expr( rho.template expr<1,1>(), se );
                auto betaExpr = expr( thermalExpansion.template expr<1,1>(), se );
                double T0 = M_BoussinesqRefTemperature;

                if ( M_useNaturalConvection )
                {
                    // TODO cst part if coeff are cst
                    bfVT +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling_fluid*rhoExpr*betaExpr*idt(t)*inner(M_gravityForce,id(u)),
                                   _geomap=this->geomap() );

                    lfV +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling_fluid*rhoExpr*(betaExpr*T0)*inner(M_gravityForce,id(u)),
                                   _geomap=this->geomap() );

                    if ( this->fluidModel()->stabilizationGLS() )
                    {
                        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>( physicData->subphysicFromType( M_fluidModel->physicType() ) );
                        auto exprAdded_lhs = -timeSteppingScaling_fluid*rhoExpr*betaExpr*idv(t)*M_gravityForce;
                        auto exprAdded_rhs = timeSteppingScaling_fluid*rhoExpr*betaExpr*T0*M_gravityForce;
                        auto exprAddedInGLSResidual_lhs = exprOptionalConcat< std::decay_t<decltype( exprAdded_lhs )> >();
                        using expr_rhs_theta_previous_type = std::decay_t<decltype( (-(1.0-timeSteppingScaling_fluid))*rhoExpr*betaExpr*(idv(t)-T0)*M_gravityForce )>;
                        auto exprAddedInGLSResidual_rhs = exprOptionalConcat< std::decay_t<decltype( exprAdded_rhs )>, expr_rhs_theta_previous_type >();
                        exprAddedInGLSResidual_lhs.expression().add( exprAdded_lhs );
                        exprAddedInGLSResidual_rhs.expression().add( exprAdded_rhs );
                        if ( !this->fluidModel()->isStationaryModel() && this->fluidModel()->timeStepping() == "Theta" )
                        {
                            auto const& mctxPrevious = mctx.additionalContext( "time-stepping.previous-model-context" );
                            auto const& sePrevious = mctxPrevious.symbolsExpr();
                            auto const& tOld = mctxPrevious.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );
                            auto rhoExprPrevious = expr( rho.template expr<1,1>(), sePrevious );
                            auto betaExprPrevious = expr( thermalExpansion.template expr<1,1>(), sePrevious );
                            auto exprAddedInRhsOld = (-(1.0-timeSteppingScaling_fluid))*rhoExprPrevious*betaExprPrevious*(idv(tOld)-T0)*M_gravityForce;
                            if ( data.hasParameterValuesInfo( "time-stepping.previous-parameter-values" ) )
                                exprAddedInRhsOld.setParameterValues( data.parameterValuesInfo( "time-stepping.previous-parameter-values" ) );
                            exprAddedInGLSResidual_rhs.expression().add( exprAddedInRhsOld );
                        }
                        M_fluidModel->updateLinearPDEStabilizationGLS( data, mctx, *physicFluidData, matProps, range,
                                                                       hana::make_tuple(exprAddedInGLSResidual_rhs),
                                                                       hana::make_tuple(std::make_tuple(XhT,blockIndexTemperature,exprAddedInGLSResidual_lhs) ) );
                    }
                }
            }
        }

    }

    this->log("HeatFluid","updateLinearPDE", "finish");
}




HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( vecCurrentSolution, this->heatModel(), this->fluidModel() );
    M_heatModel->updateLinearPDEDofElimination( data,mctx );
    M_fluidModel->updateLinearPDEDofElimination( data,mctx );
}



HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateLinear_Heat( DataUpdateLinear & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( vecCurrentSolution, this->heatModel()->startBlockSpaceIndexVector(),
                                    this->fluidModel()->algebraicBlockVectorSolution()->vectorMonolithic(), 0 );
    if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    {
        auto previousSolHeat = data.vectorInfo( "time-stepping.previous-solution");
        auto mctxPrevious = this->modelContextNoTrialSymbolsExpr(
            { { "solution", std::make_tuple( previousSolHeat, this->heatModel()->startBlockSpaceIndexVector() ) } },
            {
                { "solution", std::make_tuple( M_timeStepThetaSchemePreviousSolution, this->startSubBlockSpaceIndex("fluid") ) },
                { "velocity_extrapolated", std::make_tuple( this->fluidModel()->vectorPreviousVelocityExtrapolated(), 0 ) }
            } );

        mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    }

    M_heatModel->updateLinearPDE( data,mctx );
}





HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateLinear_Fluid( DataUpdateLinear & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("HeatFluid","updateLinear_Fluid", "start"+sc);

    auto mctx = this->modelContext( this->heatModel()->algebraicBlockVectorSolution()->vectorMonolithic(), 0,
                                    vecCurrentSolution, this->fluidModel()->startBlockSpaceIndexVector() );
    if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    {
        auto previousSolFluid = data.vectorInfo( "time-stepping.previous-solution");
        auto mctxPrevious = this->modelContextNoTrialSymbolsExpr(
            { { "solution", std::make_tuple( M_timeStepThetaSchemePreviousSolution, this->startSubBlockSpaceIndex("heat") ) } },
            {
                { "solution", std::make_tuple( previousSolFluid, this->fluidModel()->startBlockSpaceIndexVector()) },
                { "velocity_extrapolated", std::make_tuple( M_fluidModel->vectorPreviousVelocityExtrapolated(), 0 ) }
            } );
        mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    }

    M_fluidModel->updateLinearPDE( data,mctx );

    if ( !buildNonCstPart )
        return;
    if ( !M_useNaturalConvection )
        return;

    double timeSteppingScaling = 1.;
    if ( !M_fluidModel->isStationaryModel() )
        timeSteppingScaling = data.doubleInfo( prefixvm(M_fluidModel->prefix(),"time-stepping.scaling") );

    auto const& se = mctx.symbolsExpr();
    auto const& t = mctx.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );
    auto const& u = mctx.field( fluid_model_type::FieldTag::velocity(this->fluidModel().get()), "velocity" );
    //auto const& p = mctx.field( fluid_model_type::FieldTag::pressure(this->fluidModel().get()), "pressure" );

    auto XhV = M_fluidModel->functionSpaceVelocity();
    auto XhT = M_heatModel->spaceTemperature();

    auto mylfV = form1( _test=XhV, _vector=F,
                        _rowstart=M_fluidModel->rowStartInVector() );

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& matProps = this->materialsProperties()->materialProperties( matName );
            auto const& rhoHeatCapacity = this->materialsProperties()->rhoHeatCapacity( matName );

            auto const& rho = this->materialsProperties()->rho( matName );
            auto const& thermalExpansion = this->materialsProperties()->thermalExpansion( matName );
            auto const& rhoExpr = expr( rho.template expr<1,1>(), se );
            auto const& betaExpr = expr( thermalExpansion.template expr<1,1>(), se );
            double T0 = M_BoussinesqRefTemperature;

            mylfV +=
                integrate( _range=range,
                           _expr= -timeSteppingScaling*rhoExpr*betaExpr*(idv(t)-T0)*inner(M_gravityForce,id(u)),
                           _geomap=this->geomap() );

            if ( this->fluidModel()->stabilizationGLS() )
            {
                auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>( physicData->subphysicFromType( M_fluidModel->physicType() ) );
                auto exprAdded = (-timeSteppingScaling)*rhoExpr*betaExpr*(idv(t)-T0)*M_gravityForce;
                auto exprAddedInGLSResidual = exprOptionalConcat< std::decay_t<decltype( exprAdded )> >();
                exprAddedInGLSResidual.expression().add( exprAdded );
                if ( !this->fluidModel()->isStationaryModel() && this->fluidModel()->timeStepping() ==  "Theta" )
                {
                    auto const& mctxPrevious = mctx.additionalContext( "time-stepping.previous-model-context" );
                    auto const& sePrevious = mctxPrevious.symbolsExpr();
                    auto const& tOld = mctxPrevious.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );
                    auto rhoExprPrevious = expr( rho.template expr<1,1>(), sePrevious );
                    auto betaExprPrevious = expr( thermalExpansion.template expr<1,1>(), sePrevious );
                    auto exprAddedInRhsOld = (-(1.0-timeSteppingScaling))*rhoExprPrevious*betaExprPrevious*(idv(tOld)-T0)*M_gravityForce;
                    if ( data.hasParameterValuesInfo( "time-stepping.previous-parameter-values" ) )
                        exprAddedInRhsOld.setParameterValues( data.parameterValuesInfo( "time-stepping.previous-parameter-values" ) );
                    exprAddedInGLSResidual.expression().add( exprAddedInRhsOld );
                }
                M_fluidModel->updateLinearPDEStabilizationGLS( data, mctx, *physicFluidData, matProps, range, hana::make_tuple(exprAddedInGLSResidual) );
            }
        }
    }

    this->log("HeatFluid","updateLinear_Fluid", "finish"+sc);
}


} // namespace FeelModels
} // namespace Feel
