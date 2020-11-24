/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/heatfluid/heatfluid.hpp>

namespace Feel
{
namespace FeelModels
{

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool useJacobianLinearTerms = data.useJacobianLinearTerms();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("HeatFluid","updateResidual", "start"+sc);

    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false, timeSteppingEvaluateResidualWithoutTimeDerivative_fluid = false, timeSteppingEvaluateResidualWithoutTimeDerivative_heat = false;
    if ( !M_fluidModel->isStationaryModel() || !M_heatModel->isStationary() )
    {
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm(this->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        if ( !M_fluidModel->isStationaryModel() )
            timeSteppingEvaluateResidualWithoutTimeDerivative_fluid = data.hasInfo( prefixvm(M_fluidModel->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        if ( !M_heatModel->isStationary() )
            timeSteppingEvaluateResidualWithoutTimeDerivative_heat = data.hasInfo( prefixvm(M_heatModel->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
    }

    bool doAssemblyHeat = !timeSteppingEvaluateResidualWithoutTimeDerivative || timeSteppingEvaluateResidualWithoutTimeDerivative_heat;
    bool doAssemblyFluid = !timeSteppingEvaluateResidualWithoutTimeDerivative || timeSteppingEvaluateResidualWithoutTimeDerivative_fluid;

    //auto mesh = this->mesh();

    auto mctx = this->modelContext( XVec, this->heatModel()->rowStartInVector(), this->fluidModel()->rowStartInVector() );
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

    auto const& se = mctx.symbolsExpr();
    auto const& t = mctx.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );
    auto const& u = mctx.field( fluid_model_type::FieldTag::velocity(this->fluidModel().get()), "velocity" );
    auto const& p = mctx.field( fluid_model_type::FieldTag::pressure(this->fluidModel().get()), "pressure" );


    auto XhT = this->heatModel()->spaceTemperature();
    auto XhV = this->fluidModel()->functionSpaceVelocity();
    auto XhP = this->fluidModel()->functionSpacePressure();
    size_type blockIndexTemperature = this->heatModel()->rowStartInVector() + this->heatModel()->startSubBlockSpaceIndex("temperature");
    size_type blockIndexVelocity = this->fluidModel()->rowStartInVector() + this->fluidModel()->startSubBlockSpaceIndex("velocity");
    size_type blockIndexPressure = this->fluidModel()->rowStartInVector() + this->fluidModel()->startSubBlockSpaceIndex("pressure");


    if ( doAssemblyHeat )
        M_heatModel->updateResidual( data, mctx );
    if ( doAssemblyFluid )
        M_fluidModel->updateResidual( data, mctx );

    double timeSteppingScaling_fluid = 1.;
    if ( !M_fluidModel->isStationaryModel() && doAssemblyFluid )
        timeSteppingScaling_fluid = data.doubleInfo( prefixvm(M_fluidModel->prefix(),"time-stepping.scaling") );
    double timeSteppingScaling_heat = 1.;
    if ( !M_heatModel->isStationary() && doAssemblyHeat )
        timeSteppingScaling_heat = data.doubleInfo( prefixvm(M_heatModel->prefix(),"time-stepping.scaling") );

    if ( buildNonCstPart )
    {
        auto mylfT = form1( _test=XhT, _vector=R,
                            _rowstart=M_heatModel->rowStartInVector()+0 );
        auto mylfV = form1( _test=XhV, _vector=R,
                            _rowstart=M_fluidModel->rowStartInVector()+0 );


        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        {
            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
            {
                auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                //auto const& rhoHeatCapacity = this->materialsProperties()->rhoHeatCapacity( matName );
                //auto rhoHeatCapacityExpr = rhoHeatCapacity.expr();

                auto const& matProps = this->materialsProperties()->materialProperties( matName );
                auto const& rho = this->materialsProperties()->rho( matName );
                auto const& thermalExpansion = this->materialsProperties()->thermalExpansion( matName );
                auto rhoExpr = expr( rho.template expr<1,1>(), se );
                auto beta = expr( thermalExpansion.template expr<1,1>(), se );
                double T0 = M_BoussinesqRefTemperature;
                if ( doAssemblyFluid && M_useNaturalConvection )
                {
                    mylfV +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling_fluid*rhoExpr*(beta*(idv(t)-T0))*inner(M_gravityForce,id(u)),
                                   _geomap=this->geomap() );
                }

                if ( M_fluidModel->stabilizationGLS() && M_useNaturalConvection && !timeSteppingEvaluateResidualWithoutTimeDerivative )
                {
                    auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>( physicData->subphysicFromType( M_fluidModel->physicType() ) );
                    auto exprAdded = timeSteppingScaling_fluid*rhoExpr*(beta*(idv(t)-T0))*M_gravityForce;
                    auto exprAddedInGLSResidual = exprOptionalConcat< std::decay_t<decltype( exprAdded )> >();
                    exprAddedInGLSResidual.expression().add( exprAdded );
                    if ( !this->fluidModel()->isStationaryModel() && this->fluidModel()->timeStepping() ==  "Theta" )
                    {
                        auto const& mctxPrevious = mctx.additionalContext( "time-stepping.previous-model-context" );
                        auto const& sePrevious = mctxPrevious.symbolsExpr();
                        auto const& tOld = mctxPrevious.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );
                        auto rhoExprPrevious = expr( rho.template expr<1,1>(), sePrevious );
                        auto betaExprPrevious = expr( thermalExpansion.template expr<1,1>(), sePrevious );
                        auto exprAddedInRhsOld = (1.0-timeSteppingScaling_fluid)*rhoExprPrevious*(betaExprPrevious*(idv(tOld)-T0))*M_gravityForce;
                        if ( data.hasParameterValuesInfo( "time-stepping.previous-parameter-values" ) )
                            exprAddedInRhsOld.setParameterValues( data.parameterValuesInfo( "time-stepping.previous-parameter-values" ) );
                        exprAddedInGLSResidual.expression().add( exprAddedInRhsOld );
                    }
                    M_fluidModel->updateResidualStabilizationGLS( data, mctx, *physicFluidData, matProps, range, exprAddedInGLSResidual );
                }

            } //matName
        } // physic
    } // nonCstPart

    this->log("HeatFluid","updateResidual", "finish"+sc);
}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    M_heatModel->updateResidualDofElimination( data );
    M_fluidModel->updateResidualDofElimination( data );
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateResidual_Heat( DataUpdateResidual & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( vecCurrentSolution, this->heatModel()->startBlockSpaceIndexVector(),
                                    M_blockVectorSolution.vectorMonolithic(), this->startSubBlockSpaceIndex("fluid") );
    if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    {
        auto previousSolHeat = data.vectorInfo( "time-stepping.previous-solution");
        auto mctxPrevious = this->modelContextNoTrialSymbolsExpr(
            { { "solution", std::make_tuple( previousSolHeat, this->heatModel()->startBlockSpaceIndexVector() ) } },
            {
                { "solution", std::make_tuple( M_timeStepThetaSchemePreviousSolution, this->startSubBlockSpaceIndex("fluid") ) },
                { "velocity_extrapolated", std::make_tuple( M_fluidModel->vectorPreviousVelocityExtrapolated(), 0 ) }
            } );
        mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    }

    M_heatModel->updateResidual( data,mctx );
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateResidual_Fluid( DataUpdateResidual & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( this->heatModel()->blockVectorSolution().vectorMonolithic(), 0,
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

    this->log("HeatFluid","updateResidual_Fluid", "start" );


    M_fluidModel->updateResidual( data,mctx );

    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool BuildNonCstPart = !buildCstPart;

    if ( !M_useNaturalConvection )
        return;


    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false;
    if ( !M_fluidModel->isStationaryModel() )
    {
        timeSteppingScaling = data.doubleInfo( prefixvm(M_fluidModel->prefix(),"time-stepping.scaling") );
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm(M_fluidModel->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
    }

    auto XhV = M_fluidModel->functionSpaceVelocity();
    auto XhT = M_heatModel->spaceTemperature();
    auto const& u = mctx.field( fluid_model_type::FieldTag::velocity(this->fluidModel().get()), "velocity" );
    auto const& t = mctx.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );
    auto const& se = mctx.symbolsExpr();

    auto mylfV = form1( _test=XhV, _vector=R,
                        _rowstart=M_fluidModel->rowStartInVector() );

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(), matName );
            auto const& matProps = this->materialsProperties()->materialProperties( matName );
            auto rhoExpr = expr( matProps.property("density").template expr<1,1>(), se );
            auto betaExpr = expr( matProps.property("thermal-expansion").template expr<1,1>(), se );
            double T0 = M_BoussinesqRefTemperature;

            if ( buildCstPart )
            {
                mylfV +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling*rhoExpr*betaExpr*(idv(t)-T0)*inner(M_gravityForce,id(u)),
                               _geomap=this->geomap() );
            }

            if ( buildNonCstPart && this->fluidModel()->stabilizationGLS() && !timeSteppingEvaluateResidualWithoutTimeDerivative )
            {
                auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>( physicData->subphysicFromType( M_fluidModel->physicType() ) );
                auto exprAdded = timeSteppingScaling*rhoExpr*betaExpr*(idv(t)-T0)*M_gravityForce;
                auto exprAddedInGLSResidual = exprOptionalConcat< std::decay_t<decltype( exprAdded )> >();
                exprAddedInGLSResidual.expression().add( exprAdded );
                if ( !this->fluidModel()->isStationaryModel() && this->fluidModel()->timeStepping() ==  "Theta" )
                {
                    auto const& mctxPrevious = mctx.additionalContext( "time-stepping.previous-model-context" );
                    auto const& sePrevious = mctxPrevious.symbolsExpr();
                    auto const& tOld = mctxPrevious.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );
                    auto rhoExprPrevious = expr( matProps.property("density").template expr<1,1>(), sePrevious );
                    auto betaExprPrevious = expr( matProps.property("thermal-expansion").template expr<1,1>(), sePrevious );
                    auto exprAddedInRhsOld = (1.0-timeSteppingScaling)*rhoExprPrevious*betaExprPrevious*(idv(tOld)-T0)*M_gravityForce;
                    if ( data.hasParameterValuesInfo( "time-stepping.previous-parameter-values" ) )
                        exprAddedInRhsOld.setParameterValues( data.parameterValuesInfo( "time-stepping.previous-parameter-values" ) );
                    exprAddedInGLSResidual.expression().add( exprAddedInRhsOld );
                }
                M_fluidModel->updateResidualStabilizationGLS( data, mctx, *physicFluidData, matProps, range, exprAddedInGLSResidual );
            }
        }
    }
    this->log("HeatFluid","updateResidual_Fluid", "finish" );
}


} // namespace FeelModels
} // namespace Feel
