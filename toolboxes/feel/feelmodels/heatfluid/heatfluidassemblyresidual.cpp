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
        auto mctxPrevious = this->modelContextNoTrialSymbolsExpr( previousSol, this->heatModel()->rowStartInVector(), this->fluidModel()->rowStartInVector() );
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
                if ( doAssemblyFluid )
                {
                    mylfV +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling_fluid*rhoExpr*(beta*(idv(t)-T0))*inner(M_gravityForce,id(u)),
                                   _geomap=this->geomap() );
                }

#if 0 // TODO VINCENT
                if ( M_heatModel->stabilizationGLS() && !timeSteppingEvaluateResidualWithoutTimeDerivative )
                {
                    auto const& thermalConductivity = this->materialsProperties()->thermalConductivity( matName );
                    CHECK ( !thermalConductivity.isMatrix() ) << "TODO";

                    if ( M_heatModel->timeStepping() == "Theta" )
                    {
                        auto previousSol = data.vectorInfo( prefixvm( this->prefix(),"time-stepping.previous-solution") );
                        auto tOld = XhT->element( previousSol, M_heatModel->rowStartInVector() );
                        auto uOld = XhV->element( previousSol, M_fluidModel->rowStartInVector() );
                        auto exprAddedInRhsOld = (1.0-timeSteppingScaling_fluid)*rhoHeatCapacityExpr*gradv(tOld)*idv(uOld);
                        M_heatModel->updateResidualStabilizationGLS( rhoHeatCapacityExpr,thermalConductivity.expr(),idv(u),range,data, exprAddedInRhsOld );
                    }
                    else
                        M_heatModel->updateResidualStabilizationGLS( rhoHeatCapacityExpr,thermalConductivity.expr(),idv(u),range,data );
                }
#endif

                if ( M_fluidModel->stabilizationGLS() && !timeSteppingEvaluateResidualWithoutTimeDerivative )
                {
                    auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>( physicData->subphysicFromType( M_fluidModel->physicType() ) );
                    auto expraddedInGLSResidual = timeSteppingScaling_fluid*rhoExpr*(beta*(idv(t)-T0))*M_gravityForce;
                    if ( !this->fluidModel()->isStationaryModel() && this->fluidModel()->timeStepping() ==  "Theta" )
                    {
                        auto const& mctxPrevious = mctx.additionalContext( "time-stepping.previous-model-context" );
                        auto const& tOld = mctxPrevious.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );
                        auto exprAddedInRhsOld = (1.0-timeSteppingScaling_fluid)*rhoExpr*(beta*(idv(tOld)-T0))*M_gravityForce;
                        M_fluidModel->updateResidualStabilizationGLS( data, mctx, *physicFluidData, matProps, range, expraddedInGLSResidual, exprAddedInRhsOld );
                    }
                    else
                        M_fluidModel->updateResidualStabilizationGLS( data, mctx, *physicFluidData, matProps, range, expraddedInGLSResidual );
#if 0
                    // auto XhP = M_fluidModel->functionSpacePressure();
                    // auto const p = XhP->element(XVec, M_fluidModel->rowStartInVector()+1 );
                    if ( M_fluidModel->timeStepping() ==  "Theta" )
                    {

                        auto previousSol = data.vectorInfo( prefixvm( this->prefix(),"time-stepping.previous-solution") );
                        auto tOld = XhT->element( previousSol, M_heatModel->rowStartInVector() );
                        auto exprAddedInRhsOld = (1.0-timeSteppingScaling_fluid)*rhoExpr*(beta*(idv(tOld)-T0))*M_gravityForce;
                        M_fluidModel->updateResidualStabilisationGLS( data, *u, *p, rhoF, mu, matName, expraddedInGLSResidual, exprAddedInRhsOld );
                    }
                    else
                        M_fluidModel->updateResidualStabilisationGLS( data, *u, *p, rhoF, mu, matName, expraddedInGLSResidual );
#endif
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
    const vector_ptrtype& XVec = data.currentSolution();
    // auto XhT = this->heatModel()->spaceTemperature();
    // size_type blockIndexTemperature = this->heatModel()->rowStartInVector() + this->heatModel()->startSubBlockSpaceIndex("temperature");
    // auto const t = XhT->element(XVec, blockIndexTemperature );
    // auto const& u = this->fluidModel()->fieldVelocity();
    // auto const& p = this->fluidModel()->fieldPressure();
    auto mfields = this->modelFields( this->heatModel()->modelFields( XVec, this->heatModel()->rowStartInVector(), this->heatModel()->keyword() ),
                                     this->fluidModel()->modelFields( this->fluidModel()->keyword() ) );
    auto mctx = Feel::FeelModels::modelContext( std::move(mfields), this->symbolsExpr( mfields ) );

    M_heatModel->updateResidual( data,mctx );
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateResidualFluidSolver( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool BuildNonCstPart = !buildCstPart;

    if ( buildNonCstPart )
        return;

    double timeSteppingScaling = 1.;
    if ( !M_fluidModel->isStationaryModel() )
        timeSteppingScaling = data.doubleInfo( prefixvm(M_fluidModel->prefix(),"time-stepping.scaling") );

    this->log("HeatFluid","updateResidualFluidSolver", "start" );

    auto XhV = M_fluidModel->functionSpaceVelocity();
    auto const u = XhV->element(XVec, M_fluidModel->rowStartInVector()+0 );

    auto XhT = M_heatModel->spaceTemperature();
    auto const& t = M_heatModel->fieldTemperature();

    auto mylfV = form1( _test=XhV, _vector=R,
                        _rowstart=M_fluidModel->rowStartInVector() );

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(), matName );
            auto const& rhoHeatCapacity = this->materialsProperties()->rhoHeatCapacity( matName );

            auto const& rho = this->materialsProperties()->rho( matName );
            auto const& thermalExpansion = this->materialsProperties()->thermalExpansion( matName );
            auto const& rhoExpr = rho.expr();
            auto const& betaExpr = thermalExpansion.expr();
            double T0 = M_BoussinesqRefTemperature;

            mylfV +=
                integrate( _range=range,
                           _expr= timeSteppingScaling*rhoExpr*betaExpr*(idv(t)-T0)*inner(M_gravityForce,id(u)),
                           _geomap=this->geomap() );
        }
    }
    this->log("HeatFluid","updateResidualFluidSolver", "finish" );
}


} // namespace FeelModels
} // namespace Feel
