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
#if 0
    const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("HeatFluid","updateLinearPDE", "start"+sc);

    M_heatModel->updateLinearPDE( data );
    M_fluidModel->updateLinearPDE( data );

    if ( buildNonCstPart )
    {
        auto XhVP = M_fluidModel->spaceVelocityPressure();
        auto const& U = M_fluidModel->fieldVelocityPressure();
        auto u = U.template element<0>();
        auto XhT = M_heatModel->spaceTemperature();
        auto const& t = M_heatModel->fieldTemperature();

        auto mylfVP = form1( _test=XhVP, _vector=F,
                             _rowstart=M_fluidModel->rowStartInVector()+0 );

        auto mybfTT = form2( _test=XhT,_trial=XhT,_matrix=A,
                             _rowstart=M_heatModel->rowStartInMatrix(),
                             _colstart=M_heatModel->colStartInMatrix() );
        auto mybfVPT = form2( _test=XhVP,_trial=XhT,_matrix=A,
                             _rowstart=M_fluidModel->rowStartInMatrix(),
                             _colstart=M_heatModel->colStartInMatrix() );

        auto UConvection = M_fluidModel->timeStepBDF()->poly();
        auto uConvection = UConvection.template element<0>();

        for ( auto const& rangeData : this->rangeMeshElementsByMaterial() )
        {
            std::string const& matName = rangeData.first;
            auto const& range = rangeData.second;
            auto const& rhoHeatCapacity = M_heatModel->thermalProperties()->rhoHeatCapacity( matName );

            if ( rhoHeatCapacity.isConstant() )
            {
                double rhoHeatCapacityValue = rhoHeatCapacity.value();
                mybfTT +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityValue*(gradt(t)*idv(uConvection))*id(t),
                               _geomap=this->geomap() );
            }
            else
            {
                auto rhoHeatCapacityExpr = rhoHeatCapacity.expr();
                mybfTT +=
                    integrate( _range=range,
                               _expr= rhoHeatCapacityExpr*(gradt(t)*idv(uConvection))*id(t),
                               _geomap=this->geomap() );
            }
            auto const& rho = M_heatModel->thermalProperties()->rho( matName );
            auto const& thermalExpansion = M_heatModel->thermalProperties()->thermalExpansion( matName );
            CHECK( rhoHeatCapacity.isConstant() && thermalExpansion.isConstant() ) << "TODO";
            double rhoValue = rho.value();
            double beta = thermalExpansion.value();
            double T0 = M_BoussinesqRefTemperature;
            mybfVPT +=
                integrate( _range=range,
                           _expr= rhoValue*(beta*idt(t))*inner(M_gravityForce,id(u)),
                           _geomap=this->geomap() );
            mylfVP +=
                integrate( _range=range,
                           _expr= rhoValue*(beta*T0)*inner(M_gravityForce,id(u)),
                           _geomap=this->geomap() );

            if ( M_heatModel->stabilizationGLS() )
            {
                auto const& thermalConductivity = M_heatModel->thermalProperties()->thermalConductivity( matName );
                if ( thermalConductivity.isMatrix() )
                    CHECK( false ) << "TODO";
                else if ( thermalConductivity.isConstant() )
                    M_heatModel->updateLinearPDEStabilizationGLS( cst(rhoHeatCapacity.value()),cst(thermalConductivity.value()),idv(uConvection),range,data );
                else
                    CHECK( false ) << "TODO";
            }
            if ( M_fluidModel->stabilizationGLS() )
            {
                auto rhoF = idv(M_fluidModel->materialProperties()->fieldRho());
                //auto mu = Feel::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
                auto mu = idv(M_fluidModel->materialProperties()->fieldMu());
                auto expraddedInGLSResidualLF = rhoValue*beta*T0*M_gravityForce;
                auto exprAddedInGLSResidualBF = rhoValue*beta*idt(t)*M_gravityForce;
                M_fluidModel->updateLinearPDEStabilisationGLS( data, rhoF, mu, matName, hana::make_tuple(expraddedInGLSResidualLF),hana::make_tuple(std::make_pair(mybfVPT, exprAddedInGLSResidualBF)) );
            }

        }

    }

    this->log("HeatFluid","updateLinearPDE", "finish");
#else
    CHECK( false ) << "not allow here";
#endif
}




HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
#if 0
    M_heatModel->updateLinearPDEDofElimination( data );
    M_fluidModel->updateLinearPDEDofElimination( data );
#else
    CHECK( false ) << "not allow here";
#endif
}










HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateLinear_Heat( DataUpdateLinear & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( vecCurrentSolution, this->heatModel()->startBlockSpaceIndexVector(),
                                    this->fluidModel()->blockVectorSolution().vectorMonolithic(), 0 );
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
