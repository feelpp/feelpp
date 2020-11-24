/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/heatfluid/heatfluid.hpp>

namespace Feel
{
namespace FeelModels
{

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    this->log("HeatFluid","updateNewtonInitialGuess","start" );
    vector_ptrtype& U = data.initialGuess();
    auto mctx = this->modelContext( U, this->heatModel(), this->fluidModel() );
    M_heatModel->updateNewtonInitialGuess( data, mctx );
    M_fluidModel->updateNewtonInitialGuess( data, mctx );
    this->log("HeatFluid","updateNewtonInitialGuess","finish" );
}
HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("HeatFluid","updateJacobian", "start"+sc);

    auto mctx = this->modelContext( XVec, this->heatModel(), this->fluidModel() );
    auto const& se = mctx.symbolsExpr();
    auto const& t = mctx.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );
    auto const& u = mctx.field( fluid_model_type::FieldTag::velocity(this->fluidModel().get()), "velocity" );
    auto const& p = mctx.field( fluid_model_type::FieldTag::pressure(this->fluidModel().get()), "pressure" );

    auto mesh = this->mesh();

    auto XhT = this->heatModel()->spaceTemperature();
    auto XhV = this->fluidModel()->functionSpaceVelocity();
    auto XhP = this->fluidModel()->functionSpacePressure();
    size_type blockIndexTemperature = this->heatModel()->startBlockSpaceIndexVector() + this->heatModel()->startSubBlockSpaceIndex("temperature");
    size_type blockIndexVelocity = this->fluidModel()->startBlockSpaceIndexVector() + this->fluidModel()->startSubBlockSpaceIndex("velocity");
    size_type blockIndexPressure = this->fluidModel()->startBlockSpaceIndexVector() + this->fluidModel()->startSubBlockSpaceIndex("pressure");

    M_heatModel->updateJacobian( data, mctx );
    M_fluidModel->updateJacobian( data, mctx );

    if ( buildNonCstPart )
    {
        double timeSteppingScaling_fluid = 1.;
        if ( !M_fluidModel->isStationaryModel() )
            timeSteppingScaling_fluid = data.doubleInfo( prefixvm(M_fluidModel->prefix(),"time-stepping.scaling") );
        double timeSteppingScaling_heat = 1.;
        if ( !M_heatModel->isStationary() )
            timeSteppingScaling_heat = data.doubleInfo( prefixvm(M_heatModel->prefix(),"time-stepping.scaling") );

        auto bfVT = form2( _test=XhV,_trial=XhT,_matrix=J,
                            _rowstart=M_fluidModel->rowStartInMatrix(),
                            _colstart=M_heatModel->colStartInMatrix() );

        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        {
            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
            {
                auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                auto const& matProps = this->materialsProperties()->materialProperties( matName );
                auto const& rho = this->materialsProperties()->rho( matName );
                auto const& rhoHeatCapacity = this->materialsProperties()->rhoHeatCapacity( matName );
                auto const& thermalExpansion = this->materialsProperties()->thermalExpansion( matName );

                auto rhoHeatCapacityExpr = expr( rhoHeatCapacity.expr(), se );
                auto rhoExpr = expr( rho.template expr<1,1>(), se );
                auto beta = expr( thermalExpansion.template expr<1,1>(), se );

                if ( !M_fluidModel->useSemiImplicitTimeScheme() )
                {
                    form2( _test=XhT,_trial=XhV,_matrix=J,
                           _rowstart=M_heatModel->rowStartInMatrix(),
                           _colstart=M_fluidModel->colStartInMatrix() ) +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling_heat*rhoHeatCapacityExpr*(gradv(t)*idt(u))*id(t),
                                   _geomap=this->geomap() );
                }

                if ( M_useNaturalConvection )
                {
                    bfVT +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling_fluid*rhoExpr*beta*idt(t)*inner(M_gravityForce,id(u)),
                                   _geomap=this->geomap() );
                }

                if ( M_heatModel->stabilizationGLS() )
                {
#if 0
                    auto const& thermalConductivity = M_heatModel->thermalProperties()->thermalConductivity( matName );
                    if ( thermalConductivity.isMatrix() )
                        CHECK( false ) << "TODO";
                    else
                        M_heatModel->updateJacobianStabilizationGLS( rhoHeatCapacityExpr,thermalConductivity.expr(),idv(u),range,data );
#endif
#if 0
                    auto rhocp = rhoHeatCapacityExpr;
                    auto tau = idv( M_heatModel->stabilizationGLSParameter()->fieldTauPtr() );
                    if ( true ) // order==1 or supg
                    {
                        auto stab_test = rhocp*grad(t)*idv(u);
                        auto stab_residual_bilinear = rhocp*gradv(t)*idt(u);
                        form2( _test=XhT,_trial=XhV,_matrix=J,
                               _rowstart=M_heatModel->rowStartInMatrix(),
                               _colstart=M_fluidModel->colStartInMatrix() ) +=
                            integrate( _range=range,
                                       _expr=tau*stab_residual_bilinear*stab_test,
                                       _geomap=this->geomap() );
                    
                        auto stab_test2 = rhocp*gradv(t)*idt(u);
                        form2( _test=XhVP,_trial=XhV,_matrix=J,
                               _rowstart=M_fluidModel->rowStartInMatrix(),
                               _colstart=M_fluidModel->colStartInMatrix() ) +=
                            integrate( _range=range,
                                       _expr=tau*stab_residual_bilinear*stab_test2,
                                       _geomap=this->geomap() );
                        auto stab_residual_bilinear2 = rhocp*(idt(t)*M_heatModel->timeStepBdfTemperature()->polyDerivCoefficient(0) + gradt(t)*idv(u) );
                        form2( _test=XhV,_trial=XhT,_matrix=J,
                               _rowstart=M_fluidModel->rowStartInMatrix(),
                               _colstart=M_heatModel->colStartInMatrix() ) +=
                            integrate( _range=range,
                                       _expr=tau*stab_residual_bilinear2*stab_test2,
                                       _geomap=this->geomap() );
                    
                    }
                    else
                    {

                    }
#endif
                }

                if ( M_fluidModel->stabilizationGLS() && M_useNaturalConvection )
                {
                    auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>( physicData->subphysicFromType( M_fluidModel->physicType() ) );
                    auto exprAddedInGLSResidual = rhoExpr*beta*idt(t)*M_gravityForce;
                    M_fluidModel->updateJacobianStabilizationGLS( data, mctx, *physicFluidData, matProps, range, std::make_tuple(XhT,blockIndexTemperature,exprAddedInGLSResidual) );
                }
            } // matName
        } // physic
    } // nonCstPart

    this->log("HeatFluid","updateJacobian", "finish"+sc);

}

HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    M_heatModel->updateJacobianDofElimination( data );
    M_fluidModel->updateJacobianDofElimination( data );
}




HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateJacobian_Heat( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( vecCurrentSolution, this->heatModel()->startBlockSpaceIndexVector(),
                                    this->fluidModel()->blockVectorSolution().vectorMonolithic(), 0 );

    M_heatModel->updateJacobian( data,mctx );
}


HEATFLUID_CLASS_TEMPLATE_DECLARATIONS
void
HEATFLUID_CLASS_TEMPLATE_TYPE::updateJacobian_Fluid( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( this->heatModel()->blockVectorSolution().vectorMonolithic(), 0,
                                    vecCurrentSolution, this->fluidModel()->startBlockSpaceIndexVector() );

    bool currentSabDoGLSDoAssembly = M_fluidModel->stabilizationGLSDoAssembly();
    M_fluidModel->setStabilizationGLSDoAssembly( true );
    M_fluidModel->updateJacobian( data, mctx );
    M_fluidModel->setStabilizationGLSDoAssembly( currentSabDoGLSDoAssembly );
    // for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    // {
    //     for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
    //     {
    //         auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
    //         auto const& matProps = this->materialsProperties()->materialProperties( matName );
    //     }
    // }
}

} // namespace FeelModels
} // namespace Feel
