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
    M_fluidModel->updateNewtonInitialGuess( data );
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
    auto const& symbolsExpr = mctx.symbolsExpr();
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
    M_fluidModel->updateJacobian( data/*, symbolsExpr*/ );

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
                auto const& rho = this->materialsProperties()->rho( matName );
                auto const& rhoHeatCapacity = this->materialsProperties()->rhoHeatCapacity( matName );
                auto const& thermalExpansion = this->materialsProperties()->thermalExpansion( matName );

                auto rhoHeatCapacityExpr = rhoHeatCapacity.expr();
                auto rhoExpr = rho.expr();
                auto beta = thermalExpansion.expr();
#if 0
                form2( _test=XhT,_trial=XhT,_matrix=J,
                       _rowstart=M_heatModel->rowStartInMatrix(),
                       _colstart=M_heatModel->colStartInMatrix() ) +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling_heat*rhoHeatCapacityExpr*(gradt(t)*idv(u))*id(t),
                               _geomap=this->geomap() );
#endif
                form2( _test=XhT,_trial=XhV,_matrix=J,
                       _rowstart=M_heatModel->rowStartInMatrix(),
                       _colstart=M_fluidModel->colStartInMatrix() ) +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling_heat*rhoHeatCapacityExpr*(gradv(t)*idt(u))*id(t),
                               _geomap=this->geomap() );

                bfVT +=
                    integrate( _range=range,
                               _expr= timeSteppingScaling_fluid*rhoExpr*beta*idt(t)*inner(M_gravityForce,id(u)),
                               _geomap=this->geomap() );

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

                if ( M_fluidModel->stabilizationGLS() )
                {
                    auto rhoF = idv(M_fluidModel->materialProperties()->fieldRho());
                    //auto mu = Feel::FeelModels::fluidMecViscosity<2*FluidMechanicsType::nOrderVelocity>(u,p,*fluidmec.materialProperties());
                    auto mu = idv(M_fluidModel->materialProperties()->fieldMu());
                    auto exprAddedInGLSResidual = rhoExpr*beta*idt(t)*M_gravityForce;

                    auto XhP = M_fluidModel->functionSpacePressure();
                    //auto const p = XhP->element(XVec, M_fluidModel->rowStartInVector()+1 );
                    M_fluidModel->updateJacobianStabilisationGLS( data, *u, *p, rhoF, mu, matName, std::make_pair(bfVT, exprAddedInGLSResidual) );
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
    const vector_ptrtype& XVec = data.currentSolution();
    auto mfields = this->modelFields( this->heatModel()->modelFields( XVec, this->heatModel()->rowStartInVector(), this->heatModel()->keyword() ),
                                      this->fluidModel()->modelFields( this->fluidModel()->keyword() ) );
    auto mctx = Feel::FeelModels::modelContext( std::move(mfields), this->symbolsExpr( mfields ) );
    M_heatModel->updateJacobian( data,mctx );
}


} // namespace FeelModels
} // namespace Feel
