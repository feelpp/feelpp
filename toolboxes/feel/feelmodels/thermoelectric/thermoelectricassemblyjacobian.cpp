/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/thermoelectric/thermoelectric.hpp>

namespace Feel
{
namespace FeelModels
{

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    this->log("ThermoElectric","updateNewtonInitialGuess","start" );
    vector_ptrtype& U = data.initialGuess();
    auto mctx = this->modelContext( U, M_heatModel->startBlockSpaceIndexVector(), M_electricModel->startBlockSpaceIndexVector() );
    M_heatModel->updateNewtonInitialGuess( data, mctx );
    M_electricModel->updateNewtonInitialGuess( data, mctx );
    this->log("ThermoElectric","updateNewtonInitialGuess","finish" );
}
THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("ThermoElectric","updateJacobian", "start"+sc);
    size_type startBlockIndexTemperature = M_heatModel->startBlockSpaceIndexVector()+0;
    size_type startBlockIndexElectricPotential = M_electricModel->startBlockSpaceIndexVector()+0;

    auto mctx = this->modelContext( XVec, this->rowStartInVector()+startBlockIndexTemperature, this->rowStartInVector()+startBlockIndexElectricPotential ); // TODO CHECK if we need rowStartInVector here
    auto const& symbolsExpr = mctx.symbolsExpr();
    auto const& v = mctx.field( electric_model_type::FieldTag::potential(this->electricModel().get()), "electric-potential" );
    auto const& t = mctx.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );

    auto mesh = this->mesh();

    auto XhV = M_electricModel->spaceElectricPotential();
    //auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );
    auto XhT = M_heatModel->spaceTemperature();
    //auto const& t = M_heatModel->fieldTemperature();
    //auto const t = XhT->element(XVec, this->rowStartInVector()+startBlockIndexTemperature );
    //auto symbolsExpr = this->symbolsExpr(t,v);

    // auto mfield = this->modelFields( XVec, this->rowStartInVector()+startBlockIndexTemperature, this->rowStartInVector()+startBlockIndexElectricPotential );
    // auto symbolsExpr = this->symbolsExpr( mfield );

    M_heatModel->updateJacobian( data,mctx );
    M_electricModel->updateJacobian( data,mctx );

    if ( !buildCstPart )
    {
        auto mybfTT = form2( _test=XhT,_trial=XhT,_matrix=J,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexTemperature,
                             _colstart=this->colStartInMatrix()+startBlockIndexTemperature );
        auto mybfTV = form2( _test=XhT,_trial=XhV,_matrix=J,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexTemperature,
                             _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );
        auto mybfVT = form2( _test=XhV,_trial=XhT,_matrix=J,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexElectricPotential,
                             _colstart=this->colStartInMatrix()+startBlockIndexTemperature );

        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( matName );
            auto const& electricConductivity = this->materialsProperties()->electricConductivity( matName );
            std::string symbolStr = "heat_T";
            auto sigmaExpr = expr( electricConductivity.expr(), symbolsExpr );
            if ( M_modelUseJouleEffect )
            {
                mybfTV +=
                    integrate( _range=range,
                               _expr= -sigmaExpr*2*inner(gradt(v),gradv(v))*id( t ),
                               _geomap=this->geomap() );
            }

            if ( sigmaExpr.expression().hasSymbol( symbolStr ) )
            {
                auto sigmaDiffExpr = diff( sigmaExpr,symbolStr,1,"",this->worldComm(),this->repository().expr());
                if ( M_modelUseJouleEffect )
                {
                    mybfTT +=
                        integrate( _range=range,
                                   _expr= -sigmaDiffExpr*idt(t)*inner(gradv(v)/*,gradv(v)*/)*id( t ),
                                   _geomap=this->geomap() );
                }

                mybfVT +=
                    integrate( _range=range,
                               _expr= sigmaDiffExpr*idt(t)*inner(gradv(v),grad(v)),
                               _geomap=this->geomap() );
            }
        }
    }

    this->log("ThermoElectric","updateJacobian", "finish"+sc);
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    M_heatModel->updateJacobianDofElimination( data );
    M_electricModel->updateJacobianDofElimination( data );
}

} // namespace FeelModels
} // namespace Feel
