/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/thermoelectric/thermoelectric.hpp>

namespace Feel
{
namespace FeelModels
{

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool useJacobianLinearTerms = data.useJacobianLinearTerms();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("ThermoElectric","updateResidual", "start"+sc);

    size_type startBlockIndexTemperature = M_heatModel->startBlockSpaceIndexVector()+0;
    size_type startBlockIndexElectricPotential = M_electricModel->startBlockSpaceIndexVector()+0;

    auto mctx = this->modelContext( XVec, this->rowStartInVector()+startBlockIndexTemperature, this->rowStartInVector()+startBlockIndexElectricPotential );
    auto const& symbolsExpr = mctx.symbolsExpr();
    auto const& v = mctx.field( electric_model_type::FieldTag::potential(this->electricModel().get()), "electric-potential" );
    auto const& t = mctx.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );


    auto mesh = this->mesh();

    auto XhV = M_electricModel->spaceElectricPotential();
    //auto const v = XhV->element(XVec, this->rowStartInVector()+startBlockIndexElectricPotential );
    auto XhT = M_heatModel->spaceTemperature();
    //auto const t = XhT->element(XVec, this->rowStartInVector()+startBlockIndexTemperature );
    //auto symbolsExpr = this->symbolsExpr(t,v);
    // auto mfield = this->modelFields( XVec, this->rowStartInVector()+startBlockIndexTemperature, this->rowStartInVector()+startBlockIndexElectricPotential );
    // auto symbolsExpr = this->symbolsExpr( mfield );

    M_heatModel->updateResidual( data,mctx );
    M_electricModel->updateResidual( data,mctx );

    if ( !buildCstPart )
    {
        auto mylfT = form1( _test=XhT, _vector=R,
                            _rowstart=this->rowStartInVector()+startBlockIndexTemperature );
        auto mylfV = form1( _test=XhV, _vector=R,
                            _rowstart=this->rowStartInVector()+startBlockIndexElectricPotential );

        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        {
            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
            {
                auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                auto const& electricConductivity = this->materialsProperties()->electricConductivity( matName );
                auto sigmaExpr = expr( electricConductivity.expr(), symbolsExpr );
                if ( M_modelUseJouleEffect )
                {
                    mylfT +=
                        integrate( _range=range,
                                   _expr= -sigmaExpr*inner(gradv(v))*id( t ),
                                   _geomap=this->geomap() );
                }
            }
        }
    }

    this->log("ThermoElectric","updateResidual", "finish"+sc);
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    M_heatModel->updateResidualDofElimination( data );
    M_electricModel->updateResidualDofElimination( data );
}


THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateResidual_Heat( DataUpdateResidual & data ) const
{
    this->log("ThermoElectric","updateResidual_Heat","start" );

    const vector_ptrtype& XVec = data.currentSolution();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    vector_ptrtype& R = data.residual();

    auto mfields = this->modelFields( this->heatModel()->modelFields( XVec, this->heatModel()->rowStartInVector(), this->heatModel()->keyword() ),
                                     this->electricModel()->modelFields( this->electricModel()->keyword() ) );
    //auto symbolsExpr = this->symbolsExpr( mfields );
    auto mctx = Feel::FeelModels::modelContext( std::move(mfields), this->symbolsExpr( mfields ) );
    auto const& symbolsExpr = mctx.symbolsExpr();
    auto const& v = mctx.field( electric_model_type::FieldTag::potential(this->electricModel().get()), "electric-potential" );
    auto const& t = mctx.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );

    auto XhT = M_heatModel->spaceTemperature();
    //auto const& v = M_electricModel->fieldElectricPotential();
    //auto const t = XhT->element(XVec, M_heatModel->rowStartInVector());
    //auto symbolsExpr = this->symbolsExpr(t,v);

    M_heatModel->updateResidual( data,mctx );

    if ( buildNonCstPart && M_modelUseJouleEffect ) // TODO non const part only if sigma is related to Temperature
    {
        auto myLinearForm = form1( _test=XhT,_vector=R,
                                   _rowstart=M_heatModel->rowStartInVector() );

        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        {
            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
            {
                auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                auto const& electricConductivity = this->materialsProperties()->electricConductivity( matName );
                auto sigmaExpr = expr( electricConductivity.expr(), symbolsExpr );
                myLinearForm +=
                    integrate( _range=range,
                               _expr= -sigmaExpr*inner(gradv(v))*id(t),
                               _geomap=this->geomap() );
            }
        }
    }

    this->log("ThermoElectric","updateResidual_Heat","finish" );
}


} // namespace FeelModels
} // namespace Feel
