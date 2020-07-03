/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/thermoelectric/thermoelectric.hpp>

namespace Feel
{
namespace FeelModels
{

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("ThermoElectric","updateLinearPDE", "start"+sc);

    size_type startBlockIndexTemperature = M_heatModel->startBlockSpaceIndexVector()+0;
    size_type startBlockIndexElectricPotential = M_electricModel->startBlockSpaceIndexVector()+0;
    auto mesh = this->mesh();
    auto XhV = M_electricModel->spaceElectricPotential();
    auto XhT = M_heatModel->spaceTemperature();

    auto mctx = this->modelContext( vecCurrentPicardSolution, M_heatModel->startBlockSpaceIndexVector(), M_electricModel->startBlockSpaceIndexVector() );
    auto const& symbolsExpr = mctx.symbolsExpr();
    auto const& v = mctx.field( electric_model_type::FieldTag::potential(this->electricModel().get()), "electric-potential" );
    auto const& t = mctx.field( heat_model_type::FieldTag::temperature(this->heatModel().get()), "temperature" );

    //auto symbolsExpr = this->symbolsExpr( t, v );
    M_heatModel->updateLinearPDE( data,mctx );
    M_electricModel->updateLinearPDE( data,mctx );

    if ( buildNonCstPart )
    {
        auto mybfTV = form2( _test=XhT,_trial=XhV,_matrix=A,
                             _pattern=size_type(Pattern::COUPLED),
                             _rowstart=this->rowStartInMatrix()+startBlockIndexTemperature ,
                             _colstart=this->colStartInMatrix()+startBlockIndexElectricPotential );

        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(), matName );
            auto const& electricConductivity = this->materialsProperties()->electricConductivity( matName );
            if ( M_modelUseJouleEffect )
            {
                auto sigmaExpr = expr( electricConductivity.expr(), symbolsExpr );
                mybfTV +=
                    integrate( _range=range,
                               _expr= -sigmaExpr*inner(gradv(v),gradt(v))*id( t ),
                               _geomap=this->geomap() );
            }
        }
    }

    this->log("ThermoElectric","updateLinearPDE","finish" );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    const vector_ptrtype& vecCurrentSolution = data.currentSolution();
    auto mctx = this->modelContext( vecCurrentSolution, M_heatModel->startBlockSpaceIndexVector(), M_electricModel->startBlockSpaceIndexVector() );
    M_heatModel->updateLinearPDEDofElimination( data, mctx );
    M_electricModel->updateLinearPDEDofElimination( data, mctx );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateLinear_Electric( DataUpdateLinear & data ) const
{
    this->log("ThermoElectric","updateLinear_Electric","start" );

    auto mctx = this->modelContext( /*vecCurrentPicardSolution, startBlockIndexTemperature, startBlockIndexElectricPotential*/ );
    M_electricModel->updateLinearPDE( data,mctx );

    this->log("ThermoElectric","updateLinear_Electric","finish" );
}

THERMOELECTRIC_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRIC_CLASS_TEMPLATE_TYPE::updateLinear_Heat( DataUpdateLinear & data ) const
{
    this->log("ThermoElectric","updateLinear_Heat","start" );

    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    vector_ptrtype& F = data.rhs();

    auto mctx = this->modelContext( /*vecCurrentPicardSolution, startBlockIndexTemperature, startBlockIndexElectricPotential*/ );

    auto const& v = M_electricModel->fieldElectricPotential();
    auto const& t = M_heatModel->fieldTemperature();
    auto const& symbolsExpr = mctx.symbolsExpr();

    M_heatModel->updateLinearPDE( data,mctx );
    if ( buildNonCstPart && M_modelUseJouleEffect ) // TODO : not always non cst part
    {
        auto XhT = M_heatModel->spaceTemperature();
        auto myLinearForm = form1( _test=XhT,_vector=F,
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
                               _expr= sigmaExpr*inner(gradv(v))*id(t),
                               _geomap=this->geomap() );
            }
        }
    }

    this->log("ThermoElectric","updateLinear_Heat","finish" );
}


} // namespace FeelModels
} // namespace Feel
