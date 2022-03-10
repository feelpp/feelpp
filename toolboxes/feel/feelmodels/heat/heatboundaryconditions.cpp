/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/heat/heatboundaryconditions.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel
{
namespace FeelModels
{

void
HeatBoundaryConditions::HeatFlux::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
{
    if ( jarg.contains( "expr" ) )
        M_mexpr.setExpr( jarg.at( "expr" ), mparent.worldComm(), mparent.repository().expr(), indexes );
     if ( jarg.contains( "markers" ) )
     {
         ModelMarkers markers;
         markers.setup( jarg.at("markers"), indexes );
         M_markers = markers;
     }
     else
         M_markers = { M_name };
}

void
HeatBoundaryConditions::HeatFlux::updateInformationObject( nl::json & p ) const
{
    auto [exprStr,compInfo] = M_mexpr.exprInformations();
    p["expr"] = exprStr;
    p["markers"] = M_markers;
}

tabulate_informations_ptr_t
HeatBoundaryConditions::HeatFlux::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
    TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { /*"name","type",*/"expr" } );
    tabInfo.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    if ( jsonInfo.contains("markers") )
    {
        Feel::Table tabInfoMarkers = TabulateInformationTools::FromJSON::createTableFromArray( jsonInfo.at("markers") , true );
        if ( tabInfoMarkers.nRow() > 0 )
            tabInfo.add_row( { "markers", tabInfoMarkers } );
    }
    return TabulateInformations::New( tabInfo, tabInfoProp );
}
void
HeatBoundaryConditions::ConvectiveHeatFlux::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
{
    if ( jarg.contains( "h" ) )
        M_mexpr_h.setExpr( jarg.at( "h" ), mparent.worldComm(), mparent.repository().expr(), indexes );
    if ( jarg.contains( "Text" ) )
        M_mexpr_Text.setExpr( jarg.at( "Text" ), mparent.worldComm(), mparent.repository().expr(), indexes );
     if ( jarg.contains( "markers" ) )
     {
         ModelMarkers markers;
         markers.setup( jarg.at("markers"), indexes );
         M_markers = markers;
     }
     else
         M_markers = { M_name };
}

void
HeatBoundaryConditions::ConvectiveHeatFlux::updateInformationObject( nl::json & p ) const
{
    auto [exprStr_h,compInfo_h] = M_mexpr_h.exprInformations();
    p["expr_h"] = exprStr_h;
    auto [exprStr_Text,compInfo_Text] = M_mexpr_Text.exprInformations();
    p["expr_Text"] = exprStr_Text;
    p["markers"] = M_markers;
}

tabulate_informations_ptr_t
HeatBoundaryConditions::ConvectiveHeatFlux::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
    TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { /*"name","type",*/"expr_h","expr_Text" } );
    tabInfo.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    if ( jsonInfo.contains("markers") )
    {
        Feel::Table tabInfoMarkers = TabulateInformationTools::FromJSON::createTableFromArray( jsonInfo.at("markers") , true );
        if ( tabInfoMarkers.nRow() > 0 )
            tabInfo.add_row( { "markers", tabInfoMarkers } );
    }
    return TabulateInformations::New( tabInfo, tabInfoProp );
}

void
HeatBoundaryConditions::setup( nl::json const& jarg )
{
    auto tbParent = this->toolboxParent();
    ModelIndexes indexes;
    for ( std::string const& bcKeyword : { "temperature_imposed", "temperature" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_temp = jarg.at( bcKeyword );
            for ( auto const& [j_tempkey,j_tempval] : j_temp.items() )
            {
                auto bc = std::make_shared<TemperatureImposed>( j_tempkey, tbParent );
                bc->setup( j_tempval,indexes );
                M_temperatureImposed.emplace(j_tempkey, std::move( bc ) );
            }
        }
    }
    for ( std::string const& bcKeyword : { "heat_flux","flux" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_temp = jarg.at( bcKeyword );
            for ( auto const& [j_tempkey,j_tempval] : j_temp.items() )
            {
                auto bc = std::make_shared<HeatFlux>( j_tempkey );
                bc->setup( *tbParent,j_tempval,indexes );
                M_heatFlux.emplace(j_tempkey, std::move( bc ) );
            }
        }
    }
    for ( std::string const& bcKeyword : { "convective_heat_flux" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_temp = jarg.at( bcKeyword );
            for ( auto const& [j_tempkey,j_tempval] : j_temp.items() )
            {
                auto bc = std::make_shared<ConvectiveHeatFlux>( j_tempkey );
                bc->setup( *tbParent,j_tempval,indexes );
                M_convectiveHeatFlux.emplace(j_tempkey, std::move( bc ) );
            }
        }
    }
}

void
HeatBoundaryConditions::setParameterValues( std::map<std::string,double> const& paramValues )
{
    for ( auto & [bcname,bcData] : M_temperatureImposed )
        bcData->setParameterValues( paramValues );
    for ( auto & [bcname,bcData] : M_heatFlux )
        bcData->setParameterValues( paramValues );
    for ( auto & [bcname,bcData] : M_convectiveHeatFlux )
        bcData->setParameterValues( paramValues );
}

void
HeatBoundaryConditions::updateInformationObject( nl::json & p ) const
{
    if ( !M_temperatureImposed.empty() )
    {
        nl::json & pBC = p["temperature_imposed"];
        for ( auto const& [bcname,bcData] : M_temperatureImposed )
            bcData->updateInformationObject( pBC[bcname] );
    }
    if ( !M_heatFlux.empty() )
    {
        nl::json & pBC = p["heat_flux"];
        for ( auto const& [bcname,bcData] : M_heatFlux )
            bcData->updateInformationObject( pBC[bcname] );
    }
    if ( !M_convectiveHeatFlux.empty() )
    {
        nl::json & pBC = p["convective_heat_flux"];
        for ( auto const& [bcname,bcData] : M_convectiveHeatFlux )
            bcData->updateInformationObject( pBC[bcname] );
    }

}

tabulate_informations_ptr_t
HeatBoundaryConditions::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains( "temperature_imposed" ) )
    {
        auto tabInfoTemperatureImposed = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_tempkey,j_tempval]: jsonInfo.at( "temperature_imposed" ).items() )
            tabInfoTemperatureImposed->add( j_tempkey, HeatBoundaryConditions::TemperatureImposed::tabulateInformations( j_tempval, tabInfoProp ) );
        tabInfo->add( "Temperature Imposed", tabInfoTemperatureImposed );
    }
    if ( jsonInfo.contains( "heat_flux" ) )
    {
        auto tabInfoHeatFlux = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_tempkey,j_tempval]: jsonInfo.at( "heat_flux" ).items() )
            tabInfoHeatFlux->add( j_tempkey, HeatBoundaryConditions::HeatFlux::tabulateInformations( j_tempval, tabInfoProp ) );
        tabInfo->add( "Heat Flux", tabInfoHeatFlux );
    }
    if ( jsonInfo.contains( "convective_heat_flux" ) )
    {
        auto tabInfoConvectiveHeatFlux = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_tempkey,j_tempval]: jsonInfo.at( "convective_heat_flux" ).items() )
            tabInfoConvectiveHeatFlux->add( j_tempkey, HeatBoundaryConditions::ConvectiveHeatFlux::tabulateInformations( j_tempval, tabInfoProp ) );
        tabInfo->add( "Convective Heat Flux", tabInfoConvectiveHeatFlux );
    }
    return tabInfo;
}

} // namespace FeelModels
} // namespace Feel
