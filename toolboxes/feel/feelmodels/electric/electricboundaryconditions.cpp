/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/electric/electricboundaryconditions.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel
{
namespace FeelModels
{

void
ElectricBoundaryConditions::ElectricPotentialImposed::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
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
ElectricBoundaryConditions::ElectricPotentialImposed::updateInformationObject( nl::json & p ) const
{
    auto [exprStr,compInfo] = M_mexpr.exprInformations();
    p["expr"] = exprStr;
    p["markers"] = M_markers;
}

tabulate_informations_ptr_t
ElectricBoundaryConditions::ElectricPotentialImposed::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
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
ElectricBoundaryConditions::Ground::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
{
    nl::json j_bc = jarg;
    j_bc["expr"] = "0";
    super_type::setup( mparent,j_bc,indexes );
}
void
ElectricBoundaryConditions::Ground::updateInformationObject( nl::json & p ) const
{
    // auto [exprStr,compInfo] = M_mexpr.exprInformations();
    // p["expr"] = exprStr;
    p["markers"] = this->M_markers;
}

tabulate_informations_ptr_t
ElectricBoundaryConditions::Ground::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
    if ( jsonInfo.contains("markers") )
    {
        Feel::Table tabInfoMarkers = TabulateInformationTools::FromJSON::createTableFromArray( jsonInfo.at("markers") , true );
        if ( tabInfoMarkers.nRow() > 0 )
            tabInfo.add_row( { "markers", tabInfoMarkers } );
    }
    return TabulateInformations::New( tabInfo, tabInfoProp );
}

void
ElectricBoundaryConditions::SurfaceChargeDensity::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
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
ElectricBoundaryConditions::SurfaceChargeDensity::updateInformationObject( nl::json & p ) const
{
    auto [exprStr,compInfo] = M_mexpr.exprInformations();
    p["expr"] = exprStr;
    p["markers"] = M_markers;
}

tabulate_informations_ptr_t
ElectricBoundaryConditions::SurfaceChargeDensity::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
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
ElectricBoundaryConditions::setup( ModelBase const& mparent, nl::json const& jarg )
{
    ModelIndexes indexes;
    for ( std::string const& bcKeyword : { "electric_potential_imposed", "electric_potential" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<ElectricPotentialImposed>( j_bckey );
                bc->setup( mparent,j_bcval,indexes );
                M_electricPotentialImposed.emplace( std::make_pair(Type::ElectricPotentialImposed,j_bckey), std::move( bc ) );
            }
        }
    }
    for ( std::string const& bcKeyword : { "ground" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            std::string bcName = "";
            auto bc = std::make_shared<Ground>( bcName );
            bc->setup( mparent,j_bc,indexes );
            M_electricPotentialImposed.emplace( std::make_pair(Type::Ground,bcName), std::move( bc ) );
        }
    }
    for ( std::string const& bcKeyword : { "surface_charge_density","charge_density" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<SurfaceChargeDensity>( j_bckey );
                bc->setup( mparent,j_bcval,indexes );
                M_surfaceChargeDensity.emplace(j_bckey, std::move( bc ) );
            }
        }
    }
}

void
ElectricBoundaryConditions::setParameterValues( std::map<std::string,double> const& paramValues )
{
    for ( auto & [bcId,bcData] : M_electricPotentialImposed )
        bcData->setParameterValues( paramValues );
    for ( auto & [bcname,bcData] : M_surfaceChargeDensity )
        bcData->setParameterValues( paramValues );
}

void
ElectricBoundaryConditions::updateInformationObject( nl::json & p ) const
{
    if ( !M_electricPotentialImposed.empty() )
    {
        //nl::json & pBC = p["electric_potential_imposed"];
        for ( auto const& [bcId,bcData] : M_electricPotentialImposed )
        {
            if ( bcId.first == Type::ElectricPotentialImposed )
                bcData->updateInformationObject( p["electric_potential_imposed"][bcId.second] );
            else if ( bcId.first == Type::Ground )
                bcData->updateInformationObject( p["ground"][bcId.second] );
        }
    }
    if ( !M_surfaceChargeDensity.empty() )
    {
        nl::json & pBC = p["surface_charge_density"];
        for ( auto const& [bcname,bcData] : M_surfaceChargeDensity )
            bcData->updateInformationObject( pBC[bcname] );
    }
}

tabulate_informations_ptr_t
ElectricBoundaryConditions::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains( "electric_potential_imposed" ) )
    {
        auto tabInfoElectricPotentialImposed = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "electric_potential_imposed" ).items() )
            tabInfoElectricPotentialImposed->add( j_bckey, ElectricBoundaryConditions::ElectricPotentialImposed::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Electric Potential Imposed", tabInfoElectricPotentialImposed );
    }
    if ( jsonInfo.contains( "ground" ) )
    {
        auto tabInfoGround = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "ground" ).items() )
            tabInfoGround->add( j_bckey, ElectricBoundaryConditions::ElectricPotentialImposed::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Ground", tabInfoGround );
    }
    if ( jsonInfo.contains( "surface_charge_density" ) )
    {
        auto tabInfoSurfaceCharge = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "surface_charge_density" ).items() )
            tabInfoSurfaceCharge->add( j_bckey, ElectricBoundaryConditions::SurfaceChargeDensity::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Surface Charge Density", tabInfoSurfaceCharge );
    }
    return tabInfo;
}

} // namespace FeelModels
} // namespace Feel
