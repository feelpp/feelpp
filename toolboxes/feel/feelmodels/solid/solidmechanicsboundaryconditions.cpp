/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/solid/solidmechanicsboundaryconditions.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel
{
namespace FeelModels
{
template <uint16_type Dim>
void
SolidMechanicsBoundaryConditions<Dim>::NormalStress::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
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

     if ( jarg.contains( "frame" ) )
     {
         std::string const& frameStr = jarg.at( "frame" ).template get<std::string>();
         if ( frameStr == "lagrangian" )
             M_frame = Frame::Lagrangian;
         else if ( frameStr == "eulerian" )
             M_frame = Frame::Eulerian;
     }
}

template <uint16_type Dim>
void
SolidMechanicsBoundaryConditions<Dim>::NormalStress::updateInformationObject( nl::json & p ) const
{
    auto [exprStr,compInfo] = M_mexpr.exprInformations();
    p["expr"] = exprStr;
    p["markers"] = M_markers;
    if ( this->isFrameLagrangian() )
        p["frame"] = "lagrangian";
    else
        p["frame"] = "eulerianlagrangian";
}

template <uint16_type Dim>
tabulate_informations_ptr_t
SolidMechanicsBoundaryConditions<Dim>::NormalStress::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
    TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { /*"name","type",*/"frame","expr" } );
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

template <uint16_type Dim>
void
SolidMechanicsBoundaryConditions<Dim>::Robin::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
{
    if ( jarg.contains( "expr1" ) )
        M_mexpr1.setExpr( jarg.at( "expr1" ), mparent.worldComm(), mparent.repository().expr(), indexes );
    if ( jarg.contains( "expr2" ) )
        M_mexpr2.setExpr( jarg.at( "expr2" ), mparent.worldComm(), mparent.repository().expr(), indexes );
     if ( jarg.contains( "markers" ) )
     {
         ModelMarkers markers;
         markers.setup( jarg.at("markers"), indexes );
         M_markers = markers;
     }
     else
         M_markers = { M_name };
}

template <uint16_type Dim>
void
SolidMechanicsBoundaryConditions<Dim>::Robin::updateInformationObject( nl::json & p ) const
{
    auto [exprStr1,compInfo1] = M_mexpr1.exprInformations();
    p["expr1"] = exprStr1;
    auto [exprStr2,compInfo2] = M_mexpr2.exprInformations();
    p["expr2"] = exprStr2;
    p["markers"] = M_markers;
}

template <uint16_type Dim>
tabulate_informations_ptr_t
SolidMechanicsBoundaryConditions<Dim>::Robin::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
    TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { /*"name","type",*/"expr1","expr2" } );
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


template <uint16_type Dim>
void
SolidMechanicsBoundaryConditions<Dim>::setup( nl::json const& jarg )
{
    auto tbParent = this->toolboxParent();
    ModelIndexes indexes;
    for ( std::string const& bcKeyword : { "displacement_imposed", "displacement" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<DisplacementImposed>( j_bckey, tbParent );
                bc->setup( j_bcval,indexes );
                M_displacementImposed.emplace( std::make_pair(Type::DisplacementImposed, j_bckey), std::move( bc ) );
            }
        }
    }

     // TODO maybe add follower_pressure

    for ( std::string const& bcKeyword : { "normal_stress" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<NormalStress>( j_bckey );
                bc->setup( *tbParent,j_bcval,indexes );
                M_normalStress.emplace( j_bckey, std::move( bc ) );
            }
        }
    }
    for ( std::string const& bcKeyword : { "Robin" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<Robin>( j_bckey );
                bc->setup( *tbParent,j_bcval,indexes );
                M_robin.emplace( j_bckey, std::move( bc ) );
            }
        }
    }
}

template <uint16_type Dim>
void
SolidMechanicsBoundaryConditions<Dim>::setParameterValues( std::map<std::string,double> const& paramValues )
{
    for ( auto & [bcId,bcData] : M_displacementImposed )
        bcData->setParameterValues( paramValues );
    for ( auto & [bcName,bcData] : M_normalStress )
        bcData->setParameterValues( paramValues );
    for ( auto & [bcName,bcData] : M_robin )
        bcData->setParameterValues( paramValues );
}

template <uint16_type Dim>
void
SolidMechanicsBoundaryConditions<Dim>::updateInformationObject( nl::json & p ) const
{
    for ( auto const& [bcId,bcData] : M_displacementImposed )
    {
        if ( bcId.first == Type::DisplacementImposed )
            bcData->updateInformationObject( p["displacement_imposed"][bcId.second] );
    }
    for ( auto const& [bcName,bcData] : M_normalStress )
        bcData->updateInformationObject( p["normal_stress"][bcName] );
    for ( auto const& [bcName,bcData] : M_robin )
        bcData->updateInformationObject( p["robin"][bcName] );
}

template <uint16_type Dim>
tabulate_informations_ptr_t
SolidMechanicsBoundaryConditions<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );

    if ( jsonInfo.contains( "displacement_imposed" ) )
    {
        auto tabInfoBC = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "displacement_imposed" ).items() )
            tabInfoBC->add( j_bckey, DisplacementImposed::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Displacement Imposed", tabInfoBC );
    }
    if ( jsonInfo.contains( "normal_stress" ) )
    {
        auto tabInfoBC = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "normal_stress" ).items() )
            tabInfoBC->add( j_bckey, NormalStress::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Normal Stress", tabInfoBC );
    }
    if ( jsonInfo.contains( "robin" ) )
    {
        auto tabInfoBC = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "robin" ).items() )
            tabInfoBC->add( j_bckey, Robin::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Robin", tabInfoBC );
    }
    return tabInfo;
}


template class SolidMechanicsBoundaryConditions<2>;
template class SolidMechanicsBoundaryConditions<3>;

} // namespace FeelModels
} // namespace Feel
