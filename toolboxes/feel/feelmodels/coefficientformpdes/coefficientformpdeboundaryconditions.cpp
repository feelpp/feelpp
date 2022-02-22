/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpdeboundaryconditions.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim,uint8_type EquationRank>
void
CoefficientFormPDEBoundaryConditions<Dim,EquationRank>::Neumann::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
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

template <uint16_type Dim,uint8_type EquationRank>
void
CoefficientFormPDEBoundaryConditions<Dim,EquationRank>::Neumann::updateInformationObject( nl::json & p ) const
{
    auto [exprStr,compInfo] = M_mexpr.exprInformations();
    p["expr"] = exprStr;
    p["markers"] = M_markers;
}

template <uint16_type Dim,uint8_type EquationRank>
tabulate_informations_ptr_t
CoefficientFormPDEBoundaryConditions<Dim,EquationRank>::Neumann::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
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

template <uint16_type Dim,uint8_type EquationRank>
void
CoefficientFormPDEBoundaryConditions<Dim,EquationRank>::Robin::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
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

template <uint16_type Dim,uint8_type EquationRank>
void
CoefficientFormPDEBoundaryConditions<Dim,EquationRank>::Robin::updateInformationObject( nl::json & p ) const
{
    auto [exprStr1,compInfo1] = M_mexpr1.exprInformations();
    p["expr1"] = exprStr1;
    auto [exprStr2,compInfo2] = M_mexpr2.exprInformations();
    p["expr2"] = exprStr2;
    p["markers"] = M_markers;
}

template <uint16_type Dim,uint8_type EquationRank>
tabulate_informations_ptr_t
CoefficientFormPDEBoundaryConditions<Dim,EquationRank>::Robin::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
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


template <uint16_type Dim,uint8_type EquationRank>
void
CoefficientFormPDEBoundaryConditions<Dim,EquationRank>::setup( ModelBase const& mparent, nl::json const& jarg )
{

    ModelIndexes indexes;
    for ( std::string const& bcKeyword : { "Dirichlet" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<Dirichlet>( j_bckey );
                bc->setup( mparent,j_bcval,indexes );
                M_dirichlet.emplace( j_bckey, std::move( bc ) );
            }
        }
    }
    for ( std::string const& bcKeyword : { "Neumann" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<Neumann>( j_bckey );
                bc->setup( mparent,j_bcval,indexes );
                M_neumann.emplace( j_bckey, std::move( bc ) );
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
                bc->setup( mparent,j_bcval,indexes );
                M_robin.emplace(j_bckey, std::move( bc ) );
            }
        }
    }
}

template <uint16_type Dim,uint8_type EquationRank>
void
CoefficientFormPDEBoundaryConditions<Dim,EquationRank>::setParameterValues( std::map<std::string,double> const& paramValues )
{

    for ( auto & [bcName,bcData] : M_dirichlet )
        bcData->setParameterValues( paramValues );
    for ( auto & [bcName,bcData] : M_neumann )
        bcData->setParameterValues( paramValues );
    for ( auto & [bcName,bcData] : M_robin )
        bcData->setParameterValues( paramValues );
}

template <uint16_type Dim,uint8_type EquationRank>
void
CoefficientFormPDEBoundaryConditions<Dim,EquationRank>::updateInformationObject( nl::json & p ) const
{
    if ( !M_dirichlet.empty() )
    {
        nl::json & pBC = p["Dirichlet"];
        for ( auto const& [bcName,bcData] : M_dirichlet )
            bcData->updateInformationObject( pBC[bcName] );
    }
    if ( !M_neumann.empty() )
    {
        nl::json & pBC = p["Neumann"];
        for ( auto const& [bcName,bcData] : M_neumann )
            bcData->updateInformationObject( pBC[bcName] );
    }
    if ( !M_robin.empty() )
    {
        nl::json & pBC = p["Robin"];
        for ( auto const& [bcName,bcData] : M_robin )
            bcData->updateInformationObject( pBC[bcName] );
    }
}

template <uint16_type Dim,uint8_type EquationRank>
tabulate_informations_ptr_t
CoefficientFormPDEBoundaryConditions<Dim,EquationRank>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );

    if ( jsonInfo.contains( "Dirichlet" ) )
    {
        auto tabInfoDirichlet = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "Dirichlet" ).items() )
            tabInfoDirichlet->add( j_bckey, Dirichlet::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Dirichlet", tabInfoDirichlet );
    }
    if ( jsonInfo.contains( "Neumann" ) )
    {
        auto tabInfoNeumann = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "Neumann" ).items() )
            tabInfoNeumann->add( j_bckey, Neumann::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Neumann", tabInfoNeumann );
    }
    if ( jsonInfo.contains( "Robin" ) )
    {
        auto tabInfoRobin = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "Robin" ).items() )
            tabInfoRobin->add( j_bckey, Robin::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Robin", tabInfoRobin );
    }
    return tabInfo;
}

template class CoefficientFormPDEBoundaryConditions<2,0>;
template class CoefficientFormPDEBoundaryConditions<2,1>;
template class CoefficientFormPDEBoundaryConditions<3,0>;
template class CoefficientFormPDEBoundaryConditions<3,1>;


} // namespace FeelModels
} // namespace Feel
