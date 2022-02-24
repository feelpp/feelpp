/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmechanicsboundaryconditions.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel
{
namespace FeelModels
{
template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::MeshVelocityImposed::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
{
    nl::json j_bc = jarg;
    if ( Dim == 2 )
        j_bc["expr"] = "{meshes_fluid_meshmotion_v_0,meshes_fluid_meshmotion_v_1}:meshes_fluid_meshmotion_v_0:meshes_fluid_meshmotion_v_1";
    else
        j_bc["expr"] = "{meshes_fluid_meshmotion_v_0,meshes_fluid_meshmotion_v_1,meshes_fluid_meshmotion_v_2}:meshes_fluid_meshmotion_v_0:meshes_fluid_meshmotion_v_1:meshes_fluid_meshmotion_v_2";
    super_type::setup( mparent,j_bc,indexes );
}
template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::MeshVelocityImposed::updateInformationObject( nl::json & p ) const
{
    // auto [exprStr,compInfo] = M_mexpr.exprInformations();
    // p["expr"] = exprStr;
    p["markers"] = this->M_markers;
}

template <uint16_type Dim>
tabulate_informations_ptr_t
FluidMechanicsBoundaryConditions<Dim>::MeshVelocityImposed::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
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

#if 0
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
#endif

template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::setup( ModelBase const& mparent, nl::json const& jarg )
{

    ModelIndexes indexes;
    for ( std::string const& bcKeyword : { "velocity_imposed", "velocity" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<VelocityImposed>( j_bckey );
                bc->setup( mparent,j_bcval,indexes );
                M_velocityImposed.emplace( std::make_pair(Type::VelocityImposed, j_bckey), std::move( bc ) );
            }
        }
    }
    for ( std::string const& bcKeyword : { "mesh_velocity_imposed", "mesh_velocity" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            std::string bcName = "";
            auto bc = std::make_shared<MeshVelocityImposed>( bcName );
            bc->setup( mparent,j_bc,indexes );
            M_velocityImposed.emplace( std::make_pair(Type::MeshVelocityImposed,bcName), std::move( bc ) );
        }
    }
}

template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::setParameterValues( std::map<std::string,double> const& paramValues )
{
    for ( auto & [bcId,bcData] : M_velocityImposed )
        bcData->setParameterValues( paramValues );
}

template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::updateInformationObject( nl::json & p ) const
{
    if ( !M_velocityImposed.empty() )
    {
        for ( auto const& [bcId,bcData] : M_velocityImposed )
        {
            if ( bcId.first == Type::VelocityImposed )
                bcData->updateInformationObject( p["velocity_imposed"][bcId.second] );
            else if ( bcId.first == Type::MeshVelocityImposed )
                bcData->updateInformationObject( p["mesh_velocity_imposed"][bcId.second] );
        }
    }
}

template <uint16_type Dim>
tabulate_informations_ptr_t
FluidMechanicsBoundaryConditions<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains( "velocity_imposed" ) )
    {
        auto tabInfoBC = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "velocity_imposed" ).items() )
            tabInfoBC->add( j_bckey, VelocityImposed::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Velocity Imposed", tabInfoBC );
    }
    if ( jsonInfo.contains( "mesh_velocity_imposed" ) )
    {
        auto tabInfoBC = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "mesh_velocity_imposed" ).items() )
            tabInfoBC->add( j_bckey, MeshVelocityImposed::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Mesh Velocity Imposed", tabInfoBC );
    }
    return tabInfo;
}


template class FluidMechanicsBoundaryConditions<2>;
template class FluidMechanicsBoundaryConditions<3>;

} // namespace FeelModels
} // namespace Feel
