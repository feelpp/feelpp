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

template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::Inlet::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
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

     if ( jarg.contains( "shape" ) )
     {
         std::string shapeStr = jarg.at( "shape" ).template get<std::string>();
         if ( shapeStr == "constant" )
             M_shape = Shape::constant;
         else if ( shapeStr == "parabolic" )
             M_shape = Shape::parabolic;
     }

     if ( jarg.contains( "constraint" ) )
     {
         std::string constraintStr = jarg.at( "constraint" ).template get<std::string>();
         if ( constraintStr == "velocity_max" )
             M_constraint = Constraint::velocity_max;
         else if ( constraintStr == "flow_rate" )
             M_constraint = Constraint::flow_rate;
     }
}

template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::Inlet::updateInformationObject( nl::json & p ) const
{
    auto [exprStr,compInfo] = M_mexpr.exprInformations();
    p["expr"] = exprStr;
    p["markers"] = M_markers;

    std::map<Shape,std::string> shapeToStr = { { Shape::constant, "constant" }, { Shape::parabolic, "parabolic" } };
    p["shape"] = shapeToStr[M_shape];
    std::map<Constraint,std::string> constraintToStr = { { Constraint::velocity_max, "velocity_max" }, { Constraint::flow_rate, "flow_rate" } };
    p["constraint"] = constraintToStr[M_constraint];
}

template <uint16_type Dim>
tabulate_informations_ptr_t
FluidMechanicsBoundaryConditions<Dim>::Inlet::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
    TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { /*"name","type",*/"shape","constraint","expr" } );
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
FluidMechanicsBoundaryConditions<Dim>::PressureImposed::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
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

template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::PressureImposed::updateInformationObject( nl::json & p ) const
{
    auto [exprStr,compInfo] = M_mexpr.exprInformations();
    p["expr"] = exprStr;
    p["markers"] = M_markers;
}

template <uint16_type Dim>
tabulate_informations_ptr_t
FluidMechanicsBoundaryConditions<Dim>::PressureImposed::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
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

template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::BodyInterface::setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes )
{
    // if ( jarg.contains( "expr" ) )
    //     M_mexpr.setExpr( jarg.at( "expr" ), mparent.worldComm(), mparent.repository().expr(), indexes );
     if ( jarg.contains( "markers" ) )
     {
         ModelMarkers markers;
         markers.setup( jarg.at("markers"), indexes );
         M_markers = markers;
     }
     else
         M_markers = { M_name };

     if ( jarg.contains("materials") )
         M_jsonMaterials = jarg.at("materials");

      if ( jarg.contains( "translational-velocity" ) )
          M_mexprTranslationalVelocity.setExpr( jarg.at( "translational-velocity"), mparent.worldComm(), mparent.repository().expr(), indexes );
      if ( jarg.contains( "angular-velocity" ) )
          M_mexprAngularVelocity.setExpr( jarg.at( "angular-velocity"), mparent.worldComm(), mparent.repository().expr(), indexes );
}

template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::BodyInterface::updateInformationObject( nl::json & p ) const
{
    // auto [exprStr,compInfo] = M_mexpr.exprInformations();
    // p["expr"] = exprStr;
    p["markers"] = M_markers;
}

template <uint16_type Dim>
tabulate_informations_ptr_t
FluidMechanicsBoundaryConditions<Dim>::BodyInterface::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
#if 0
    TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { /*"name","type",*/"expr" } );
#endif
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
    for ( std::string const& bcKeyword : { "inlet" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<Inlet>( j_bckey );
                bc->setup( mparent,j_bcval,indexes );
                M_inlet.emplace( j_bckey, std::move( bc ) );
            }
        }
    }
    for ( std::string const& bcKeyword : { "pressure_imposed", "pressure" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<PressureImposed>( j_bckey );
                bc->setup( mparent,j_bcval,indexes );
                M_pressureImposed.emplace( j_bckey, std::move( bc ) );
            }
        }
    }
    for ( std::string const& bcKeyword : { "body_interface", "body" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<BodyInterface>( j_bckey );
                bc->setup( mparent,j_bcval,indexes );
                M_bodyInterface.emplace( j_bckey, std::move( bc ) );
            }
        }
    }
}

template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::setParameterValues( std::map<std::string,double> const& paramValues )
{
    for ( auto & [bcId,bcData] : M_velocityImposed )
        bcData->setParameterValues( paramValues );
    for ( auto & [bcName,bcData] : M_inlet )
        bcData->setParameterValues( paramValues );
    for ( auto & [bcName,bcData] : M_pressureImposed )
        bcData->setParameterValues( paramValues );
    for ( auto & [bcName,bcData] : M_bodyInterface )
        bcData->setParameterValues( paramValues );
}

template <uint16_type Dim>
void
FluidMechanicsBoundaryConditions<Dim>::updateInformationObject( nl::json & p ) const
{
    for ( auto const& [bcId,bcData] : M_velocityImposed )
    {
        if ( bcId.first == Type::VelocityImposed )
            bcData->updateInformationObject( p["velocity_imposed"][bcId.second] );
        else if ( bcId.first == Type::MeshVelocityImposed )
            bcData->updateInformationObject( p["mesh_velocity_imposed"][bcId.second] );
    }
    for ( auto const& [bcName,bcData] : M_inlet )
        bcData->updateInformationObject( p["inlet"][bcName] );
    for ( auto const& [bcName,bcData] : M_pressureImposed )
        bcData->updateInformationObject( p["pressure_imposed"][bcName] );
    for ( auto const& [bcName,bcData] : M_bodyInterface )
        bcData->updateInformationObject( p["body_interface"][bcName] );
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
    if ( jsonInfo.contains( "inlet" ) )
    {
        auto tabInfoBC = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "inlet" ).items() )
            tabInfoBC->add( j_bckey, Inlet::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Inlet", tabInfoBC );
    }
    if ( jsonInfo.contains( "pressure_imposed" ) )
    {
        auto tabInfoBC = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "pressure_imposed" ).items() )
            tabInfoBC->add( j_bckey, PressureImposed::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Pressure Imposed", tabInfoBC );
    }
    if ( jsonInfo.contains( "body_interface" ) )
    {
        auto tabInfoBC = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "body_interface" ).items() )
            tabInfoBC->add( j_bckey, BodyInterface::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Body Interface", tabInfoBC );
    }
    return tabInfo;
}


template class FluidMechanicsBoundaryConditions<2>;
template class FluidMechanicsBoundaryConditions<3>;

} // namespace FeelModels
} // namespace Feel
