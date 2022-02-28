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
SolidMechanicsBoundaryConditions<Dim>::setup( ModelBase const& mparent, nl::json const& jarg )
{
    ModelIndexes indexes;
    for ( std::string const& bcKeyword : { "displacement_imposed", "displacement" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_bcval] : j_bc.items() )
            {
                auto bc = std::make_shared<DisplacementImposed>( j_bckey );
                bc->setup( mparent,j_bcval,indexes );
                M_displacementImposed.emplace( std::make_pair(Type::DisplacementImposed, j_bckey), std::move( bc ) );
            }
        }
    }
#if 0
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
#endif
}

template <uint16_type Dim>
void
SolidMechanicsBoundaryConditions<Dim>::setParameterValues( std::map<std::string,double> const& paramValues )
{
    for ( auto & [bcId,bcData] : M_displacementImposed )
        bcData->setParameterValues( paramValues );
    // for ( auto & [bcName,bcData] : M_inlet )
    //     bcData->setParameterValues( paramValues );
    // for ( auto & [bcName,bcData] : M_pressureImposed )
    //     bcData->setParameterValues( paramValues );
    // for ( auto & [bcName,bcData] : M_bodyInterface )
    //     bcData->setParameterValues( paramValues );
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
#if 0
    for ( auto const& [bcName,bcData] : M_inlet )
        bcData->updateInformationObject( p["inlet"][bcName] );
    for ( auto const& [bcName,bcData] : M_pressureImposed )
        bcData->updateInformationObject( p["pressure_imposed"][bcName] );
    for ( auto const& [bcName,bcData] : M_bodyInterface )
        bcData->updateInformationObject( p["body_interface"][bcName] );
#endif
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
#if 0
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
#endif
    return tabInfo;
}


template class SolidMechanicsBoundaryConditions<2>;
template class SolidMechanicsBoundaryConditions<3>;

} // namespace FeelModels
} // namespace Feel
