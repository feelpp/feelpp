/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/magnetic/magneticboundaryconditions.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
void
MagneticBoundaryConditions<Dim>::setup( nl::json const& jarg )
{
    auto tbParent = this->toolboxParent();
    ModelIndexes indexes;
    for ( std::string const& bcKeyword : { "magnetic_potential_imposed", "magnetic_potential" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_temp = jarg.at( bcKeyword );
            for ( auto const& [j_tempkey,j_tempval] : j_temp.items() )
            {
                auto bc = std::make_shared<MagneticPotentialImposed>( j_tempkey, tbParent );
                bc->setup( j_tempval,indexes );
                M_magneticPotentialImposed.emplace(j_tempkey, std::move( bc ) );
            }
        }
    }

    for ( std::string const& bcKeyword : { "magnetic_insulation" } )
    {
        if ( jarg.contains( bcKeyword ) )
        {
            auto const& j_bc = jarg.at( bcKeyword );
            for ( auto const& [j_bckey,j_tempval] : j_bc.items() )
            {
                nl::json j_bcvalUpdated = j_tempval;
                j_bcvalUpdated["expr"] = "{0,0,0}";
                auto bc = std::make_shared<MagneticInsulation>( j_bckey, tbParent );
                bc->setup( j_bcvalUpdated,indexes );
                M_magneticInsulation.emplace(j_bckey, std::move( bc ) );
            }
        }
    }
}

template <uint16_type Dim>
void
MagneticBoundaryConditions<Dim>::setParameterValues( std::map<std::string,double> const& paramValues )
{
    for ( auto & [bcname,bcData] : M_magneticPotentialImposed )
        bcData->setParameterValues( paramValues );
    for ( auto & [bcname,bcData] : M_magneticInsulation )
        bcData->setParameterValues( paramValues );
}

template <uint16_type Dim>
void
MagneticBoundaryConditions<Dim>::updateInformationObject( nl::json & p ) const
{
    if ( !M_magneticPotentialImposed.empty() )
    {
        nl::json & pBC = p["magnetic_potential_imposed"];
        for ( auto const& [bcname,bcData] : M_magneticPotentialImposed )
            bcData->updateInformationObject( pBC[bcname] );
    }

    if ( !M_magneticInsulation.empty() )
    {
        nl::json & pBC = p["magnetic_insulation"];
        for ( auto const& [bcname,bcData] : M_magneticInsulation )
            bcData->updateInformationObject( pBC[bcname] );
    }
}

template <uint16_type Dim>
tabulate_informations_ptr_t
MagneticBoundaryConditions<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains( "magnetic_potential_imposed" ) )
    {
        auto tabInfoMagneticPotentialImposed = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "magnetic_potential_imposed" ).items() )
            tabInfoMagneticPotentialImposed->add( j_bckey, MagneticBoundaryConditions::MagneticPotentialImposed::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Magnetic Potential Imposed", tabInfoMagneticPotentialImposed );
    }
    if ( jsonInfo.contains( "magnetic_insulation" ) )
    {
        auto tabInfoMagneticInsulation = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "magnetic_insulation" ).items() )
            tabInfoMagneticInsulation->add( j_bckey, MagneticBoundaryConditions::MagneticInsulation::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Magnetic Insulation", tabInfoMagneticInsulation );
    }
    return tabInfo;
}

template class MagneticBoundaryConditions<2>;
template class MagneticBoundaryConditions<3>;

} // namespace FeelModels
} // namespace Feel
