/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/solid/solidmechanics1dreducedboundaryconditions.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel
{
namespace FeelModels
{
template <uint16_type Dim>
void
SolidMechanics1dReducedBoundaryConditions<Dim>::setup( nl::json const& jarg )
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
}

template <uint16_type Dim>
void
SolidMechanics1dReducedBoundaryConditions<Dim>::setParameterValues( std::map<std::string,double> const& paramValues )
{
    for ( auto & [bcId,bcData] : M_displacementImposed )
        bcData->setParameterValues( paramValues );
}

template <uint16_type Dim>
void
SolidMechanics1dReducedBoundaryConditions<Dim>::updateInformationObject( nl::json & p ) const
{
    for ( auto const& [bcId,bcData] : M_displacementImposed )
    {
        if ( bcId.first == Type::DisplacementImposed )
            bcData->updateInformationObject( p["displacement_imposed"][bcId.second] );
    }
}

template <uint16_type Dim>
tabulate_informations_ptr_t
SolidMechanics1dReducedBoundaryConditions<Dim>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );

    if ( jsonInfo.contains( "displacement_imposed" ) )
    {
        auto tabInfoBC = TabulateInformationsSections::New( tabInfoProp );
        for ( auto const& [j_bckey,j_bcval]: jsonInfo.at( "displacement_imposed" ).items() )
            tabInfoBC->add( j_bckey, DisplacementImposed::tabulateInformations( j_bcval, tabInfoProp ) );
        tabInfo->add( "Displacement Imposed", tabInfoBC );
    }
    return tabInfo;
}


template class SolidMechanics1dReducedBoundaryConditions<2>;
template class SolidMechanics1dReducedBoundaryConditions<3>;

} // namespace FeelModels
} // namespace Feel
