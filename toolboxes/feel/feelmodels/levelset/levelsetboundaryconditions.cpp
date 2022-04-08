/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/levelset/levelsetboundaryconditions.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel
{
namespace FeelModels
{

void
LevelSetBoundaryConditions::setup( nl::json const& jarg )
{
    auto tbParent = this->toolboxParent();
    ModelIndexes indexes;
    //TODO
}

void
LevelSetBoundaryConditions::setParameterValues( std::map<std::string,double> const& paramValues )
{
    //TODO
}

void
LevelSetBoundaryConditions::updateInformationObject( nl::json & p ) const
{
    //TODO
}

tabulate_informations_ptr_t
LevelSetBoundaryConditions::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    //TODO
    return tabInfo;
}

} // namespace FeelModels
} // namespace Feel
