/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/solid/solidmechanics.hpp>

namespace Feel
{
namespace FeelModels
{

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    this->updateLinearPDE( data, this->modelContext() );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    this->updateLinearPDEDofElimination( data, this->modelContext() );
}

} // FeelModels

} // Feel


