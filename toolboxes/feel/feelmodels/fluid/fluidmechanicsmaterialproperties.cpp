/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmechanicsmaterialproperties.hpp>

namespace Feel
{
namespace FeelModels
{

void
DynamicViscosityLaw::setup( nl::json const& jarg )
{
    bool isMultifluid = false;

    if ( jarg.is_string() )
    {
        this->setLaw( jarg.template get<std::string>() );
    }
    else if ( jarg.is_object() )
    {
        CHECK( jarg.contains( "law" ) && jarg.at( "law" ).is_string() ) << "invalid law";
        this->setLaw( jarg.at( "law" ).template get<std::string>() );
    }
}

const bool dynamicviscosity_law_register = ModelMaterial::ModelMaterialLawFactory::instance().registerProduct(
        "dynamic-viscosity",
        []() { return std::make_unique<DynamicViscosityLaw>( "newtonian" ); }
        );

}
}
