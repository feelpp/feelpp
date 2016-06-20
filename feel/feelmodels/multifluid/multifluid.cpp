/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/multifluid.hpp>

namespace Feel {
namespace FeelModels {

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
MULTIFLUID_CLASS_TEMPLATE_TYPE::MultiFluid(
        std::string const& prefix,
        WorldComm const& wc,
        std::string const& subPrefix,
        std::string const& rootRepository )
: super_type( prefix, wc, subPrefix, self_type::expandStringFromSpec( rootRepository ) ) 
{}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
std::string
MULTIFLUID_CLASS_TEMPLATE_TYPE::expandStringFromSpec( std::string const& s )
{
    std::string res = s;
    res = fluid_type::expandStringFromSpec( res );
    return res;
}

MULTIFLUID_CLASS_TEMPLATE_DECLARATIONS
void
MULTIFLUID_CLASS_TEMPLATE_TYPE::build( uint16_type nLevelSets )
{
    this->log("MultiFluid", "build", "start");

    M_fluid = fluid_ptrtype( 
            new fluid_type("fluid", false, this->worldComm(), "", this->rootRepositoryWithoutNumProc() ) 
            ); 
    M_fluid->build();

    M_levelsets.resize( nLevelSets );
    for( uint16_type i = 0; i < M_levelsets.size(); ++i )
    {
        M_levelsets[i] = levelset_ptrtype(
                new levelset_type( (boost::format("levelset%1%")%i).str(), this->worldComm(), "", this->rootRepositoryWithoutNumProc() )
                );
        M_levelsets[i]->build( M_fluid->mesh() );
    }

    this->log("MultiFluid", "build", "finish");
}

} // namespace FeelModels
} // namespace Feel
