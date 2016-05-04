
#include <feel/feelmodels/advection/advection.hpp>

namespace Feel {
namespace FeelModels {

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
ADVECTION_CLASS_TEMPLATE_TYPE::Advection( 
        space_advection_ptrtype const& space,
        std::string const& prefix,
        WorldComm const& worldComm,
        std::string const& subPrefix,
        std::string const& rootRepository )
: super_type(space, prefix, worldComm, subPrefix, rootRepository)
{
}

ADVECTION_CLASS_TEMPLATE_DECLARATIONS
void
ADVECTION_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );
}

} // namespace FeelModels
} // namespace Feel
