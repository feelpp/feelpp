#ifndef _ADVECTION_HPP
#define _ADVECTION_HPP 1

#include <feel/feelmodels/advection/advectionbase.hpp>

namespace Feel {
namespace FeelModels {

template< typename ConvexType, typename BasisAdvectionType >
class Advection
    : public AdvectionBase<ConvexType, BasisAdvectionType>
    , public boost::enable_shared_from_this< Advection<ConvexType, BasisAdvectionType> >
{
public:
    typedef AdvectionBase<ConvexType, BasisAdvectionType> super_type;

    typedef Advection<ConvexType, BasisAdvectionType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef typename super_type::space_advection_ptrtype space_advection_ptrtype;

    //--------------------------------------------------------------------//
    // Constructor
    Advection( 
            space_advection_ptrtype const& space,
            std::string const& prefix,
            WorldComm const& _worldComm = Environment::worldComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );
    //--------------------------------------------------------------------//
    // Initialization
    void init( bool buildModelAlgebraicFactory = true );
};
    

} // namespace FeelModels
} // namespace Feel

#endif
