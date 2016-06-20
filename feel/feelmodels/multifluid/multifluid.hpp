#ifndef _MULTIFLUID_HPP
#define _MULTIFLUID_HPP 1

#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelmodels/levelset/levelset.hpp>

namespace Feel {
namespace FeelModels {

template< typename FluidType, typename LevelSetType>
class MultiFluid : public ModelNumerical
{
public:
    // Typedefs
    //--------------------------------------------------------------------//
    // Class
    typedef ModelNumerical super_type;
    typedef MultiFluid< FluidType, LevelSetType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef FluidType fluid_type;
    typedef boost::shared_ptr<fluid_type> fluid_ptrtype;

    typedef LevelSetType levelset_type;
    typedef boost::shared_ptr<levelset_type> levelset_ptrtype;

    typedef std::vector<levelset_ptrtype> vector_levelset_ptrtype;
    
    //--------------------------------------------------------------------//
    // Mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//

    //--------------------------------------------------------------------//
    // Constructor
    MultiFluid(
            std::string const& prefix,
            WorldComm const& wc = Environment::wordlComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );
    MultiFluid( self_type const& M ) = default;

    static self_ptrtype New(
            std::string const& prefix,
            WorldComm const& wc = Environment::wordlComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );

    static std::string expandStringFromSpec( std::string const& expr );

    //--------------------------------------------------------------------//
    // Initialization
    void build( uint16_type nLevelSets );

    void init();

    //--------------------------------------------------------------------//
    // Solve
    void solve();

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private:
    
    fluid_ptrtype M_fluid;
    vector_levelset_ptrtype M_levelsets;

};
        

} // namespace FeelModels
} // namespace Feel

#endif
