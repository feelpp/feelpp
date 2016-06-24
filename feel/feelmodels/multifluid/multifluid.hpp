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

    //--------------------------------------------------------------------//
    // Function spaces
    typedef typename levelset_type::space_levelset_ptrtype space_levelset_ptrtype;
    typedef typename levelset_type::space_levelset_vectorial_ptrtype space_levelset_vectorial_ptrtype;
    typedef typename levelset_type::space_markers_ptrtype space_levelset_markers_ptrtype;

    typedef typename levelset_type::element_levelset_type element_levelset_type;
    typedef typename levelset_type::element_levelset_ptrtype element_levelset_ptrtype; 
    typedef typename levelset_type::element_levelset_vectorial_type element_levelset_vectorial_type;
    typedef typename levelset_type::element_levelset_vectorial_ptrtype element_levelset_vectorial_ptrtype; 
    
    //--------------------------------------------------------------------//
    // Mesh
    typedef typename fluid_type::convex_type convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //--------------------------------------------------------------------//
    // Density/viscosity
    typedef typename fluid_type::densityviscosity_model_type densityviscosity_model_type;
    typedef typename fluid_type::densityviscosity_model_ptrtype densityviscosity_model_ptrtype;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//

    //--------------------------------------------------------------------//
    // Constructor
    MultiFluid(
            std::string const& prefix,
            WorldComm const& wc = Environment::worldComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );
    MultiFluid( self_type const& M ) = default;

    static self_ptrtype New(
            std::string const& prefix,
            WorldComm const& wc = Environment::worldComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );

    static std::string expandStringFromSpec( std::string const& expr );

    //--------------------------------------------------------------------//
    // Initialization
    void build( uint16_type nLevelSets );

    void init();

    virtual void loadParametersFromOptionsVm();

    //--------------------------------------------------------------------//
    // Function spaces
    space_levelset_ptrtype const& functionSpaceLevelset() const { return M_globalLevelset->functionSpace(); }
    space_levelset_vectorial_ptrtype const& functionSpaceLevelsetVectorial() const { return M_globalLevelset->functionSpaceVectorial(); }
    space_levelset_markers_ptrtype const& functionSpaceLevelsetMarkers() const { return M_globalLevelset->functionSpaceMarkers(); }
    //--------------------------------------------------------------------//
    // Mesh
    mesh_ptrtype const& mesh() const { return M_fluid->mesh(); }
    //--------------------------------------------------------------------//
    bool hasSurfaceTension() const { return M_enableSurfaceTension; }
    bool hasInterfaceForces() const;

    //--------------------------------------------------------------------//
    // Solve
    void solve();
    //--------------------------------------------------------------------//
    // Time step
    void updateTimeStep();

protected:
    void updateGlobalLevelset();

    void updateFluidDensityViscosity();
    void updateInterfaceForces();
    void advectLevelsets();

private:
    //--------------------------------------------------------------------//
    fluid_ptrtype M_fluid;
    levelset_ptrtype M_globalLevelset;
    std::vector<levelset_ptrtype> M_levelsets;

    uint16_type M_nFluids;

    //--------------------------------------------------------------------//
    // Parameters
    densityviscosity_model_ptrtype M_fluidDensityViscosityModel;
    std::vector<densityviscosity_model_ptrtype> M_levelsetDensityViscosityModels;
    //--------------------------------------------------------------------//
    // Forces
    bool M_enableSurfaceTension;
    ublas::symmetric_matrix<double, ublas::upper> M_surfaceTensionCoeff;

    element_levelset_vectorial_ptrtype M_interfaceForces; 
};
        

} // namespace FeelModels
} // namespace Feel

#endif
