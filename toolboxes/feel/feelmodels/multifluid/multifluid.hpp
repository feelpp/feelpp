#ifndef _MULTIFLUID_HPP
#define _MULTIFLUID_HPP 1

#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelmodels/levelset/levelset.hpp>
#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>

namespace Feel {
namespace FeelModels {

template< typename FluidType, typename LevelSetType>
class MultiFluid 
: public FluidType
{
public:
    // Typedefs
    //--------------------------------------------------------------------//
    // Class
    typedef FluidType super_type;
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
    typedef typename fluid_type::component_space_fluid_velocity_type component_space_fluid_velocity_type;

    typedef typename levelset_type::element_levelset_type element_levelset_type;
    typedef typename levelset_type::element_levelset_ptrtype element_levelset_ptrtype; 
    typedef typename levelset_type::element_levelset_vectorial_type element_levelset_vectorial_type;
    typedef typename levelset_type::element_levelset_vectorial_ptrtype element_levelset_vectorial_ptrtype; 

    //--------------------------------------------------------------------//
    typedef typename super_type::DataUpdateJacobian DataUpdateJacobian;
    typedef typename super_type::DataUpdateResidual DataUpdateResidual;
    typedef typename super_type::DataUpdateLinear DataUpdateLinear;
    
    //--------------------------------------------------------------------//
    // Mesh
    typedef typename fluid_type::convex_type convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    //--------------------------------------------------------------------//
    // Lagrange P1 iso-Pn
    typedef OperatorLagrangeP1<component_space_fluid_velocity_type> op_lagrangeP1_type;
    typedef boost::shared_ptr<op_lagrangeP1_type> op_lagrangeP1_ptrtype;

    //--------------------------------------------------------------------//
    // Density/viscosity
    typedef typename fluid_type::densityviscosity_model_type densityviscosity_model_type;
    typedef typename fluid_type::densityviscosity_model_ptrtype densityviscosity_model_ptrtype;
    //--------------------------------------------------------------------//
    // Interface forces model
    typedef InterfaceForcesModel<levelset_type> interfaceforces_model_type;
    typedef boost::shared_ptr<interfaceforces_model_type> interfaceforces_model_ptrtype;
    typedef Singleton<Feel::Factory<interfaceforces_model_type, std::string>> interfaceforces_factory_type;

    //--------------------------------------------------------------------//
    // Inextensibility
    typedef typename super_type::basis_fluid_p_type basis_fluid_p_type;
    typedef FunctionSpace< mesh_type, bases<basis_fluid_p_type> > space_inextensibilitylm_type;
    typedef boost::shared_ptr<space_inextensibilitylm_type> space_inextensibilitylm_ptrtype;

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
    void build();
    mesh_ptrtype createMesh();

    void init();

    void loadParametersFromOptionsVm();

    std::string const& prefix() const { return M_prefix; }
    std::string const& fluidPrefix() const { return super_type::prefix(); }

    boost::shared_ptr<std::ostringstream> getInfo() const override;

    //--------------------------------------------------------------------//
    boost::shared_ptr<self_type> shared_from_this() { return boost::static_pointer_cast<self_type>(super_type::shared_from_this()); }
    boost::shared_ptr<self_type const> shared_from_this() const { return boost::static_pointer_cast<self_type const>(super_type::shared_from_this()); }

    //--------------------------------------------------------------------//
    // Function spaces
    space_levelset_ptrtype const& functionSpaceLevelset() const { return M_globalLevelset->functionSpace(); }
    space_levelset_vectorial_ptrtype const& functionSpaceLevelsetVectorial() const { return M_globalLevelset->functionSpaceVectorial(); }
    space_levelset_markers_ptrtype const& functionSpaceLevelsetMarkers() const { return M_globalLevelset->functionSpaceMarkers(); }
    space_inextensibilitylm_ptrtype const& functionSpaceInextensibilityLM() const;
    //--------------------------------------------------------------------//
    // Mesh
    //mesh_ptrtype const& mesh() const { return M_mesh; }
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"MultiFluidMesh.path"); }
    //--------------------------------------------------------------------//
    // Models
    fluid_ptrtype fluidModel() { return this->shared_from_this(); }
    levelset_ptrtype const& levelsetModel(uint16_type n) const { return M_levelsets.at(n); }
    uint16_type nLevelsets() const { return M_levelsets.size(); }

    levelset_ptrtype const& globalLevelset() const { return M_globalLevelset; }

    densityviscosity_model_ptrtype const& fluidDensityViscosityModel( uint16_type n = 0 ) const;

    //--------------------------------------------------------------------//
    // Algebraic data
    int nBlockMatrixGraph() const override;
    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    size_type nLocalDof() const override;

    //--------------------------------------------------------------------//
    bool hasInextensibility( uint16_type n = 0 ) const { return M_hasInextensibility.at(n); }
    std::string const& inextensibilityMethod( uint16_type n = 0 ) const { return M_inextensibilityMethod.at(n); }
    //--------------------------------------------------------------------//
    bool hasInterfaceForces() const;

    void addInterfaceForce( interfaceforces_model_ptrtype model, std::string const& name = "" );
    interfaceforces_model_ptrtype const& interfaceForce( std::string const& name ) const;
    std::map<std::string, interfaceforces_model_ptrtype> const& interfaceForces() const;

    //--------------------------------------------------------------------//
    // Solve
    void solve();

    //--------------------------------------------------------------------//
    // Time step
    //boost::shared_ptr<TSBase> timeStepBase() const { return M_fluid->timeStepBase(); }
    boost::shared_ptr<TSBase> fluidTimeStepBase() const { return this->timeStepBase(); }
    boost::shared_ptr<TSBase> levelsetTimeStepBase(uint16_type n) const { return this->levelsetModel(n)->timeStepBase(); }
    void updateTime( double time );
    void updateTimeStep();

    //--------------------------------------------------------------------//
    // Export
    void exportResults( double time ) { this->exportResultsImpl( time ); }
    void exportResults() { this->exportResults( this->currentTime() ); }

protected:
    size_type initStartBlockIndexFieldsInMatrix() override;
    int initBlockVector() override;

    void updateGlobalLevelset();

    void updateFluidDensityViscosity();
    void updateInterfaceForces();
    void solveFluid();
    void advectLevelsets();
    virtual void solveImpl();

    void setRebuildMatrixVector( bool b = true ) { M_doRebuildMatrixVector = b; }
    bool rebuildMatrixVector() const { return M_doRebuildMatrixVector; }
    // Linear solve
    void updateLinearPDEAdditional( DataUpdateLinear & data ) const override;
    // Non-linear solve
    void updateJacobianAdditional( DataUpdateJacobian & data ) const override;
    void updateResidualAdditional( DataUpdateResidual & data ) const override;
    //--------------------------------------------------------------------//
    // Export
    virtual void exportResultsImpl( double time );

    //--------------------------------------------------------------------//
    uint16_type M_nFluids;

private:
    std::string M_prefix;
    //--------------------------------------------------------------------//
    // Lagrange P1 iso-Pn
    bool M_useLagrangeP1iso;
    op_lagrangeP1_ptrtype M_opLagrangeP1iso;
    //--------------------------------------------------------------------//
    //mesh_ptrtype M_mesh;
    levelset_ptrtype M_globalLevelset;
    std::vector<levelset_ptrtype> M_levelsets;

    //--------------------------------------------------------------------//
    // Solve
    bool M_doRebuildMatrixVector;
    bool M_usePicardIterations;

    //--------------------------------------------------------------------//
    // Parameters
    densityviscosity_model_ptrtype M_fluidDensityViscosityModel;
    std::vector<densityviscosity_model_ptrtype> M_levelsetDensityViscosityModels;
    std::vector<std::map<std::string, interfaceforces_model_ptrtype>> M_levelsetInterfaceForcesModels;
    std::map<std::string, interfaceforces_model_ptrtype> M_additionalInterfaceForcesModel;
    //--------------------------------------------------------------------//
    // Forces
    bool M_hasInterfaceForcesModel;

    element_levelset_vectorial_ptrtype M_interfaceForces; 

    //--------------------------------------------------------------------//
    // Inextensibility
    std::vector<bool> M_hasInextensibility;
    bool M_enableInextensibility;
    std::vector<std::string> M_inextensibilityMethod;

    mutable space_inextensibilitylm_ptrtype M_spaceInextensibilityLM;
    mutable bool M_doRebuildSpaceInextensibilityLM;
    // Penalty method gamma
    std::vector<double> M_inextensibilityGamma;
    //--------------------------------------------------------------------//
    // Reinitialization
    std::vector<int> M_levelsetReinitEvery;
    std::vector<int> M_levelsetReinitSmoothEvery;
};
        

} // namespace FeelModels
} // namespace Feel

#endif
