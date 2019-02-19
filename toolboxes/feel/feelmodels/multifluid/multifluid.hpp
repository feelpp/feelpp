#ifndef _MULTIFLUID_HPP
#define _MULTIFLUID_HPP 1

#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelmodels/levelset/levelset.hpp>
#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>
#include <feel/feelmodels/levelset/globallevelsetexpr.hpp>

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
    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef FluidType fluid_type;
    typedef std::shared_ptr<fluid_type> fluid_ptrtype;

    typedef LevelSetType levelset_type;
    typedef std::shared_ptr<levelset_type> levelset_ptrtype;

    //--------------------------------------------------------------------//
    // Mesh
    typedef typename fluid_type::convex_type convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    //--------------------------------------------------------------------//
    // Range types
    typedef typename MeshTraits<mesh_type>::element_reference_wrapper_const_iterator element_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
    typedef elements_reference_wrapper_t<mesh_type> range_elements_type;

    typedef typename MeshTraits<mesh_type>::face_reference_wrapper_const_iterator face_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_type faces_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype faces_reference_wrapper_ptrtype;
    typedef faces_reference_wrapper_t<mesh_type> range_faces_type;

    //--------------------------------------------------------------------//
    // Function spaces
    typedef typename levelset_type::space_levelset_ptrtype space_levelset_ptrtype;
    typedef typename levelset_type::space_vectorial_ptrtype space_levelset_vectorial_ptrtype;
    typedef typename levelset_type::space_markers_ptrtype space_levelset_markers_ptrtype;
    typedef typename levelset_type::space_advection_velocity_type space_levelset_advection_velocity_type;
    typedef typename levelset_type::space_advection_velocity_ptrtype space_levelset_advection_velocity_ptrtype;
    typedef typename fluid_type::component_space_fluid_velocity_type component_space_fluid_velocity_type;

    typedef typename fluid_type::space_fluid_velocity_type space_fluid_velocity_type;
    typedef typename fluid_type::space_fluid_velocity_ptrtype space_fluid_velocity_ptrtype;

    typedef typename levelset_type::element_levelset_type element_levelset_type;
    typedef typename levelset_type::element_levelset_ptrtype element_levelset_ptrtype; 
    typedef typename levelset_type::element_vectorial_type element_levelset_vectorial_type;
    typedef typename levelset_type::element_vectorial_ptrtype element_levelset_vectorial_ptrtype; 

    // Levelset function space manager
    typedef typename levelset_type::levelset_space_manager_type levelset_space_manager_type;
    typedef typename levelset_type::levelset_space_manager_ptrtype levelset_space_manager_ptrtype;
    // Levelset tool manager
    typedef typename levelset_type::levelset_tool_manager_type levelset_tool_manager_type;
    typedef typename levelset_type::levelset_tool_manager_ptrtype levelset_tool_manager_ptrtype;

    //--------------------------------------------------------------------//
    typedef typename super_type::DataUpdateJacobian DataUpdateJacobian;
    typedef typename super_type::DataUpdateResidual DataUpdateResidual;
    typedef typename super_type::DataUpdateLinear DataUpdateLinear;
    
    //--------------------------------------------------------------------//
    // Lagrange P1 iso-Pn
    typedef OperatorLagrangeP1<component_space_fluid_velocity_type> op_lagrangeP1_type;
    typedef std::shared_ptr<op_lagrangeP1_type> op_lagrangeP1_ptrtype;

    //--------------------------------------------------------------------//
    // Density/viscosity
    typedef typename fluid_type::material_properties_type material_properties_type;
    typedef typename fluid_type::material_properties_ptrtype material_properties_ptrtype;
    //--------------------------------------------------------------------//
    // Interface forces model
    typedef InterfaceForcesModel<levelset_type, fluid_type> interfaceforces_model_type;
    typedef std::shared_ptr<interfaceforces_model_type> interfaceforces_model_ptrtype;
    typedef Singleton<Feel::Factory<interfaceforces_model_type, std::string>> interfaceforces_factory_type;

    //--------------------------------------------------------------------//
    // Inextensibility
    typedef typename super_type::basis_fluid_p_type basis_fluid_p_type;
    typedef FunctionSpace< mesh_type, bases<basis_fluid_p_type> > space_inextensibilitylm_type;
    typedef std::shared_ptr<space_inextensibilitylm_type> space_inextensibilitylm_ptrtype;
    //--------------------------------------------------------------------//
    // Exporter
    typedef Exporter<mesh_type, nOrderGeo> exporter_type;
    typedef std::shared_ptr<exporter_type> exporter_ptrtype;

    //--------------------------------------------------------------------//
    // Measures
    using force_type = Eigen::Matrix<typename super_type::value_type, nDim, 1, Eigen::ColMajor>;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
public:
    //--------------------------------------------------------------------//
    // Constructor
    MultiFluid(
            std::string const& prefix,
            worldcomm_ptr_t const& wc = Environment::worldCommPtr(),
            std::string const& subPrefix = "",
            ModelBaseRepository const& modelRep = ModelBaseRepository() );
    MultiFluid( self_type const& M ) = default;

    static self_ptrtype New(
            std::string const& prefix,
            worldcomm_ptr_t const& wc = Environment::worldCommPtr(),
            std::string const& subPrefix = "",
            ModelBaseRepository const& modelRep = ModelBaseRepository() );

    void loadParametersFromOptionsVm();

    static std::string expandStringFromSpec( std::string const& expr );

    //--------------------------------------------------------------------//
    // Initialization
    void init();

    //--------------------------------------------------------------------//
    std::string const& prefix() const { return M_prefix; }
    std::string const& fluidPrefix() const { return super_type::prefix(); }
    std::string globalLevelsetPrefix() const { return prefixvm( this->prefix(), "levelset"); }

    static std::string levelsetName( uint16_type n ) { return (boost::format( "levelset%1%" ) %(n+1)).str(); }

    std::shared_ptr<std::ostringstream> getInfo() const override;

    //--------------------------------------------------------------------//
    std::shared_ptr<self_type> shared_from_this() { return std::static_pointer_cast<self_type>(super_type::shared_from_this()); }
    std::shared_ptr<self_type const> shared_from_this() const { return std::static_pointer_cast<self_type const>(super_type::shared_from_this()); }

    //--------------------------------------------------------------------//
    // Function spaces
    space_levelset_ptrtype const& functionSpaceLevelset() const { return M_levelsetSpaceManager->functionSpaceScalar(); }
    space_levelset_vectorial_ptrtype const& functionSpaceLevelsetVectorial() const { return M_levelsetSpaceManager->functionSpaceVectorial(); }
    space_inextensibilitylm_ptrtype const& functionSpaceInextensibilityLM() const;
    //--------------------------------------------------------------------//
    // Mesh
    //mesh_ptrtype const& mesh() const { return M_mesh; }
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"MultiFluidMesh.path"); }
    //--------------------------------------------------------------------//
    // Operator Lagrange P1
    bool useLagrangeP1iso() const { return M_useLagrangeP1iso; }
    op_lagrangeP1_ptrtype opLagrangeP1() const { return M_opLagrangeP1iso; }
    //--------------------------------------------------------------------//
    // Models
    fluid_ptrtype fluidModel() { return this->shared_from_this(); }
    std::map<std::string, levelset_ptrtype> levelsetModels() const { return M_levelsets; }
    levelset_ptrtype const& levelsetModel(std::string const& name) const { return M_levelsets.at(name); }
    uint16_type nLevelsets() const { return M_levelsets.size(); }

    // Global levelset
    auto globalLevelsetExpr() const {
        std::vector< element_levelset_ptrtype > levelsets;
        std::transform( M_levelsets.begin(), M_levelsets.end(), std::back_inserter(levelsets),
                [](std::pair<std::string, levelset_ptrtype> const& l) { return l.second->phi(); }
                );
        return Feel::FeelModels::globalLevelsetExpr( levelsets );
    }
    element_levelset_ptrtype const& globalLevelsetElt() const;

    // Fluid density-viscosity model
    material_properties_ptrtype const& fluidMaterialProperties() const { return M_fluidMaterialProperties; }
    material_properties_ptrtype const& levelsetMaterialProperties( std::string const& name ) const { return M_levelsetsMaterialProperties.at(name); }

    //--------------------------------------------------------------------//
    // Algebraic data
    int nBlockMatrixGraph() const override;
    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    size_type nLocalDof() const override;

    //--------------------------------------------------------------------//
    double globalLevelsetThicknessInterface() const { return M_globalLevelsetThicknessInterface; }
    //--------------------------------------------------------------------//
    bool hasInextensibility( std::string const& name ) const { return M_hasInextensibility.at(name); }
    std::string const& inextensibilityMethod( std::string const& name ) const { return M_inextensibilityMethods.at(name); }
    bool hasInextensibilityLM() const;
    void updateInextensibilityLM();
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
    //std::shared_ptr<TSBase> timeStepBase() const { return M_fluid->timeStepBase(); }
    std::shared_ptr<TSBase> fluidTimeStepBase() const { return this->timeStepBase(); }
    std::shared_ptr<TSBase> levelsetTimeStepBase( std::string const& name) const { return this->levelsetModel(name)->timeStepBase(); }
    void updateTime( double time );
    void updateTimeStep();

    //--------------------------------------------------------------------//
    // Export
    void exportResults( double time ) { this->exportResultsImpl( time ); }
    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportMeasures( double time );
    // Measures
    force_type computeLevelsetForce( std::string const& name ) const;

protected:
    //--------------------------------------------------------------------//
    // Initialization
    mesh_ptrtype createMesh();
    void createLevelsets();
    void initPostProcess();

    size_type initStartBlockIndexFieldsInMatrix() override;
    int initBlockVector() override;

    //--------------------------------------------------------------------//
    void updateFluidDensityViscosity();
    void updateInterfaceForces();
    void solveFluid();
    void advectLevelsets();
    virtual void solveImpl();

    void setRebuildMatrixVector( bool b = true ) { M_doRebuildMatrixVector = b; }
    bool rebuildMatrixVector() const { return M_doRebuildMatrixVector; }
    // Linear solve
    void updateLinearPDEAdditional( DataUpdateLinear & data ) const;
    // Non-linear solve
    void updateJacobianAdditional( DataUpdateJacobian & data ) const;
    void updateResidualAdditional( DataUpdateResidual & data ) const;
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
    levelset_space_manager_ptrtype M_levelsetSpaceManager;
    levelset_tool_manager_ptrtype M_levelsetToolManager;
    std::map<std::string, levelset_ptrtype> M_levelsets;
    mutable element_levelset_ptrtype M_globalLevelsetElt;
    mutable bool M_doUpdateGlobalLevelset;
    exporter_ptrtype M_globalLevelsetExporter;

    //--------------------------------------------------------------------//
    // Solve
    bool M_doRebuildMatrixVector;
    bool M_usePicardIterations;

    //--------------------------------------------------------------------//
    // Parameters
    material_properties_ptrtype M_fluidMaterialProperties;
    std::map<std::string, material_properties_ptrtype> M_levelsetsMaterialProperties;
    std::map<std::string, std::map<std::string, interfaceforces_model_ptrtype>> M_levelsetInterfaceForcesModels;
    std::map<std::string, interfaceforces_model_ptrtype> M_additionalInterfaceForcesModel;
    //--------------------------------------------------------------------//
    // Global levelset parameters
    double M_globalLevelsetThicknessInterface;
    //--------------------------------------------------------------------//
    // Forces
    bool M_hasInterfaceForcesModel;

    element_levelset_vectorial_ptrtype M_interfaceForces; 

    //--------------------------------------------------------------------//
    // Inextensibility
    bool M_enableInextensibility;
    std::map<std::string, bool> M_hasInextensibility;
    bool M_hasInextensibilityLM;
    std::map<std::string, std::string> M_inextensibilityMethods;
    std::vector<element_levelset_ptrtype> M_inextensibleLevelsets;

    mutable range_elements_type M_rangeInextensibilityLM;
    mutable space_inextensibilitylm_ptrtype M_spaceInextensibilityLM;
    mutable bool M_doUpdateInextensibilityLM;
    // Penalty method gamma
    std::map<std::string, double> M_inextensibilityGamma;
    //--------------------------------------------------------------------//
    // Reinitialization
    std::map<std::string, int> M_levelsetReinitEvery;

    //--------------------------------------------------------------------//
    // Post-process
    // levelset forces
    std::vector< std::string > M_postProcessMeasuresLevelsetForces;
};
        

} // namespace FeelModels
} // namespace Feel

#endif
