#ifndef _MULTIFLUID_HPP
#define _MULTIFLUID_HPP 1

#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelmodels/levelset/levelset.hpp>
#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>
#include <feel/feelmodels/levelset/globallevelsetexpr.hpp>

namespace Feel {
namespace FeelModels {

template< typename FluidType, typename LevelSetType>
class MultiFluid : 
    public ModelNumerical,
    public std::enable_shared_from_this< MultiFluid<FluidType, LevelSetType> >
{
public:
    // Typedefs
    //--------------------------------------------------------------------//
    // Class
    typedef ModelNumerical super_type;
    using size_type = typename super_type::size_type;
    typedef MultiFluid< FluidType, LevelSetType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef FluidType fluid_model_type;
    typedef std::shared_ptr<fluid_model_type> fluid_model_ptrtype;

    typedef LevelSetType levelset_model_type;
    typedef std::shared_ptr<levelset_model_type> levelset_model_ptrtype;

    //--------------------------------------------------------------------//
    // Mesh
    typedef typename fluid_model_type::convex_type convex_type;
    static inline const uint16_type nDim = convex_type::nDim;
    static inline const uint16_type nOrderGeo = convex_type::nOrder;
    static inline const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    //--------------------------------------------------------------------//
    // Range types
    typedef typename MeshTraits<mesh_type>::element_reference_wrapper_const_iterator element_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
    typedef Range<mesh_type,MESH_ELEMENTS> range_elements_type;

    typedef typename MeshTraits<mesh_type>::face_reference_wrapper_const_iterator face_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_type faces_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype faces_reference_wrapper_ptrtype;
    typedef Range<mesh_type,MESH_FACES> range_faces_type;

    //--------------------------------------------------------------------//
    // Function spaces
    typedef typename levelset_model_type::space_levelset_ptrtype space_levelset_ptrtype;
    typedef typename levelset_model_type::space_vectorial_ptrtype space_levelset_vectorial_ptrtype;
    typedef typename levelset_model_type::space_markers_ptrtype space_levelset_markers_ptrtype;
    typedef typename levelset_model_type::space_advection_velocity_type space_levelset_advection_velocity_type;
    typedef typename levelset_model_type::space_advection_velocity_ptrtype space_levelset_advection_velocity_ptrtype;
    typedef typename fluid_model_type::component_space_velocity_type component_space_fluid_velocity_type;

    typedef typename fluid_model_type::space_velocity_type space_fluid_velocity_type;
    typedef typename fluid_model_type::space_velocity_ptrtype space_fluid_velocity_ptrtype;

    typedef typename levelset_model_type::element_levelset_type element_levelset_type;
    typedef typename levelset_model_type::element_levelset_ptrtype element_levelset_ptrtype; 
    typedef typename levelset_model_type::element_vectorial_type element_levelset_vectorial_type;
    typedef typename levelset_model_type::element_vectorial_ptrtype element_levelset_vectorial_ptrtype; 

    // Levelset function space manager
    typedef typename levelset_model_type::levelset_space_manager_type levelset_space_manager_type;
    typedef typename levelset_model_type::levelset_space_manager_ptrtype levelset_space_manager_ptrtype;
    // Levelset tool manager
    typedef typename levelset_model_type::levelset_tool_manager_type levelset_tool_manager_type;
    typedef typename levelset_model_type::levelset_tool_manager_ptrtype levelset_tool_manager_ptrtype;

    //--------------------------------------------------------------------//
    // Lagrange P1 iso-Pn
    typedef OperatorLagrangeP1<component_space_fluid_velocity_type> op_lagrangeP1_type;
    typedef std::shared_ptr<op_lagrangeP1_type> op_lagrangeP1_ptrtype;

    //--------------------------------------------------------------------//
    // Density/viscosity
    typedef typename fluid_model_type::material_properties_type material_properties_type;
    typedef typename fluid_model_type::material_properties_ptrtype material_properties_ptrtype;
    //--------------------------------------------------------------------//
    // Interface forces model
    typedef InterfaceForcesModel<levelset_model_type, fluid_model_type> interfaceforces_model_type;
    typedef std::shared_ptr<interfaceforces_model_type> interfaceforces_model_ptrtype;
    typedef Singleton<Feel::Factory<interfaceforces_model_type, std::string>> interfaceforces_factory_type;

    //--------------------------------------------------------------------//
    // Inextensibility
    typedef typename fluid_model_type::basis_fluid_p_type basis_fluid_p_type;
    typedef FunctionSpace< mesh_type, bases<basis_fluid_p_type> > space_inextensibilitylm_type;
    typedef std::shared_ptr<space_inextensibilitylm_type> space_inextensibilitylm_ptrtype;
    //--------------------------------------------------------------------//
    // Exporter
    typedef Exporter<mesh_type, nOrderGeo> exporter_type;
    typedef std::shared_ptr<exporter_type> exporter_ptrtype;

    //--------------------------------------------------------------------//
    // Exporter
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr<model_algebraic_factory_type> model_algebraic_factory_ptrtype;

    //--------------------------------------------------------------------//
    // Measures
    using force_type = Eigen::Matrix<typename fluid_model_type::value_type, nDim, 1, Eigen::ColMajor>;

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
    void init( bool buildModelAlgebraicFactory = true );
    void initAlgebraicFactory();

    //--------------------------------------------------------------------//
    std::string globalLevelsetPrefix() const { return prefixvm( this->prefix(), "levelset"); }

    static std::string levelsetName( uint16_type n ) { return (boost::format( "levelset%1%" ) %(n+1)).str(); }

    std::shared_ptr<std::ostringstream> getInfo() const override;

    //--------------------------------------------------------------------//
    // Function spaces
    space_levelset_ptrtype const& functionSpaceLevelset() const { return M_levelsetSpaceManager->functionSpaceScalar(); }
    space_levelset_vectorial_ptrtype const& functionSpaceLevelsetVectorial() const { return M_levelsetSpaceManager->functionSpaceVectorial(); }
    space_inextensibilitylm_ptrtype const& functionSpaceInextensibilityLM() const;
    //--------------------------------------------------------------------//
    // Mesh
    mesh_ptrtype const& mesh() const { return M_mesh; }
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"MultiFluidMesh.path"); }
    //--------------------------------------------------------------------//
    // Operator Lagrange P1
    bool useLagrangeP1iso() const { return M_useLagrangeP1iso; }
    op_lagrangeP1_ptrtype opLagrangeP1() const { return M_opLagrangeP1iso; }
    //--------------------------------------------------------------------//
    // Models
    fluid_model_ptrtype const& fluidModel() const { return M_fluidModel; }
    fluid_model_ptrtype fluidModel() { return M_fluidModel; }

    std::map<std::string, levelset_model_ptrtype> const& levelsetModels() const { return M_levelsets; }
    std::map<std::string, levelset_model_ptrtype> & levelsetModels() { return M_levelsets; }
    levelset_model_ptrtype const& levelsetModel(std::string const& name) const { return M_levelsets.at(name); }
    uint16_type nLevelsets() const { return M_levelsets.size(); }

    // Global levelset
    auto globalLevelsetExpr() const {
        std::vector< element_levelset_ptrtype > levelsets;
        std::transform( M_levelsets.begin(), M_levelsets.end(), std::back_inserter(levelsets),
                [](std::pair<std::string, levelset_model_ptrtype> const& l) { return l.second->phi(); }
                );
        return Feel::FeelModels::globalLevelsetExpr( levelsets );
    }
    element_levelset_ptrtype const& globalLevelsetElt( bool up = true ) const;
    void updateGlobalLevelsetElt( element_levelset_ptrtype & globalLevelsetElt, bool & doUpdateGlobalLevelset ) const;

    // Fluid density-viscosity model
    material_properties_ptrtype const& fluidMaterialProperties() const { return M_fluidMaterialProperties; }
    material_properties_ptrtype const& levelsetMaterialProperties( std::string const& name ) const { return M_levelsetsMaterialProperties.at(name); }

    //--------------------------------------------------------------------//
    // Accessors interfaces
    //decltype(auto) functionSpaceVelocityPressure() const { return this->fluidModel()->functionSpace(); }
    decltype(auto) functionSpaceVelocity() const { return this->fluidModel()->functionSpaceVelocity(); }
    decltype(auto) functionSpacePressure() const { return this->fluidModel()->functionSpacePressure(); }
    //decltype(auto) fieldVelocityPressurePtr() const { return this->fluidModel()->fieldVelocityPressurePtr(); }
    //decltype(auto) fieldVelocityPressure() const { return this->fluidModel()->fieldVelocityPressure(); }
    decltype(auto) fieldVelocity() const { return this->fluidModel()->fieldVelocity(); }
    decltype(auto) fieldPressure() const { return this->fluidModel()->fieldPressure(); }

    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    struct FieldTag
    {
        static auto globalLevelsetElt( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
    };
    auto modelFields( std::string const& prefix = "" ) const
        {
            auto mfieldsFluid = M_fluidModel->modelFields( prefixvm( prefix, M_fluidModel->keyword() ) );

            using mfields_levelset_type = std::decay_t<decltype(M_levelsets.begin()->second->modelFields( "" ) )>;
            mfields_levelset_type mfieldsLevelsets;
            // for ( auto const& [lsName,lsObject] : M_levelsets )
            //     mfieldsLevelsets = Feel::FeelModels::modelFields( mfieldsLevelsets, lsObject->modelFields(  prefixvm( prefix,  lsObject->keyword() ) ) );
            for ( auto it = M_levelsets.begin() ; it != M_levelsets.end() ; ++it )
                mfieldsLevelsets = Feel::FeelModels::modelFields( mfieldsLevelsets, it->second->modelFields(  prefixvm( prefix,  it->second->keyword() ) ) );

            return Feel::FeelModels::modelFields( mfieldsFluid, mfieldsLevelsets,
                                                  modelField<FieldCtx::ID>( FieldTag::globalLevelsetElt(this), prefixvm( prefix, "global-levelset"), "phi", this->globalLevelsetElt(false), "phi", this->keyword(),
                                                                            std::bind( &self_type::updateGlobalLevelsetElt, this, std::placeholders::_1, std::ref(M_doUpdateGlobalLevelset) ) )
                                                  );
        }

    //___________________________________________________________________________________//
    // symbols expression
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            auto seFluid = this->fluidModel()->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            //auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr();
            return Feel::vf::symbolsExpr( seFluid,seParam,seFields );
        }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

#if 0
    //--------------------------------------------------------------------//
    // Symbols expr
    auto symbolsExpr() const { 
        return this->symbolsExpr( M_fluidModel->fieldVelocity(), M_fluidModel->fieldPressure() ); 
        // TODO add levelsets symbols
    }

    template <typename FieldVelocityType, typename FieldPressureType>
    auto symbolsExpr( FieldVelocityType const& u, FieldPressureType const& p ) const
    {
        auto seFluid = this->fluidModel()->symbolsExprToolbox( u,p );
        auto seParam = this->symbolsExprParameter();
        //auto symbolExprMaterial = Feel::vf::symbolsExpr( M_fluidModel->symbolsExprMaterial( Feel::vf::symbolsExpr( symbolExprField, symbolExprFit ) ) );
        return Feel::vf::symbolsExpr( seFluid, seParam );
    }
#endif
    //--------------------------------------------------------------------//
    // Algebraic data
    int nBlockMatrixGraph() const;
    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    size_type nLocalDof() const;
    void buildBlockVector();

    backend_ptrtype const& backend() const { return M_backend; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    BlocksBaseVector<double> & blockVectorSolution() { return M_blockVectorSolution; }

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
    virtual void solveExplicitCoupling();
    virtual void solveImplicitCoupling();
    void solvePicard();

    void updateLinearPDE( DataUpdateLinear & data ) const override;
    void updateLinearPDEInterfaceForces( DataUpdateLinear & data ) const;
    void updateLinearPDEInextensibility( DataUpdateLinear & data ) const;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;

    //void updateNewtonInitialGuess( vector_ptrtype& U ) const;
    void updateJacobian( DataUpdateJacobian & data ) const override;
    void updateJacobianInterfaceForces( DataUpdateJacobian & data ) const;
    void updateJacobianInextensibility( DataUpdateJacobian & data ) const;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
    void updateResidual( DataUpdateResidual & data ) const override;
    void updateResidualInterfaceForces( DataUpdateResidual & data ) const;
    void updateResidualInextensibility( DataUpdateResidual & data ) const;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;
    
    void updateParameterValues();

    //--------------------------------------------------------------------//
    // Time step
    std::shared_ptr<TSBase> timeStepBase() const { return M_fluidModel->timeStepBase(); }
    std::shared_ptr<TSBase> fluidTimeStepBase() const { return this->timeStepBase(); }
    std::shared_ptr<TSBase> levelsetTimeStepBase( std::string const& name) const { return this->levelsetModel(name)->timeStepBase(); }
    //void updateTime( double time );
    void startTimeStep();
    void updateTimeStep();

    //--------------------------------------------------------------------//
    // Export
    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );
    void exportMeasures( double time );
    // Measures
    force_type computeLevelsetForce( std::string const& name ) const;

protected:
    //--------------------------------------------------------------------//
    // Initialization
    void initMesh();
    void initLevelsets();
    void initPostProcess() override;

    virtual int initBlockVector();
    bool useImplicitCoupling() const;

    //--------------------------------------------------------------------//
    void updateFluidDensityViscosity();
    void updateInterfaceForces();
    void solveFluid();
    void advectLevelsets();

    void setRebuildMatrixVector( bool b = true ) { M_doRebuildMatrixVector = b; }
    bool rebuildMatrixVector() const { return M_doRebuildMatrixVector; }

    //--------------------------------------------------------------------//
    uint16_type M_nFluids;

private:
    std::string M_prefix;
    //--------------------------------------------------------------------//
    // Lagrange P1 iso-Pn
    bool M_useLagrangeP1iso;
    op_lagrangeP1_ptrtype M_opLagrangeP1iso;
    //--------------------------------------------------------------------//
    mesh_ptrtype M_mesh;
    fluid_model_ptrtype M_fluidModel;
    levelset_space_manager_ptrtype M_levelsetSpaceManager;
    levelset_tool_manager_ptrtype M_levelsetToolManager;
    std::map<std::string, levelset_model_ptrtype> M_levelsets;
    mutable element_levelset_ptrtype M_globalLevelsetElt;
    mutable bool M_doUpdateGlobalLevelset;
    exporter_ptrtype M_globalLevelsetExporter;

    //--------------------------------------------------------------------//
    // Solve
    bool M_doRebuildMatrixVector;
    bool M_usePicardIterations;

    //--------------------------------------------------------------------//
    // Algebraic data
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;

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
    // Redistanciation
    std::map<std::string, int> M_levelsetRedistEvery;

    //--------------------------------------------------------------------//
    // Post-process
    exporter_ptrtype M_exporter;
    // fluid exports
    std::set<std::string> M_postProcessFieldsExportedFluid;
    // levelset forces
    std::vector< std::string > M_postProcessMeasuresLevelsetForces;
};
        

} // namespace FeelModels
} // namespace Feel

#endif
