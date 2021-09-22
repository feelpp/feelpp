#ifndef _MULTIFLUID_HPP
#define _MULTIFLUID_HPP 1

#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelmodels/levelset/levelset.hpp>
#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>
#include <feel/feelmodels/levelset/globallevelsetexpr.hpp>

#include <feel/feelmodels/modelcore/modelfields.hpp>

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
    typedef MultiFluid< FluidType, LevelSetType > self_type;
    typedef std::shared_ptr< self_type > self_ptrtype;

    typedef FluidType fluid_model_type;
    typedef std::shared_ptr<fluid_model_type> fluid_model_ptrtype;

    typedef LevelSetType levelset_model_type;
    typedef std::shared_ptr<levelset_model_type> levelset_model_ptrtype;

    using size_type = typename super_type::size_type;

    //--------------------------------------------------------------------//
    // Mesh
    typedef typename fluid_model_type::mesh_type mesh_fluid_type;
    typedef typename levelset_model_type::mesh_type mesh_levelset_type;
    typedef mesh_fluid_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    static const uint16_type nDim = mesh_type::nDim;
    static const uint16_type nOrderGeo = mesh_type::nOrder;
    static const uint16_type nRealDim = mesh_type::nRealDim;

    //--------------------------------------------------------------------//
    // Materials properties
    typedef MaterialsProperties<mesh_type::nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;

    ////--------------------------------------------------------------------//
    //// Range types
    //typedef typename MeshTraits<mesh_type>::element_reference_wrapper_const_iterator element_reference_wrapper_const_iterator;
    //typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
    //typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
    //typedef elements_reference_wrapper_t<mesh_type> range_elements_type;

    //typedef typename MeshTraits<mesh_type>::face_reference_wrapper_const_iterator face_reference_wrapper_const_iterator;
    //typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_type faces_reference_wrapper_type;
    //typedef typename MeshTraits<mesh_type>::faces_reference_wrapper_ptrtype faces_reference_wrapper_ptrtype;
    //typedef faces_reference_wrapper_t<mesh_type> range_faces_type;

    //--------------------------------------------------------------------//
    // Function spaces
    typedef typename levelset_model_type::space_levelset_ptrtype space_levelset_ptrtype;
    typedef typename levelset_model_type::space_vectorial_ptrtype space_levelset_vectorial_ptrtype;
    typedef typename levelset_model_type::space_markers_ptrtype space_levelset_markers_ptrtype;
    typedef typename fluid_model_type::component_space_velocity_type component_space_fluid_velocity_type;

    typedef typename fluid_model_type::space_velocity_type space_fluid_velocity_type;
    typedef typename fluid_model_type::space_velocity_ptrtype space_fluid_velocity_ptrtype;

    typedef typename levelset_model_type::element_scalar_type element_levelset_scalar_type;
    typedef typename levelset_model_type::element_scalar_ptrtype element_levelset_scalar_ptrtype; 
    typedef typename levelset_model_type::element_vectorial_type element_levelset_vectorial_type;
    typedef typename levelset_model_type::element_vectorial_ptrtype element_levelset_vectorial_ptrtype; 

    // Levelset function space manager
    typedef typename levelset_model_type::levelset_space_manager_type levelset_space_manager_type;
    typedef typename levelset_model_type::levelset_space_manager_ptrtype levelset_space_manager_ptrtype;
    // Levelset tool manager
    typedef typename levelset_model_type::levelset_tool_manager_type levelset_tool_manager_type;
    typedef typename levelset_model_type::levelset_tool_manager_ptrtype levelset_tool_manager_ptrtype;

    ////--------------------------------------------------------------------//
    //// Lagrange P1 iso-Pn
    //typedef OperatorLagrangeP1<component_space_fluid_velocity_type> op_lagrangeP1_type;
    //typedef std::shared_ptr<op_lagrangeP1_type> op_lagrangeP1_ptrtype;

    ////--------------------------------------------------------------------//
    //// Density/viscosity
    //typedef typename fluid_model_type::materialsproperties_type materials_properties_type;
    //typedef typename fluid_model_type::materialsproperties_ptrtype materials_properties_ptrtype;
    ////--------------------------------------------------------------------//
    //// Interface forces model
    //typedef InterfaceForcesModel<levelset_model_type, fluid_model_type> interfaceforces_model_type;
    //typedef std::shared_ptr<interfaceforces_model_type> interfaceforces_model_ptrtype;
    //typedef Singleton<Feel::Factory<interfaceforces_model_type, std::string>> interfaceforces_factory_type;

    ////--------------------------------------------------------------------//
    //// Inextensibility
    //typedef typename fluid_model_type::basis_fluid_p_type basis_fluid_p_type;
    //typedef FunctionSpace< mesh_type, bases<basis_fluid_p_type> > space_inextensibilitylm_type;
    //typedef std::shared_ptr<space_inextensibilitylm_type> space_inextensibilitylm_ptrtype;

    //--------------------------------------------------------------------//
    // Cached fields
    using cached_levelset_scalar_field_type = CachedModelField<element_levelset_scalar_type, InplaceUpdatePolicy>;
    using cached_levelset_vectorial_field_type = CachedModelField<element_levelset_vectorial_type, InplaceUpdatePolicy>;
    //--------------------------------------------------------------------//
    // Exporter
    typedef Exporter<mesh_type, nOrderGeo> exporter_type;
    typedef std::shared_ptr<exporter_type> exporter_ptrtype;

    //--------------------------------------------------------------------//
    // Algebraic factory
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr<model_algebraic_factory_type> model_algebraic_factory_ptrtype;

    //--------------------------------------------------------------------//
    // Measures
    using force_type = Eigen::Matrix<typename fluid_model_type::value_type, nDim, 1, Eigen::ColMajor>;

public:
    //--------------------------------------------------------------------//
    // Constructor
    MultiFluid(
            std::string const& prefix,
            std::string const& keyword = "multifluid",
            worldcomm_ptr_t const& wc = Environment::worldCommPtr(),
            std::string const& subPrefix = "",
            ModelBaseRepository const& modelRep = ModelBaseRepository() );
    MultiFluid( self_type const& M ) = default;

    static self_ptrtype New(
            std::string const& prefix,
            worldcomm_ptr_t const& wc = Environment::worldCommPtr(),
            std::string const& subPrefix = "",
            ModelBaseRepository const& modelRep = ModelBaseRepository() );

    //--------------------------------------------------------------------//
    // Initialization
    void init( bool buildModelAlgebraicFactory = true );
    void initAlgebraicFactory();

    void loadParametersFromOptionsVm();

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
    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }
    //std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"MultiFluidMesh.path"); }
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
        std::vector< element_levelset_scalar_ptrtype > levelsets;
        std::transform( this->levelsetModels().begin(), this->levelsetModels().end(), std::back_inserter(levelsets),
                []( std::pair<std::string, levelset_model_ptrtype> const& l ) { return l.second->phi(); }
                );
        return Feel::FeelModels::globalLevelsetExpr( levelsets );
    }
    element_levelset_ptrtype globalLevelsetElt() const { return M_globalLevelset.fieldPtr(); }

    // Physical parameters
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }
    //// Fluid density-viscosity model
    //material_properties_ptrtype const& fluidMaterialProperties() const { return M_fluidMaterialProperties; }
    //material_properties_ptrtype const& levelsetMaterialProperties( std::string const& name ) const { return M_levelsetsMaterialProperties.at(name); }

    //--------------------------------------------------------------------//
    // Time stepping
    std::shared_ptr<TSBase> timeStepBase() { return this->fluidModel()->timeStepBase(); }
    std::shared_ptr<TSBase> timeStepBase() const { return this->fluidModel()->timeStepBase(); }
    void startTimeStep();
    void updateTimeStep();

    ////--------------------------------------------------------------------//
    //// Accessors interfaces
    //decltype(auto) functionSpaceVelocity() const { return this->fluidModel()->functionSpaceVelocity(); }
    //decltype(auto) functionSpacePressure() const { return this->fluidModel()->functionSpacePressure(); }
    //decltype(auto) fieldVelocity() const { return this->fluidModel()->fieldVelocity(); }
    //decltype(auto) fieldPressure() const { return this->fluidModel()->fieldPressure(); }

    //--------------------------------------------------------------------//
    // Fields
    struct FieldTag
    {
        static auto globalLevelsetElt( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
    };
    auto modelFieldsLevelsets( std::string const& prefix = "" ) const
    {
        using mfields_levelset_type = std::decay_t<decltype(this->levelsetModels().begin()->second->modelFields( "" ) )>;
        mfields_levelset_type mfieldsLevelsets;
        //for ( auto it = M_levelsets.begin() ; it != M_levelsets.end() ; ++it )
            //mfieldsLevelsets = Feel::FeelModels::modelFields( 
                    //mfieldsLevelsets, 
                    //it->second->modelFields(  prefixvm( prefix,  it->second->keyword() ) ) 
                    //);
        for ( auto const& [lsName,lsModel] : this->levelsetModels() )
            mfieldsLevelsets = Feel::FeelModels::modelFields( 
                    mfieldsLevelsets, 
                    lsModel->modelFields( lsModel->keyword() )
                    );

        //return Feel::FeelModels::modelFields( 
                //mfieldsLevelsets,
                //modelField<FieldCtx::ID>( FieldTag::globalLevelsetElt(this), prefix, "global-levelset", M_globalLevelset, "phi", this->keyword() )
                //);
        return mfieldsLevelsets;
    }
    auto modelFieldsLevelsets( vector_ptrtype sol, std::vector<size_type> rowStartInVectorLevelsets, std::string const& prefix = "" ) const
    {
        using mfields_levelset_type = std::decay_t<decltype(this->levelsetModels().begin()->second->modelFields( "" ) )>;
        mfields_levelset_type mfieldsLevelsets;
        for ( size_type i = 0; i < this->levelsetModels()->size(); ++i )
            mfieldsLevelsets = Feel::FeelModels::modelFields(
                    mfieldsLevelsets,
                    this->levelsetModel(i)->modelFields( sol, rowStartInVectorLevelsets[i], this->levelsetModel(i)->keyword() )
                    );

        //return Feel::FeelModels::modelFields( 
                //mfieldsLevelsets,
                //modelField<FieldCtx::ID>( FieldTag::globalLevelsetElt(this), prefix, "global-levelset", M_globalLevelset, "phi", this->keyword() )
                //);
        return mfieldsLevelsets;
    }
    auto modelFieldsLevelsets( std::vector<vector_ptrtype> sols, std::vector<size_type> rowStartInVectorLevelsets, std::string const& prefix = "" ) const
    {
        using mfields_levelset_type = std::decay_t<decltype(this->levelsetModels().begin()->second->modelFields( "" ) )>;
        mfields_levelset_type mfieldsLevelsets;
        for ( size_type i = 0; i < this->levelsetModels()->size(); ++i )
            mfieldsLevelsets = Feel::FeelModels::modelFields(
                    mfieldsLevelsets,
                    this->levelsetModel(i)->modelFields( sol[i], rowStartInVectorLevelsets[i], this->levelsetModel(i)->keyword() )
                    );

        //return Feel::FeelModels::modelFields( 
                //mfieldsLevelsets,
                //modelField<FieldCtx::ID>( FieldTag::globalLevelsetElt(this), prefix, "global-levelset", M_globalLevelset, "phi", this->keyword() )
                //);
        return mfieldsLevelsets;
    }

    auto modelFields( std::string const& prefix = "" ) const
    {
        return Feel::FeelModels::modelFields( 
                this->fluidModel()->modelFields( this->fluidModel()->keyword() ),
                this->modelFieldsLevelsets( prefix )
                );
    }
    auto modelFields( vector_ptrtype sol, size_type rowStartInVectorFluid, std::vector<size_type> rowStartInVectorLevelsets, std::string const& prefix = "" ) const
    {
        return Feel::FeelModels::modelFields( 
                this->fluidModel()->modelFields( sol, rowStartInVectorFluid, this->fluidModel()->keyword() ),
                this->modelFieldsLevelsets( sol, rowStartInVectorLevelsets, prefix )
                );
    }
    auto modelFields( vector_ptrtype solFluid, size_type rowStartInVectorFluid, std::vector<vector_ptrtype> solLevelsets, std::vector<size_type> rowStartInVectorLevelsets, std::string const& prefix = "" ) const
    {
        return Feel::FeelModels::modelFields( 
                this->fluidModel()->modelFields( solFluid, rowStartInVectorFluid, this->fluidModel()->keyword() ),
                this->modelFieldsLevelsets( solLevelsets, rowStartInVectorLevelsets, prefix )
                );
    }
    //auto modelFields( 
            //std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorDataFluid,
            //std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorDataLevelsets,
            //std::string const& prefix = "" ) const
    //{
        ////TODO
    //}

    auto trialSelectorModelFieldsLevelsets( std::vector<size_type> startBlockSpaceIndexLevelsets ) const
    {
        using tsmfields_levelset_type = std::decay_t<decltype(this->levelsetModel()->trialSelectorModelFields( "" ) )>;
        tsmfields_levelset_type tsmfieldsLevelsets;
        for ( size_type i = 0; i < this->levelsetModels()->size(); ++i )
            tsmfieldsLevelsets = Feel::FeelModels::selectorModelFields(
                    tsmfieldsLevelsets,
                    this->levelsetModel(i)->trialSelectorModelFields( startBlockSpaceIndexLevelsets[i] )
                    );

        return tsmfieldsLevelsets;
    }
    auto trialSelectorModelFields( size_type startBlockSpaceIndexFluid, std::vector<size_type> startBlockSpaceIndexLevelsets ) const
    {
        return Feel::FeelModels::selectorModelFields( 
                this->fluidModel()->trialSelectorModelFields( startBlockSpaceIndexFluid ),
                this->trialSelectorModelFieldsLevelsets( startBlockSpaceIndexLevelsets )
                );
    }

    //--------------------------------------------------------------------//
    // Symbols expression
    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
    {
        auto seFluid = this->fluidModel()->symbolsExprToolbox( mfields );
        // TODO: add levelset symbols
        //auto seLevelsets = this->levelsetModels()->symbolsExprToolbox( mfields );
        auto seParam = this->symbolsExprParameter();
        auto seMeshes = this->template symbolsExprMeshes<mesh_type>();
        auto seMat = this->materialsProperties()->symbolsExpr();
        auto seFields = mfields.symbolsExpr();
        auto sePhysics = this->symbolsExprPhysics( this->physics() );
        return Feel::vf::symbolsExpr( seFluid,seParam,seFields );
    }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

    //--------------------------------------------------------------------//
    // Model context
    auto modelContext( std::string const& prefix = "" ) const
    {
        auto mfields = this->modelFields( prefix );
        auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
        return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
    }
    auto modelContext( vector_ptrtype sol, size_type startBlockSpaceIndexFluid, std::vector<size_type> startBlockSpaceIndexLevelsets, std::string const& prefix = "" ) const
    {
        return this->modelContext( sol, startBlockSpaceIndexFluid, sol, startBlockSpaceIndexLevelsets, prefix );
    }
    auto modelContext( vector_ptrtype solFluid, size_type startBlockSpaceIndexFluid, std::vector<vector_ptrtype> solLevelsets, std::vector<size_type> startBlockSpaceIndexLevelsets, std::string const& prefix = "" ) const
    {
        auto mfields = this->modelFields( solFluid, startBlockSpaceIndexFluid, solLevelsets, startBlockSpaceIndexLevelsets, prefix );
        auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
        auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( startBlockSpaceIndexFluid, startBlockSpaceIndexLevelsets ) );
        return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ), std::move( tse ) );
    }
    auto modelContextNoTrialSymbolsExpr( vector_ptrtype sol, size_type startBlockSpaceIndexFluid, std::vector<size_type> startBlockSpaceIndexLevelsets, std::string const& prefix = "" ) const
    {
        return this->modelContextNoTrialSymbolsExpr( sol, startBlockSpaceIndexFluid, sol, startBlockSpaceIndexLevelsets, prefix );
    }
    auto modelContextNoTrialSymbolsExpr( vector_ptrtype solHeat, size_type startBlockSpaceIndexHeat, vector_ptrtype solFluid, size_type startBlockSpaceIndexFluid, std::string const& prefix = "" ) const
    {
        // auto mfields = this->modelFields( solHeat, startBlockSpaceIndexHeat, solFluid, startBlockSpaceIndexFluid, prefix );
        // auto se = this->symbolsExpr( mfields );
        // return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        return this->modelContextNoTrialSymbolsExpr( { { "solution", std::make_tuple( solHeat, startBlockSpaceIndexHeat ) } },
                { { "solution", std::make_tuple( solFluid, startBlockSpaceIndexFluid ) } },
                prefix );
    }
    auto modelContextNoTrialSymbolsExpr( std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorDataHeat,
                                         std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorDataFluid,
                                         std::string const& prefix = "" ) const
    {
        auto mfields = this->modelFields( vectorDataHeat, vectorDataFluid, prefix );
        auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
        return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
    }
    auto modelContext( vector_ptrtype sol, heat_model_ptrtype const& heatModel, fluid_model_ptrtype const& fluidModel, std::string const& prefix = "" ) const
    {
        return this->modelContext( sol, heatModel->startBlockSpaceIndexVector(), fluidModel->startBlockSpaceIndexVector(), prefix );
    }

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
    // Assembly and solve
    void solve();
    virtual void solveSemiImplicitCoupling();
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
    void updateGlobalLevelset( element_levelset_scalar_ptrtype & globalLevelset ) const;
    void updateFluidDensityViscosity();
    void updateInterfaceForces();
    void advectLevelsets();

    void setRebuildMatrixVector( bool b = true ) { M_doRebuildMatrixVector = b; }
    bool rebuildMatrixVector() const { return M_doRebuildMatrixVector; }

    //--------------------------------------------------------------------//
    uint16_type M_nFluids;

private:
    void updateLinear_Fluid( DataUpdateLinear & data ) const;
    void updateResidual_Fluid( DataUpdateResidual & data ) const;
    void updateJacobian_Fluid( DataUpdateJacobian & data ) const;

    void updateLinear_Levelset( levelset_model_ptrtype const& lsModel, DataUpdateLinear & data ) const;
    void updateResidual_Levelset( levelset_model_ptrtype const& lsModel, DataUpdateResidual & data ) const;
    void updateJacobian_Levelset( levelset_model_ptrtype const& lsModel, DataUpdateJacobian & data ) const;

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
    cached_levelset_scalar_field_type M_globalLevelset;
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
