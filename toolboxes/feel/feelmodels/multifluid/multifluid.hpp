#ifndef _MULTIFLUID_HPP
#define _MULTIFLUID_HPP 1

#include <algorithm>
#include <fmt/core.h>

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
    public ModelPhysics<FluidType::convex_type::nDim>,
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
    typedef MaterialsProperties<mesh_type::nRealDim> materials_properties_type;
    typedef std::shared_ptr<materials_properties_type> materials_properties_ptrtype;

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

    //--------------------------------------------------------------------//
    // Lagrange P1 iso-Pn
    typedef OperatorLagrangeP1<component_space_fluid_velocity_type> op_lagrangeP1_type;
    typedef std::shared_ptr<op_lagrangeP1_type> op_lagrangeP1_ptrtype;

    ////--------------------------------------------------------------------//
    //// Density/viscosity
    //typedef typename fluid_model_type::materialsproperties_type materials_properties_type;
    //typedef typename fluid_model_type::materialsproperties_ptrtype materials_properties_ptrtype;
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

    //--------------------------------------------------------------------//
    // Initialization
    void init( bool buildModelAlgebraicFactory = true );
    void initAlgebraicFactory();

    void loadParametersFromOptionsVm();

    template <typename SymbolsExprType>
    void updateInitialConditions( SymbolsExprType const& se );

    //--------------------------------------------------------------------//
    std::string globalLevelsetPrefix() const { return prefixvm( this->prefix(), "levelset"); }

    std::shared_ptr<std::ostringstream> getInfo() const override;

    //--------------------------------------------------------------------//
    // Mesh
    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }

    //--------------------------------------------------------------------//
    // Models
    fluid_model_ptrtype const& fluidModel() const { return M_fluidModel; }
    fluid_model_ptrtype fluidModel() { return M_fluidModel; }

    std::vector<levelset_model_ptrtype> const& levelsetModels() const { return M_levelsetModels; }
    std::vector<levelset_model_ptrtype> & levelsetModels() { return M_levelsetModels; }
    levelset_model_ptrtype const& levelsetModel( index_type i ) const { return M_levelsetModels.at(i); }
    size_type nLevelsets() const { return M_levelsetModels.size(); }

    // Global levelset
    auto globalLevelsetExpr() const {
        std::vector< element_levelset_scalar_ptrtype > levelsets;
        std::transform( this->levelsetModels().begin(), this->levelsetModels().end(), std::back_inserter(levelsets),
                []( levelset_model_ptrtype const& lsModel ) { return lsModel->phi(); }
                );
        return Feel::FeelModels::globalLevelsetExpr( levelsets );
    }
    element_levelset_scalar_ptrtype globalLevelsetElt() const { return M_globalLevelset.fieldPtr(); }

    //--------------------------------------------------------------------//
    // Function spaces
    space_levelset_ptrtype const& functionSpaceLevelset() const { return M_levelsetSpaceManager->functionSpaceScalar(); }
    space_levelset_vectorial_ptrtype const& functionSpaceLevelsetVectorial() const { return M_levelsetSpaceManager->functionSpaceVectorial(); }
    space_inextensibilitylm_ptrtype const& functionSpaceInextensibilityLM() const;
    //--------------------------------------------------------------------//
    // Operator Lagrange P1
    bool useLagrangeP1iso() const { return M_useLagrangeP1iso; }
    op_lagrangeP1_ptrtype opLagrangeP1() const { return M_opLagrangeP1iso; }

    // Physical parameters
    materials_properties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materials_properties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materials_properties_ptrtype mp ) { M_materialsProperties = mp; }
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
    // Model fields
    struct FieldTag
    {
        static auto globalLevelsetElt( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
    };
    auto modelFieldsLevelsets( std::string const& prefix = "" ) const
    {
        using mfields_levelset_type = std::decay_t<decltype(this->levelsetModel(0)->modelFields() )>;
        mfields_levelset_type mfieldsLevelsets;
        //for ( auto it = M_levelsetModels.begin() ; it != M_levelsetModels.end() ; ++it )
            //mfieldsLevelsets = Feel::FeelModels::modelFields( 
                    //mfieldsLevelsets, 
                    //it->second->modelFields(  prefixvm( prefix,  it->second->keyword() ) ) 
                    //);
        for ( levelset_model_ptrtype const& lsModel : this->levelsetModels() )
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
        using mfields_levelset_type = std::decay_t<decltype(this->levelsetModel(0)->modelFields( sol ) )>;
        mfields_levelset_type mfieldsLevelsets;
        for ( index_type i = 0; i < this->levelsetModels().size(); ++i )
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
        using mfields_levelset_type = std::decay_t<decltype(this->levelsetModel(0)->modelFields( sols[0] ) )>;
        mfields_levelset_type mfieldsLevelsets;
        for ( index_type i = 0; i < this->levelsetModels().size(); ++i )
            mfieldsLevelsets = Feel::FeelModels::modelFields(
                    mfieldsLevelsets,
                    this->levelsetModel(i)->modelFields( sols[i], rowStartInVectorLevelsets[i], this->levelsetModel(i)->keyword() )
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
        using tsmfields_levelset_type = std::decay_t<decltype(this->levelsetModel(0)->trialSelectorModelFields() )>;
        tsmfields_levelset_type tsmfieldsLevelsets;
        for ( index_type i = 0; i < this->levelsetModels().size(); ++i )
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
        return Feel::vf::symbolsExpr( seFluid, seParam, seMeshes, seMat, seFields, sePhysics );
    }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

    template <typename ModelFieldsType, typename TrialSelectorModelFieldsType>
    auto trialSymbolsExpr( ModelFieldsType const& mfields, TrialSelectorModelFieldsType const& tsmf ) const
    {
        return mfields.trialSymbolsExpr( tsmf );
    }

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
        auto mfields = this->modelFields( sol, startBlockSpaceIndexFluid, startBlockSpaceIndexLevelsets, prefix );
        auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
        auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( startBlockSpaceIndexFluid, startBlockSpaceIndexLevelsets ) );
        return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ), std::move( tse ) );
    }
    auto modelContext( vector_ptrtype solFluid, size_type startBlockSpaceIndexFluid, std::vector<vector_ptrtype> solLevelsets, std::vector<size_type> startBlockSpaceIndexLevelsets, std::string const& prefix = "" ) const
    {
        auto mfields = this->modelFields( solFluid, startBlockSpaceIndexFluid, solLevelsets, startBlockSpaceIndexLevelsets, prefix );
        auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
        auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( startBlockSpaceIndexFluid, startBlockSpaceIndexLevelsets ) );
        return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ), std::move( tse ) );
    }
    //auto modelContextNoTrialSymbolsExpr( vector_ptrtype sol, size_type startBlockSpaceIndexFluid, std::vector<size_type> startBlockSpaceIndexLevelsets, std::string const& prefix = "" ) const
    //{
        //return this->modelContextNoTrialSymbolsExpr( sol, startBlockSpaceIndexFluid, sol, startBlockSpaceIndexLevelsets, prefix );
    //}
    //auto modelContextNoTrialSymbolsExpr( vector_ptrtype solHeat, size_type startBlockSpaceIndexHeat, vector_ptrtype solFluid, size_type startBlockSpaceIndexFluid, std::string const& prefix = "" ) const
    //{
        //// auto mfields = this->modelFields( solHeat, startBlockSpaceIndexHeat, solFluid, startBlockSpaceIndexFluid, prefix );
        //// auto se = this->symbolsExpr( mfields );
        //// return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        //return this->modelContextNoTrialSymbolsExpr( { { "solution", std::make_tuple( solHeat, startBlockSpaceIndexHeat ) } },
                //{ { "solution", std::make_tuple( solFluid, startBlockSpaceIndexFluid ) } },
                //prefix );
    //}
    //auto modelContextNoTrialSymbolsExpr( std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorDataHeat,
                                         //std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorDataFluid,
                                         //std::string const& prefix = "" ) const
    //{
        //auto mfields = this->modelFields( vectorDataHeat, vectorDataFluid, prefix );
        //auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
        //return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
    //}
    auto modelContext( vector_ptrtype sol, fluid_model_ptrtype const& fluidModel, std::vector<levelset_model_ptrtype> const& levelsetModels, std::string const& prefix = "" ) const
    {
        std::vector<size_type> startBlockSpaceIndexLevelsets;
        std::transform( levelsetModels.begin(), levelsetModels.end(), std::back_inserter( startBlockSpaceIndexLevelsets ),
                []( levelset_model_ptrtype const& lsModel ) { return lsModel->startBlockSpaceIndexVector(); } );
        return this->modelContext( sol, fluidModel->startBlockSpaceIndexVector(), startBlockSpaceIndexLevelsets, prefix );
    }

    //--------------------------------------------------------------------//
    // Algebraic data
    int nBlockMatrixGraph() const;
    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;

    //--------------------------------------------------------------------//
    double globalLevelsetThicknessInterface() const { return M_globalLevelsetThicknessInterface; }
    //--------------------------------------------------------------------//
    bool hasInextensibility( index_type i ) const { return M_hasInextensibility.at(i); }
    std::string const& inextensibilityMethod( index_type i ) const { return M_inextensibilityMethod.at(i); }
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

    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
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

    void buildBlockVectorSolution();
    virtual int initBlockVectorSolution();

    bool useImplicitCoupling() const;

    //--------------------------------------------------------------------//
    void updateGlobalLevelset( element_levelset_scalar_ptrtype & globalLevelset ) const;
    void updateInterfaceForces();

private:
    void updateLinear_Fluid( DataUpdateLinear & data ) const;
    void updateResidual_Fluid( DataUpdateResidual & data ) const;
    void updateJacobian_Fluid( DataUpdateJacobian & data ) const;

    void updateLinear_Levelset( size_type lsModelIndex, DataUpdateLinear & data ) const;
    void updateResidual_Levelset( size_type lsModelIndex, DataUpdateResidual & data ) const;
    void updateJacobian_Levelset( size_type lsModelIndex, DataUpdateJacobian & data ) const;

private:
    //--------------------------------------------------------------------//
    // Lagrange P1 iso-Pn
    bool M_useLagrangeP1iso;
    op_lagrangeP1_ptrtype M_opLagrangeP1iso;
    //--------------------------------------------------------------------//
    mesh_ptrtype M_mesh;

    fluid_model_ptrtype M_fluidModel;
    size_type M_nLevelsets;
    std::vector<levelset_model_ptrtype> M_levelsetModels;

    levelset_space_manager_ptrtype M_levelsetSpaceManager;
    levelset_tool_manager_ptrtype M_levelsetToolManager;

    //--------------------------------------------------------------------//
    materials_properties_ptrtype M_materialsProperties;

    //--------------------------------------------------------------------//
    cached_levelset_scalar_field_type M_globalLevelset;
    mutable bool M_doUpdateGlobalLevelset;

    exporter_ptrtype M_globalLevelsetExporter;

    //--------------------------------------------------------------------//
    // Solve
    bool M_usePicardIterations;

    //--------------------------------------------------------------------//
    // Parameters
    materials_properties_ptrtype M_fluidMaterialProperties;
    std::map<std::string, materials_properties_ptrtype> M_levelsetsMaterialProperties;
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
    std::vector<bool> M_hasInextensibility;
    bool M_hasInextensibilityLM;
    std::vector<std::string> M_inextensibilityMethod;
    std::vector<element_levelset_scalar_ptrtype> M_inextensibleLevelsets;

    mutable range_elements_type M_rangeInextensibilityLM;
    mutable space_inextensibilitylm_ptrtype M_spaceInextensibilityLM;
    mutable bool M_doUpdateInextensibilityLM;
    // Penalty method gamma
    std::vector<double> M_inextensibilityGamma;
    //--------------------------------------------------------------------//
    // Redistanciation
    std::vector<int> M_levelsetRedistEvery;

    //--------------------------------------------------------------------//
    std::vector<vector_ptrtype> M_algebraicBlockVectorSolutionLevelsets;

    //--------------------------------------------------------------------//
    // Post-process
    exporter_ptrtype M_exporter;
    // fluid exports
    std::set<std::string> M_postProcessFieldsExportedFluid;
    // levelset forces
    std::vector< std::string > M_postProcessMeasuresLevelsetForces;
};

template< typename FluidType, typename LevelSetType>
template <typename SymbolsExprType>
void
MultiFluid<FluidType, LevelSetType>::updateInitialConditions( SymbolsExprType const& se )
{
    M_fluidModel->updateInitialConditions( se );
    for( levelset_model_ptrtype const& lsModel: this->levelsetModels() )
        lsModel->updateInitialConditions( /*se*/ );
}
        

} // namespace FeelModels
} // namespace Feel

#endif
