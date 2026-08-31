//!

#ifndef FEELPP_TOOLBOXES_ELECTROMAGNETIC_HPP
#define FEELPP_TOOLBOXES_ELECTROMAGNETIC_HPP 1

#include <feel/feelmodels/electric/electric.hpp>
#include <feel/feelmodels/magnetic/magnetic.hpp>


namespace Feel
{
namespace FeelModels
{

template< typename ElectricType, typename MagneticType>
class Electromagnetic : public ModelNumerical,
                        public ModelPhysics<ElectricType::convex_type::nDim>
{
    using super_physics_type = ModelPhysics<ElectricType::convex_type::nDim>;
public:
    using super_type = ModelNumerical;
    using self_type = Electromagnetic<ElectricType,MagneticType>;
    using self_ptrtype = std::shared_ptr<self_type>;

    using electric_model_type = ElectricType;
    using electric_model_ptrtype = std::shared_ptr<electric_model_type>;

    using magnetic_model_type = MagneticType;
    using magnetic_model_ptrtype = std::shared_ptr<magnetic_model_type>;

    // mesh
    using mesh_electric_type = typename electric_model_type::mesh_type;
    using mesh_magnetic_type = typename magnetic_model_type::mesh_type;
    using mesh_type = mesh_electric_type;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

    using materialsproperties_type = MaterialsProperties<mesh_type::nRealDim>;
    using materialsproperties_ptrtype = std::shared_ptr<materialsproperties_type>;

    // exporter
    using export_type = Exporter<mesh_type,mesh_type::nOrder>;
    using export_ptrtype = std::shared_ptr<export_type>;

    //___________________________________________________________________________________//

    template <typename ... Ts>
    static self_ptrtype New( Ts && ... v )
        {
            auto args = NA::make_arguments( std::forward<Ts>(v)... );
            std::string const& prefix = args.get(_prefix);
            std::string const& keyword = args.get_else(_keyword,"magnetic");
            worldcomm_ptr_t worldcomm = args.get_else(_worldcomm,Environment::worldCommPtr());
            auto && repository = args.get_else(_repository,ModelBaseRepository{});
            auto && vm = args.get_else(_vm, create_program_options( prefix ) );
            return std::make_shared<self_type>( prefix, keyword, worldcomm, repository, ModelBaseCommandLineOptions{vm} );
        }

    static Feel::po::options_description create_program_options( std::string const& prefix = "electromagnetic" ) { return electromagnetic_options( prefix );}

    // constructor
    Electromagnetic( std::string const& prefix,
                     std::string const& keyword = "electromagnetic",
                     worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                     ModelBaseRepository const& modelRep = ModelBaseRepository(),
                     ModelBaseCommandLineOptions const& modelOptions = ModelBaseCommandLineOptions{} );

    std::shared_ptr<self_type> shared_from_this() { return std::dynamic_pointer_cast<self_type>( super_type::shared_from_this() ); }

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;


private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initPostProcess() override;

    void updatePhysics( typename super_physics_type::PhysicsTreeNode & physicsTree, ModelModels const& models ) override;

public :
    // update for use
    void init( bool buildModelAlgebraicFactory = true );

    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    //int nBlockMatrixGraph() const;

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );

    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& paramValues );

    //___________________________________________________________________________________//

    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }
    //Range<mesh_type,MESH_ELEMENTS> const& rangeMeshElements() const { return M_rangeMeshElements; }

    electric_model_ptrtype electricModel() const { return M_electricModel; }
    magnetic_model_ptrtype magneticModel() const { return M_magneticModel; }

    materialsproperties_ptrtype materialsProperties() const { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

    //___________________________________________________________________________________//
    // time step scheme (TODO)
    std::shared_ptr<TSBase> timeStepBase() const { return {}; }
    // std::shared_ptr<TSBase> timeStepBase() { return this->heatModel()->timeStepBase(); }
    // std::shared_ptr<TSBase> timeStepBase() const { return this->heatModel()->timeStepBase(); }
    void startTimeStep();
    void updateTimeStep();

    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    auto modelFields( std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( this->electricModel()->modelFields( prefixvm( prefix, this->electricModel()->keyword() ) ),
                                                  this->magneticModel()->modelFields( prefixvm( prefix, this->magneticModel()->keyword() ) ),
                                                  this->template modelFieldsMeshes<mesh_type>( prefix ) );
        }
    auto modelFields( vector_ptrtype sol, size_type rowStartInVectorElectric, size_type rowStartInVectorMagnetic, std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( this->electricModel()->modelFields( sol, rowStartInVectorElectric, prefixvm( prefix,this->electricModel()->keyword() ) ),
                                                  this->magneticModel()->modelFields( sol, rowStartInVectorMagnetic, prefixvm( prefix,this->magneticModel()->keyword() ) ),
                                                  this->template modelFieldsMeshes<mesh_type>( prefix ) );
        }
    template <typename ModelFieldsElectricType,typename ModelFieldsMagneticType>
    auto modelFields( ModelFieldsElectricType const& mfieldsElectric, ModelFieldsMagneticType const& mfieldsMagnetic, std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( mfieldsElectric, mfieldsMagnetic, this->template modelFieldsMeshes<mesh_type>( prefix ) );
        }

    auto trialSelectorModelFields( size_type startBlockSpaceIndexElectric, size_type startBlockSpaceIndexMagnetic ) const
        {
            return Feel::FeelModels::selectorModelFields( this->electricModel()->trialSelectorModelFields( startBlockSpaceIndexElectric ),
                                                          this->magneticModel()->trialSelectorModelFields( startBlockSpaceIndexMagnetic ) );
        }


    //___________________________________________________________________________________//
    // model context helper
    //___________________________________________________________________________________//

    auto modelContext( std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        }
    auto modelContext( vector_ptrtype sol, size_type startBlockSpaceIndexElectric, size_type startBlockSpaceIndexMagnetic, std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( sol, startBlockSpaceIndexElectric, startBlockSpaceIndexMagnetic, prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            auto tse =  this->trialSymbolsExpr( mfields, trialSelectorModelFields( startBlockSpaceIndexElectric, startBlockSpaceIndexMagnetic ) );
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ), std::move( tse ) );
        }

    //___________________________________________________________________________________//
    // symbols expressions
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            auto seElectric = this->electricModel()->symbolsExprToolbox( mfields );
            auto seMagnetic = this->magneticModel()->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            auto seMeshes = this->template symbolsExprMeshes<mesh_type,false>();
            auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr();
            auto sePhysics = this->symbolsExprPhysics( this->physics() );
            return Feel::vf::symbolsExpr( seElectric,seMagnetic,seParam,seMeshes,seMat,seFields,sePhysics );
        }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

    template <typename ModelFieldsType, typename TrialSelectorModelFieldsType>
    auto trialSymbolsExpr( ModelFieldsType const& mfields, TrialSelectorModelFieldsType const& tsmf ) const
        {
            return mfields.trialSymbolsExpr( tsmf );
        }

    //___________________________________________________________________________________//
    // apply assembly and solver
    void solve();


    void updateLinearPDE( DataUpdateLinear & data ) const override;
    template <typename ModelContextType>
    void updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mfields ) const;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;
    template <typename ModelContextType>
    void updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mfields ) const;

#if 0
    void updateLinearPDE( DataUpdateLinear & data ) const override;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;

    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
    void updateJacobian( DataUpdateJacobian & data ) const override;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
    void updateResidual( DataUpdateResidual & data ) const override;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;
#endif
    //___________________________________________________________________________________//

    bool checkResults() const override
        {
            // several calls (not do in on line) to be sure that all check have been run
            bool checkElectromagnetic = super_type::checkResults();
            bool checkElectric = this->electricModel()->checkResults();
            bool checkMagnetic = this->magneticModel()->checkResults();
            return checkElectromagnetic && checkElectric && checkMagnetic;
        }

private :
    electric_model_ptrtype M_electricModel;
    magnetic_model_ptrtype M_magneticModel;

    // physical parameter
    // std::string M_modelName;
    // bool M_modelUseJouleEffect;
    materialsproperties_ptrtype M_materialsProperties;

    // solver
    std::string M_solverName;
    // bool M_solverNewtonInitialGuessUseLinearElectromagnetic,M_solverNewtonInitialGuessUseLinearHeat,M_solverNewtonInitialGuessUseLinearElectric;

    // post-process
    export_ptrtype M_exporter;
};

} // namespace FeelModels
} // namespace Feel

#include <feel/feelmodels/electromagnetic/electromagneticassembly.hpp>

#endif // FEELPP_TOOLBOXES_ELECTROMAGNETIC_HPP
