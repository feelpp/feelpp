//!

#ifndef FEELPP_TOOLBOXES_MAGNETIC_HPP
#define FEELPP_TOOLBOXES_MAGNETIC_HPP

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
//#include <feel/feelvf/vf.hpp>
//#include <feel/feelts/bdf.hpp>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/modelphysics.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelcore/remeshinterpolation.hpp>

#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>

#include <feel/feelmodels/magnetic/magneticboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisVectorPotentialType>
class Magnetic : public ModelNumerical,
                 public ModelPhysics<ConvexType::nDim>
{
    typedef ModelPhysics<ConvexType::nDim> super_physics_type;
public:
    typedef ModelNumerical super_type;
    using size_type = typename super_type::size_type;
    typedef Magnetic<ConvexType,BasisVectorPotentialType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // basis
    static const uint16_type nOrderVectorPotential = BasisVectorPotentialType::nOrder;
    static const uint16_type nOrderPoly = nOrderVectorPotential;
    typedef BasisVectorPotentialType basis_vector_potential_type;
    // function space magnetic potential
    typedef FunctionSpace<mesh_type, bases<basis_vector_potential_type> > space_vector_potential_type;
    typedef std::shared_ptr<space_vector_potential_type> space_vector_potential_ptrtype;
    typedef typename space_vector_potential_type::element_type element_vector_potential_type;
    typedef std::shared_ptr<element_vector_potential_type> element_vector_potential_ptrtype;
    // materials properties
    typedef MaterialsProperties<nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;
    // exporter
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    struct FieldTag
    {
        static auto magneticPotential( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
    };

    template <typename ... Ts>
    static self_ptrtype New( Ts && ... v )
        {
            auto args = NA::make_arguments( std::forward<Ts>(v)... );
            std::string const& prefix = args.get(_prefix);
            std::string const& keyword = args.get_else(_keyword,"magnetic");
            worldcomm_ptr_t worldcomm = args.get_else(_worldcomm,Environment::worldCommPtr());
            auto && repository = args.get_else(_repository,ModelBaseRepository{});
            auto && vm = args.get_else(_vm, create_program_options( prefix ) );
            return std::make_shared<self_type>( prefix, keyword, worldcomm, "", repository, ModelBaseCommandLineOptions{vm} );
        }

    static Feel::po::options_description create_program_options( std::string const& prefix = "magnetic" ) { return magnetic_options( prefix );}

    Magnetic( std::string const& prefix,
              std::string const& keyword = "magnetic",
              worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
              ModelBaseRepository const& modelRep = ModelBaseRepository(),
              ModelBaseCommandLineOptions const& modelOptions = ModelBaseCommandLineOptions{} );

    std::shared_ptr<self_type> shared_from_this() { return std::dynamic_pointer_cast<self_type>( super_type::shared_from_this() ); }

    //___________________________________________________________________________________//
    // mesh, function space, element
    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }
    Range<mesh_type,MESH_ELEMENTS> const& rangeMeshElements() const { return M_rangeMeshElements; }

    void applyRemesh( mesh_ptrtype oldMesh, mesh_ptrtype newMesh, std::shared_ptr<RemeshInterpolation> remeshInterp = std::make_shared<RemeshInterpolation>() );

    space_vector_potential_ptrtype const& spaceVectorPotential() const { return M_spaceVectorPotential; }
    element_vector_potential_ptrtype const& fieldVectorPotentialPtr() const { return M_fieldVectorPotential; }
    element_vector_potential_type const& fieldVectorPotential() const { return *M_fieldVectorPotential; }

    //___________________________________________________________________________________//
    // physical parameters
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

    //___________________________________________________________________________________//

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initMaterialProperties();
    void initFunctionSpaces();
    void initBoundaryConditions();
    void initTimeStep();
    void initPostProcess() override;

    void initAlgebraicModel();
    void updateAlgebraicDofEliminationIds();

public :
    void initAlgebraicFactory();

    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    int nBlockMatrixGraph() const { return 1; }
    void init( bool buildModelAlgebraicFactory=true );

    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& paramValues );


    template <typename SymbolsExprType>
    void updateInitialConditions( SymbolsExprType const& se );
    //___________________________________________________________________________________//
    // execute post-processing
    //___________________________________________________________________________________//

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );

    template <typename ModelFieldsType,typename SymbolsExpr,typename ExportsExprType>
    void exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr );

    template <typename SymbolsExpr>
    void exportResults( double time, SymbolsExpr const& symbolsExpr )
        {
            return this->exportResults( time, this->modelFields(), symbolsExpr, this->exprPostProcessExports( symbolsExpr ) );
        }

    template <typename ModelFieldsType,typename SymbolsExpr, typename ModelMeasuresQuantitiesType>
    void executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ModelMeasuresQuantitiesType const& mquantities );

    bool checkResults() const override;
    //___________________________________________________________________________________//
    // export expressions
    //___________________________________________________________________________________//

    template <typename SymbExprType>
    auto exprPostProcessExportsToolbox( SymbExprType const& se, std::string const& prefix ) const
        {
#if 0
            using _expr_velocity_convection_type = std::decay_t<decltype( std::declval<ModelPhysicHeat<nDim>>().convection().expr( se ) )>;
            std::map<std::string,std::vector<std::tuple< _expr_velocity_convection_type, elements_reference_wrapper_t<mesh_type>, std::string > > > mapExprVelocityConvection;

            for ( auto const& [physicId,physicData] : this->physicsFromCurrentType() )
            {
                auto physicHeatData = std::static_pointer_cast<ModelPhysicHeat<nDim>>(physicData);
                for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicId ) )
                {
                    auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                    if ( physicHeatData->hasConvectionEnabled() )
                    {
                        auto velocityConvectionExpr = physicHeatData->convection().expr( se );
                        mapExprVelocityConvection[prefixvm(prefix,"velocity-convection")].push_back( std::make_tuple( velocityConvectionExpr, range, "nodal" ) );
                    }
                }
            }
            return hana::make_tuple( mapExprVelocityConvection );
#else
            return hana::make_tuple();
#endif
        }
    template <typename SymbExprType>
    auto exprPostProcessExports( SymbExprType const& se, std::string const& prefix = "" ) const
        {
            return hana::concat( this->materialsProperties()->exprPostProcessExports( this->mesh(),this->physicsAvailable(),se ),
                                 this->exprPostProcessExportsToolbox( se, prefix ) );
        }
    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    auto modelFields( std::string const& prefix = "" ) const
        {
            return this->modelFields( this->fieldVectorPotentialPtr(), prefix );
        }
    auto modelFields( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto field_t = this->spaceVectorPotential()->elementPtr( *sol, rowStartInVector + this->startSubBlockSpaceIndex( "vector_potential" ) );
            return this->modelFields( field_t, prefix );
        }
    auto modelFields( std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorData, std::string const& prefix = "" ) const
        {
            auto itFindSolution = vectorData.find( "solution" );
            CHECK( itFindSolution != vectorData.end() ) << "require solution data";
            vector_ptrtype sol = std::get<0>( itFindSolution->second );
            size_type rowStartInVector =  std::get<1>( itFindSolution->second );
            auto field_t = this->spaceVectorPotential()->elementPtr( *sol, rowStartInVector + this->startSubBlockSpaceIndex( "vector_potential" ) );
            return this->modelFields( field_t, prefix );
        }
    template <typename MagneticVectorPotentialFieldType>
    auto modelFields( MagneticVectorPotentialFieldType const& field_t, std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( modelField<FieldCtx::FULL>( FieldTag::magneticPotential(this), prefix, "vector_potential", field_t, "A", this->keyword() ) );
        }

    auto trialSelectorModelFields( size_type startBlockSpaceIndex = 0 ) const
        {
            return Feel::FeelModels::selectorModelFields( selectorModelField( FieldTag::magneticPotential(this), "vector_potential", startBlockSpaceIndex ) );
        }

    //___________________________________________________________________________________//
    // symbols expressions
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            auto seToolbox = this->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            auto seMeshes = this->template symbolsExprMeshes<mesh_type>();
            auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr(); // generate symbols heat_T, heat_grad_T(_x,_y,_z), heat_dn_T
            return Feel::vf::symbolsExpr( seToolbox, seParam, seMeshes, seMat, seFields );
        }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

    template <typename ModelFieldsType>
    auto symbolsExprToolbox( ModelFieldsType const& mfields ) const
        {
            return symbols_expression_empty_t{};
        }

    template <typename ModelFieldsType, typename TrialSelectorModelFieldsType>
    auto trialSymbolsExpr( ModelFieldsType const& mfields, TrialSelectorModelFieldsType const& tsmf ) const
        {
            return mfields.trialSymbolsExpr( tsmf );
        }

    //___________________________________________________________________________________//
    // model context helper
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto modelContext( ModelFieldsType const& mfields, std::string const& prefix = "" ) const
        {
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            return Feel::FeelModels::modelContext( mfields, std::move( se ) );
        }
    auto modelContext( std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        }
    auto modelContext( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( sol, rowStartInVector, prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( rowStartInVector ) );
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ), std::move( tse ) );
        }
    auto modelContextNoTrialSymbolsExpr( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( sol, rowStartInVector, prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        }

    //___________________________________________________________________________________//
    // toolbox expressions
    //___________________________________________________________________________________//
#if 0
    template <typename FieldTemperatureType, typename SymbolsExpr = symbols_expression_empty_t>
    auto normalHeatFluxExpr( FieldTemperatureType const& t, bool isOutward = true, SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
        {
            double signFlux = isOutward? -1.0 : 1.0;
            auto kappa = this->materialsProperties()->thermalConductivityExpr( symbolsExpr );
            if constexpr ( std::decay_t<decltype(kappa)>::template evaluator_t<typename mesh_type::element_type>::shape::is_scalar )
                return signFlux*kappa*gradv(t)*N();
            else
                return signFlux*inner(kappa*trans(gradv(t)),N());
        }
#endif
    //___________________________________________________________________________________//
    // apply assembly and solver
    //___________________________________________________________________________________//

    void solve();

    void updateLinearPDE( DataUpdateLinear & data ) const override;
    template <typename ModelContextType>
    void updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mfields ) const;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;
    template <typename ModelContextType>
    void updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mfields ) const;

    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
    template <typename ModelContextType>
    void updateNewtonInitialGuess( DataNewtonInitialGuess & data, ModelContextType const& mfields ) const;
    void updateJacobian( DataUpdateJacobian & data ) const override;
    template <typename ModelContextType>
    void updateJacobian( DataUpdateJacobian & data, ModelContextType const& mfields ) const;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
    void updateResidual( DataUpdateResidual & data ) const override;
    template <typename ModelContextType>
    void updateResidual( DataUpdateResidual & data, ModelContextType const& mfields ) const;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;

    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
private :
    void updateTimeStepCurrentResidual();

    auto modelMeasuresQuantities( std::string const& prefix = "" ) const
        {
            return model_measures_quantities_empty_t{};
        }

protected :

    Range<mesh_type,MESH_ELEMENTS> M_rangeMeshElements;

    space_vector_potential_ptrtype M_spaceVectorPotential;
    element_vector_potential_ptrtype M_fieldVectorPotential;

    std::map<std::string,double> M_currentParameterValues;

    // physical parameter
    materialsproperties_ptrtype M_materialsProperties;

    // boundary conditions
    using boundary_conditions_type = MagneticBoundaryConditions<nRealDim>;
    std::shared_ptr<boundary_conditions_type> M_boundaryConditions;

    std::string M_solverName;

    // post-process
    export_ptrtype M_exporter;
};


template< typename ConvexType, typename BasisVectorPotentialType>
template <typename SymbolsExprType>
void
Magnetic<ConvexType,BasisVectorPotentialType>::updateInitialConditions( SymbolsExprType const& se )
{
#if 0
    if ( !this->doRestart() )
    {
        std::vector<element_temperature_ptrtype> icTemperatureFields;
        std::map<int, double> icPriorTimes;
        if ( this->isStationary() )
        {
            icTemperatureFields = { this->fieldTemperaturePtr() };
            icPriorTimes = {{0,0}};
        }
        else
        {
            icTemperatureFields = M_bdfTemperature->unknowns();
            icPriorTimes = M_bdfTemperature->priorTimes();
        }

        super_type::updateInitialConditions( "temperature", M_rangeMeshElements, se, icTemperatureFields, icPriorTimes );

        if ( Environment::vm().count( prefixvm(this->prefix(),"initial-solution.temperature").c_str() ) )
        {
            auto myexpr = expr( soption(_prefix=this->prefix(),_name="initial-solution.temperature"),
                                "",this->worldComm(),this->repository().expr() );
            icTemperatureFields[0]->on(_range=M_rangeMeshElements,_expr=myexpr);
            for ( int k=1;k<icTemperatureFields.size();++k )
                *icTemperatureFields[k] = *icTemperatureFields[0];
        }

        if ( !this->isStationary() )
            *this->fieldTemperaturePtr() = M_bdfTemperature->unknown(0);
    }
#endif
}

template< typename ConvexType, typename BasisVectorPotentialType>
template <typename ModelFieldsType, typename SymbolsExpr, typename ExportsExprType>
void
Magnetic<ConvexType,BasisVectorPotentialType>::exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr )
{
    this->log("Magnetic","exportResults", "start");
    this->timerTool("PostProcessing").start();

    if ( M_exporter && M_exporter->exporterGeometry() == EXPORTER_GEOMETRY_CHANGE ) // TODO mv this code
        M_exporter->defaultTimeSet()->setMesh( this->mesh() );
    this->executePostProcessExports( M_exporter, time, mfields, symbolsExpr, exportsExpr );
    this->executePostProcessMeasures( time, mfields, symbolsExpr, this->modelMeasuresQuantities() );
#if 0
    this->executePostProcessSave( (this->isStationary())? invalid_uint32_type_value : M_bdfTemperature->iteration(), mfields );
#else
    this->executePostProcessSave( invalid_uint32_type_value, mfields );
#endif
    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("Magnetic","exportResults", "finish");
}

template< typename ConvexType, typename BasisVectorPotentialType>
template <typename ModelFieldsType, typename SymbolsExpr, typename ModelMeasuresQuantitiesType>
void
Magnetic<ConvexType,BasisVectorPotentialType>::executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ModelMeasuresQuantitiesType const& mquantities )
{
#if 0 // NOT COMPILE??
    // execute common post process and save measures
    super_type::executePostProcessMeasures( time, this->mesh(), M_rangeMeshElements, symbolsExpr, mfields, mquantities );
#endif
}

} // namespace FeelModels
} // namespace Feel

#include <feel/feelmodels/magnetic/magneticassembly.hpp>


#endif /* FEELPP_TOOLBOXES_HEAT_HPP */
