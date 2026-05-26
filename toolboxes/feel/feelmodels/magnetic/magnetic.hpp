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
    using self_type = Magnetic<ConvexType,BasisVectorPotentialType>;
    using self_ptrtype = std::shared_ptr<self_type>;
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static constexpr uint16_type nDim = convex_type::nDim;
    static constexpr uint16_type nOrderGeo = convex_type::nOrder;
    static constexpr uint16_type nRealDim = convex_type::nRealDim;
    using mesh_type = Mesh<convex_type>;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    // basis
    typedef BasisVectorPotentialType basis_vector_potential_type;
    // function space magnetic vector potential
    using space_vector_potential_type = FunctionSpace<mesh_type, bases<basis_vector_potential_type>>;
    using space_vector_potential_ptrtype = std::shared_ptr<space_vector_potential_type>;
    using element_vector_potential_type = typename space_vector_potential_type::element_type;
    using element_vector_potential_ptrtype = std::shared_ptr<element_vector_potential_type>;
    // function space lagrange multiplier for Coulomb gauge
    using space_lm_coulombgauge_type = FunctionSpace<mesh_type, bases<Lagrange<1,Scalar,Continuous,PointSetFekete>> >;
    using space_lm_coulombgauge_ptrtype = std::shared_ptr<space_lm_coulombgauge_type>;
    using element_lm_coulombgauge_type = typename space_lm_coulombgauge_type::element_type;
    using element_lm_coulombgauge_ptrtype = std::shared_ptr<element_lm_coulombgauge_type>;

    // materials properties
    using materialsproperties_type = MaterialsProperties<nRealDim>;
    using materialsproperties_ptrtype = std::shared_ptr<materialsproperties_type>;
    // exporter
    using export_type = Exporter<mesh_type,nOrderGeo>;
    using export_ptrtype = std::shared_ptr<export_type>;

    struct FieldTag
    {
        static auto vectorPotential( self_type const* t ) { return ModelFieldTag<self_type,0,"vector_potential">( t ); }
        static auto lagrangeMultiplierCoulombGauge( self_type const* t ) { return ModelFieldTag<self_type,1,"lagrange_multiplier_CoulombGauge">( t ); }
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
            return std::make_shared<self_type>( prefix, keyword, worldcomm, repository, ModelBaseCommandLineOptions{vm} );
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

    space_lm_coulombgauge_ptrtype spaceLagrangeMultiplierCoulombGauge() const { return M_spaceLagrangeMultiplierCoulombGauge; }
    element_lm_coulombgauge_ptrtype fieldLagrangeMultiplierCoulombGaugePtr() const { return M_fieldLagrangeMultiplierCoulombGauge; }
    element_lm_coulombgauge_type const& fieldLagrangeMultiplierCoulombGauge() const { return *M_fieldLagrangeMultiplierCoulombGauge; }

    //___________________________________________________________________________________//
    // physical parameters
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

    //___________________________________________________________________________________//
    // time step scheme (TODO)
    std::shared_ptr<TSBase> timeStepBase() const { return {}; }
    void startTimeStep() {}
    void updateTimeStep() {}

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
    //int nBlockMatrixGraph() const { return 1; }
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
            auto const& A = this->fieldVectorPotential();
            auto B = this->fluxDensityExpr();
            using _expr_flux_density_type = std::decay_t<decltype( B )>;
            using _expr_field_intensity_type = std::decay_t<decltype( this->fieldIntensityExpr( A, "", se ) )>;
            using _range_mesh_type = Range<mesh_type,MESH_ELEMENTS>;
            std::map<std::string,std::vector<std::tuple< _expr_flux_density_type, _range_mesh_type, std::string > > > mapExprFluxDensity;
            std::map<std::string,std::vector<std::tuple< _expr_field_intensity_type, _range_mesh_type, std::string > > > mapExprFieldIntensity;
            mapExprFluxDensity[prefixvm(prefix,"flux-density")].push_back( std::make_tuple( B, M_rangeMeshElements, "element" ) );

            for ( auto const& [physicId,physicData] : this->physicsFromCurrentType() )
            {
                auto physicMagneticData = std::static_pointer_cast<ModelPhysicMagnetic<nDim>>(physicData);
                auto mu_0 = physicMagneticData->vacuumPermeabilityExpr();
                for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicId ) )
                {
                    auto const& matRange = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                    mapExprFieldIntensity[prefixvm(prefix,"field-intensity")].push_back( std::make_tuple( this->fieldIntensityExpr( A, matName, se ), matRange, "element" ) );
                }
            }

            return hana::make_tuple( mapExprFluxDensity, mapExprFieldIntensity );
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
            return this->modelFields( this->fieldVectorPotentialPtr(), this->fieldLagrangeMultiplierCoulombGaugePtr(), prefix );
        }
    auto modelFields( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto field_A = this->spaceVectorPotential()->elementPtr( *sol, rowStartInVector + this->startSubBlockSpaceIndex( FieldTag::vectorPotential(this).identifier() ) );
            element_lm_coulombgauge_ptrtype field_lmcg;
            if ( this->spaceLagrangeMultiplierCoulombGauge() && M_nullSpaceMethod == "saddle-point" )
                    field_lmcg = this->spaceLagrangeMultiplierCoulombGauge()->elementPtr(
                        *sol, rowStartInVector + this->startSubBlockSpaceIndex( FieldTag::lagrangeMultiplierCoulombGauge(this).identifier() ) );
            return this->modelFields( field_A, field_lmcg, prefix );
        }
    auto modelFields( std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorData, std::string const& prefix = "" ) const
        {
            auto itFindSolution = vectorData.find( "solution" );
            CHECK( itFindSolution != vectorData.end() ) << "require solution data";
            vector_ptrtype sol = std::get<0>( itFindSolution->second );
            size_type rowStartInVector =  std::get<1>( itFindSolution->second );
            auto field_A = this->spaceVectorPotential()->elementPtr( *sol, rowStartInVector + this->startSubBlockSpaceIndex( FieldTag::vectorPotential(this).identifier() ) );
            element_lm_coulombgauge_ptrtype field_lmcg;
            if ( this->spaceLagrangeMultiplierCoulombGauge() && M_nullSpaceMethod == "saddle-point" )
                    field_lmcg = this->spaceLagrangeMultiplierCoulombGauge()->elementPtr(
                        *sol, rowStartInVector + this->startSubBlockSpaceIndex( FieldTag::lagrangeMultiplierCoulombGauge(this).identifier() ) );
            return this->modelFields( field_A, field_lmcg, prefix );
        }
    template <typename MagneticVectorPotentialFieldType,typename LagrangeMultiplierCoulombGaugeFieldType>
    auto modelFields( MagneticVectorPotentialFieldType const& field_A, LagrangeMultiplierCoulombGaugeFieldType const& field_lmcg,
                      std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields(
                modelField<FieldCtx::FULL>( FieldTag::vectorPotential(this), prefix, FieldTag::vectorPotential(this).identifierString(), field_A, "A", this->keyword() ),
                modelField<FieldCtx::ID>( FieldTag::lagrangeMultiplierCoulombGauge(this), prefix, FieldTag::lagrangeMultiplierCoulombGauge(this).identifierString(), field_lmcg, "lmcg", this->keyword() )
                                                 );
        }

    auto trialSelectorModelFields( size_type startBlockSpaceIndex = 0 ) const
        {
            return Feel::FeelModels::selectorModelFields( selectorModelField( FieldTag::vectorPotential(this), FieldTag::vectorPotential(this).identifierString(), startBlockSpaceIndex ) );
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

    auto fluxDensityExpr() const { return this->fluxDensityExpr( this->fieldVectorPotential() ); }

    template <typename FieldVectorPotentialType>
    auto fluxDensityExpr( FieldVectorPotentialType const& A ) const
        {
          return curlv(A);
        }

    template <typename SymbolsExpr = symbols_expression_empty_t>
    auto reluctivityExpr( std::string const& matName, SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
    {
        auto const& magneticRelativePermeability = this->materialsProperties()->materialProperty( matName, "magnetic-relative-permeability" );
        double mu_0 = ModelPhysicMagnetic<nDim>::vacuumPermeabilityConstant();
#if 0
        if ( magneticRelativePermeability.isMatrix() )
        {
          if constexpr (  nDim == 3 )
            {
              auto mu_r = expr( magneticRelativePermeability.template expr<nDim,nDim>(), symbolsExpr );
              // TODO
            }
        }

        auto mu_r = expr( magneticRelativePermeability.expr(), symbolsExpr );

        return 1./(mu_0*mu_r);
#else
        if constexpr ( nDim == 2 && true )
        {
            // case scalar only
            auto mu_r = expr( magneticRelativePermeability.expr(), symbolsExpr );
            return 1./(mu_0*mu_r);
        }
        else
        {
            auto Id = eye<nDim,nDim>();
            using expr_mu_r_scalar_type = std::decay_t<decltype( expr( magneticRelativePermeability.expr(), symbolsExpr )*Id )>;
            using expr_mu_r_matrix_type = std::decay_t<decltype( expr( magneticRelativePermeability.template expr<nDim,nDim>(), symbolsExpr ) )>;
            auto mu_r = exprOptionalConcat<expr_mu_r_scalar_type,expr_mu_r_matrix_type>();
            if ( magneticRelativePermeability.isMatrix() )
                mu_r.expression().add( expr( magneticRelativePermeability.template expr<nDim,nDim>(), symbolsExpr ) );
            else
                mu_r.expression().add( expr( magneticRelativePermeability.expr(), symbolsExpr )*Id );
            return (1./mu_0)*inv(mu_r);
        }
#endif
    }
    template <typename SymbolsExpr = symbols_expression_empty_t>
    auto reluctivityExpr( SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
        {
            double mu_0 = ModelPhysicMagnetic<nDim>::vacuumPermeabilityConstant();
            static constexpr bool is2DPlanarTransverseMagnetic = nDim == 2 && false; // unknown is scalar
            static constexpr bool is2DPlanarTransverseElectric = nDim == 2 && true; // unknown is vector
            if constexpr ( is2DPlanarTransverseElectric )
            {
                // case scalar only
                auto mu_r = this->materialsProperties()->template materialPropertyExpr<1,1>( "magnetic-relative-permeability", symbolsExpr );
                return 1./(mu_0*mu_r);
            }
            else
            {
                // scalar or matrix, so return always matrix shape (scalar is mutliplied by Identity matrix)
                auto mu_r = this->materialsProperties()->template materialPropertyExprScalarOrMatrix<nDim>( "magnetic-relative-permeability", symbolsExpr );
                return (1./mu_0)*inv(mu_r);
            }
        }


    template <typename SymbolsExpr = symbols_expression_empty_t>
    auto fieldIntensityExpr( std::string const& matName, SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
        {
            return this->fieldIntensityExpr( this->fieldVectorPotential(), matName, symbolsExpr );
        }
    template <typename FieldVectorPotentialType, typename SymbolsExpr = symbols_expression_empty_t>
    auto fieldIntensityExpr( FieldVectorPotentialType const& A, std::string const& matName, SymbolsExpr const& symbolsExpr = symbols_expression_empty_t{} ) const
        {
          return this->reluctivityExpr( matName, symbolsExpr )*this->fluxDensityExpr( A );
        }

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

    void initInHousePreconditioner();
    void updateInHousePreconditioner( DataUpdateLinear & data ) const override;
    void updateInHousePreconditioner( DataUpdateJacobian & data ) const override;

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

    space_lm_coulombgauge_ptrtype M_spaceLagrangeMultiplierCoulombGauge;
    element_lm_coulombgauge_ptrtype M_fieldLagrangeMultiplierCoulombGauge;

    std::map<std::string,double> M_currentParameterValues;

    // physical parameter
    materialsproperties_ptrtype M_materialsProperties;

    // boundary conditions
    using boundary_conditions_type = MagneticBoundaryConditions<nRealDim>;
    std::shared_ptr<boundary_conditions_type> M_boundaryConditions;

    std::string M_solverName;
    std::string M_nullSpaceMethod = "regularized-formulation"; // "regularized-formulation", "saddle-point"
    bool M_preconditionerAttachAms = false;
    sparse_matrix_ptrtype M_preconditionerAmsMatrixG;
    std::array<vector_ptrtype,nRealDim> M_preconditionerAmsVectorOnes;

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
    // execute common post process and save measures
    super_type::executePostProcessMeasures( time, this->mesh(), M_rangeMeshElements, symbolsExpr, mfields, mquantities );
}

} // namespace FeelModels
} // namespace Feel

#include <feel/feelmodels/magnetic/magneticassembly.hpp>


#endif /* FEELPP_TOOLBOXES_HEAT_HPP */
