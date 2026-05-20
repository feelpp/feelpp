//!

#include <feel/feelmodels/magnetic/magnetic.hpp>
#include <feel/feeldiscr/pch.hpp>

namespace Feel
{
namespace FeelModels
{

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
MAGNETIC_CLASS_TEMPLATE_TYPE::Magnetic( std::string const& prefix,
                                        std::string const& keyword,
                                        worldcomm_ptr_t const& worldComm,
                                        ModelBaseRepository const& modelRep,
                                        ModelBaseCommandLineOptions const& modelOptions )
    :
    super_type( prefix, keyword, worldComm, "", modelRep, modelOptions ),
    ModelPhysics<nDim>( "magnetic" ),
    ModelBase( prefix, keyword, worldComm, "", modelRep, modelOptions )
{
    this->log("Magnetic","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MagneticConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MagneticSolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MagneticPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MagneticTimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    this->log("Magnetic","constructor", "finish");

}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    M_solverName = soption(_name="solver",_prefix=this->prefix(),_vm=this->clovm());
    M_nullSpaceMethod = soption(_name="null-space.method",_prefix=this->prefix(),_vm=this->clovm());
    if ( M_nullSpaceMethod != "regularized-formulation" && M_nullSpaceMethod != "saddle-point" && M_nullSpaceMethod != "none" )
        throw std::runtime_error( "null-space.method should be regularized-formulation, saddle-point or none" );
    M_preconditionerAttachAms = boption(_name="preconditioner.attach-ams",_prefix=this->prefix(),_vm=this->clovm());
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::initMesh()
{
    this->log("Magnetic","initMesh", "start");
    this->timerTool("Constructor").start();

    if ( this->modelProperties().jsonData().contains("Meshes") )
        super_type::super_model_meshes_type::setup( this->modelProperties().jsonData().at("Meshes"), {this->keyword()} );
    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    super_type::super_model_meshes_type::modelMesh( this->keyword() ).setFunctionApplyRemesh(
        [this]( typename super_type::super_model_meshes_type::mesh_base_ptrtype mnew,
                typename super_type::super_model_meshes_type::mesh_base_ptrtype mold ) { this->applyRemesh( std::dynamic_pointer_cast<mesh_type>( mold ),
                                                                                                            std::dynamic_pointer_cast<mesh_type>( mnew ) ); }
                                                                                             );

    CHECK( this->mesh() ) << "mesh generation fail";
    this->log("Magnetic","initMesh", fmt::format("mesh numGlobalElements : {}", this->mesh()->numGlobalElements()));

    double tElpased = this->timerTool("Constructor").stop("initMesh");
    this->log("Magnetic","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::initMaterialProperties()
{
    this->log("Magnetic","initMaterialProperties", "start" );
    this->timerTool("Constructor").start();

    if ( !M_materialsProperties )
    {
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    double tElpased = this->timerTool("Constructor").stop("initMaterialProperties");
    this->log("Magnetic","initMaterialProperties",(boost::format("finish in %1% s")%tElpased).str() );
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::initFunctionSpaces()
{
    this->log("Magnetic","initFunctionSpaces", "start" );
    this->timerTool("Constructor").start();

    auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
    // functionspace
    if ( mom->isDefinedOnWholeMesh( this->physicsAvailableFromCurrentType() ) )
    {
        M_rangeMeshElements = elements(this->mesh());
        M_spaceVectorPotential = space_vector_potential_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm() );
    }
    else
    {
        M_rangeMeshElements = markedelements(this->mesh(), mom->markers( this->physicsAvailableFromCurrentType() ));
        M_spaceVectorPotential = space_vector_potential_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),_range=M_rangeMeshElements );
    }
    this->log("Magnetic","initFunctionSpaces", fmt::format("vector_potential space ndof : {}",M_spaceVectorPotential->nDof()) );

    M_fieldVectorPotential = M_spaceVectorPotential->elementPtr( FieldTag::vectorPotential(this).identifierString() );

    if ( M_nullSpaceMethod == "saddle-point" )
    {
        M_spaceLagrangeMultiplierCoulombGauge = space_lm_coulombgauge_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),_range=this->rangeMeshElements() );
        this->log("Magnetic","initFunctionSpaces", fmt::format("lagrange_multiplier_CoulombGauge space ndof : {}",M_spaceLagrangeMultiplierCoulombGauge->nDof()) );
        M_fieldLagrangeMultiplierCoulombGauge = M_spaceLagrangeMultiplierCoulombGauge->elementPtr( FieldTag::lagrangeMultiplierCoulombGauge(this).identifierString() );
    }

    double tElpased = this->timerTool("Constructor").stop("initFunctionSpaces");
    this->log("Magnetic","initFunctionSpaces",(boost::format("finish in %1% s")%tElpased).str() );
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("Magnetic","init", "start" );
    this->timerTool("Constructor").start();

    this->initModelProperties();

    // physics
    this->initPhysics( this->shared_from_this(), this->modelProperties().models() );

    this->initMaterialProperties();

    this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

#if 0
    for ( auto & [physicId,physicObj] : this->physicsFromCurrentType() )
        std::static_pointer_cast<ModelPhysicMagnetic<nDim>>(physicObj)->updateForUse( this->materialsProperties(), this->mesh() );
#endif

    this->initFunctionSpaces();

    this->initBoundaryConditions();

    // start or restart time step scheme
    if ( !this->isStationary() )
        this->initTimeStep();

    // post-process
    this->initPostProcess();

    // update constant parameters into expressions
    this->updateParameterValues();

    // update initial conditions
    this->updateInitialConditions( this->symbolsExpr() );

    // automatic solver selection
    if ( M_solverName == "automatic" )
    {
        bool isNonLinear = false;
#if 0
        auto mfields = this->modelFields();
        auto se = this->symbolsExpr( mfields );
        auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( 0/*rowStartInVector*/ ) );
        auto trialSymbolNames = tse.names();
        for ( std::string tsName : trialSymbolNames )
        {
            if ( this->materialsProperties()->hasThermalConductivityDependingOnSymbol( tsName ) )
            {
                isNonLinear = true;
                break;
            }
            for ( auto const& [bcName,bcData] : M_boundaryConditions->heatFlux() )
            {
                auto neumannExpr = bcData->expr();
                if ( neumannExpr.hasSymbolDependency( tsName, se ) )
                {
                    isNonLinear = true;
                    break;
                }
            }

            if ( isNonLinear )
                break;
        }
#endif
        M_solverName = isNonLinear? "Newton" : "Linear";
    }

    // algebraic model
    this->initAlgebraicModel();


    // // update constant parameters into (second time because some parameters can be defined from initial conditions or post-process measures)
    // this->updateParameterValues();


    // mesh adaptation at event after_init
    using mesh_adaptation_type = typename super_type::super_model_meshes_type::mesh_adaptation_type;
    this->template updateMeshAdaptation<mesh_type>( this->keyword(),
                                                    mesh_adaptation_type::createEvent<mesh_adaptation_type::Event::Type::after_init>(),
                                                    this->symbolsExpr() );

    // algebraic solver
    if ( buildModelAlgebraicFactory )
        this->initAlgebraicFactory();

    this->setIsUpdatedForUse( true );

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("Magnetic","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::initAlgebraicModel()
{
    // backend
    this->initAlgebraicBackend();

    // subspaces index
    size_type currentStartIndex = 0;
    this->setStartSubBlockSpaceIndex( FieldTag::vectorPotential(this).identifier(), currentStartIndex++ );
    if ( M_nullSpaceMethod == "saddle-point" )
        this->setStartSubBlockSpaceIndex( FieldTag::lagrangeMultiplierCoulombGauge(this).identifier(), currentStartIndex++ );
    size_type nBlock = this->startSubBlockSpaceIndices().size();

    this->updateAlgebraicDofEliminationIds();

    // vector solution
    auto bvs = this->initAlgebraicBlockVectorSolution( nBlock );
    bvs->operator()( this->startSubBlockSpaceIndex( FieldTag::vectorPotential(this).identifier() ) ) = this->fieldVectorPotentialPtr();
    if ( M_nullSpaceMethod == "saddle-point" )
        bvs->operator()( this->startSubBlockSpaceIndex( FieldTag::lagrangeMultiplierCoulombGauge(this).identifier() ) ) = this->fieldLagrangeMultiplierCoulombGaugePtr();
    // init petsc vector associated to the block
    bvs->buildVector( this->backend() );


    // InHousePreconditioner : Hypre-AMS
    this->initInHousePreconditioner();
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
MAGNETIC_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->startSubBlockSpaceIndices().size();//this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    size_type startVectorPotential = this->startSubBlockSpaceIndex( FieldTag::vectorPotential(this).identifier() );
    this->log("Magnetic","buildBlockMatrixGraph", fmt::format("start with nBlock: {} et startVectorPotential{}", nBlock, startVectorPotential ) );

    myblockGraph(startVectorPotential,startVectorPotential) = stencil(_test=this->spaceVectorPotential(),
                                                                      _trial=this->spaceVectorPotential() )->graph();
    if ( M_nullSpaceMethod == "saddle-point" )
    {
        size_type startLmCoulombGauge = this->startSubBlockSpaceIndex( FieldTag::lagrangeMultiplierCoulombGauge(this).identifier() );
        myblockGraph(startVectorPotential,startLmCoulombGauge) = stencil(_test=this->spaceVectorPotential(),
                                                                         _trial=this->spaceLagrangeMultiplierCoulombGauge() )->graph();
        myblockGraph(startLmCoulombGauge,startVectorPotential) = stencil(_test=this->spaceLagrangeMultiplierCoulombGauge(),
                                                                         _trial=this->spaceVectorPotential() )->graph();
    }
    myblockGraph.close();

    this->log("Magnetic","buildBlockMatrixGraph", "finish" );
    return myblockGraph;
}


MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::initInHousePreconditioner()
{
    if ( !M_preconditionerAttachAms )
        return;

    auto Xh = this->spaceVectorPotential();
    auto XhL = Pch<1>( Xh->mesh(), this->rangeMeshElements() );
    M_preconditionerAmsMatrixG = Grad( _domainSpace=XhL, _imageSpace=Xh).matPtr();

#if 1
    for ( int k=0 ; k<nRealDim ; ++k )
    {
        M_preconditionerAmsVectorOnes[k] = this->backend()->newVector(Xh);
        auto oneField = Xh->element( M_preconditionerAmsVectorOnes[k] );
        oneField.on(_range=this->rangeMeshElements(),_expr=one<nRealDim>(k),_close=true);
    }
#else

    if constexpr (nRealDim == 2 )
    {
        auto ozz = Xh->element();
        auto zoz = Xh->element();
        //auto zzo = Xh->element();
        ozz.on(_range=elements(Xh->mesh()),_expr=vec(cst(1),cst(0)/*,cst(0)*/));
        zoz.on(_range=elements(Xh->mesh()),_expr=vec(cst(0),cst(1)/*,cst(0)*/));
        //zzo.on(_range=elements(Xh->mesh()),_expr=vec(cst(0),cst(0),cst(1)));
        M_preconditionerAmsVectorOnes[0] = this->backend()->newVector(Xh); *M_preconditionerAmsVectorOnes[0] = ozz; M_preconditionerAmsVectorOnes[0]->close();
        M_preconditionerAmsVectorOnes[1] = this->backend()->newVector(Xh); *M_preconditionerAmsVectorOnes[1] = zoz; M_preconditionerAmsVectorOnes[1]->close();
    }
    else if constexpr (nRealDim == 3 )
    {
        auto ozz = Xh->element();
        auto zoz = Xh->element();
        auto zzo = Xh->element();
        ozz.on(_range=elements(Xh->mesh()),_expr=vec(cst(1),cst(0),cst(0)));
        zoz.on(_range=elements(Xh->mesh()),_expr=vec(cst(0),cst(1),cst(0)));
        zzo.on(_range=elements(Xh->mesh()),_expr=vec(cst(0),cst(0),cst(1)));
        M_preconditionerAmsVectorOnes[0] = this->backend()->newVector(Xh); *M_preconditionerAmsVectorOnes[0] = ozz; M_preconditionerAmsVectorOnes[0]->close();
        M_preconditionerAmsVectorOnes[1] = this->backend()->newVector(Xh); *M_preconditionerAmsVectorOnes[1] = zoz; M_preconditionerAmsVectorOnes[1]->close();
        M_preconditionerAmsVectorOnes[2] = this->backend()->newVector(Xh); *M_preconditionerAmsVectorOnes[2] = zzo; M_preconditionerAmsVectorOnes[2]->close();
    }

    // WARNING TODO
#endif
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::updateInHousePreconditioner( DataUpdateLinear & data ) const
{
#if 0
    if ( !M_preconditionerAttachPMM && !M_preconditionerAttachPCD )
        return;
    vector_ptrtype const& vecSol = data.currentSolution();
    this->updateInHousePreconditioner( data, this->modelContext( vecSol, this->rowStartInVector() ) );
#endif
}
MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::updateInHousePreconditioner( DataUpdateJacobian & data ) const
{
#if 0
    if ( !M_preconditionerAttachPMM && !M_preconditionerAttachPCD )
        return;
    vector_ptrtype const& vecSol = data.currentSolution();
    this->updateInHousePreconditioner( data, this->modelContext( vecSol, this->rowStartInVector() ) );
#endif
}


MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::applyRemesh( mesh_ptrtype oldMesh, mesh_ptrtype newMesh, std::shared_ptr<RemeshInterpolation> remeshInterp )
{
    this->log("Magnetic","applyRemesh", "start" );

    this->log("Magnetic","applyRemesh", "finish" );
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::initTimeStep()
{
#if 0
    this->log("Magnetic","initTimeStep", "start" );
    this->timerTool("Constructor").start();

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix

    int bdfOrder = 1;
    if ( M_timeStepping == "BDF" )
        bdfOrder = ioption(_prefix=this->prefix(),_name="bdf.order",_vm=this->clovm());
    int nConsecutiveSave = std::max( 3, bdfOrder ); // at least 3 is required when restart with theta scheme

    M_bdfTemperature = this->createBdf( this->spaceTemperature(),"temperature", bdfOrder, nConsecutiveSave, myFileFormat );

    if (!this->doRestart())
    {
        // up current time
        this->updateTime( M_bdfTemperature->timeInitial() );
    }
    else
    {
        // start time step
        double tir = M_bdfTemperature->restart();
        // load a previous solution as current solution
        *this->fieldTemperaturePtr() = M_bdfTemperature->unknown(0);
        // up initial time
        this->setTimeInitial( tir );
        // up current time
        this->updateTime( tir );
    }

    double tElapsed = this->timerTool("Constructor").stop("initTimeStep");
    this->log("Magnetic","initTimeStep", (boost::format("finish in %1% s") %tElapsed).str() );
#endif
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("Magnetic","initPostProcess", "start");
    this->timerTool("Constructor").start();

    this->setPostProcessExportsAllFieldsAvailable( { FieldTag::vectorPotential(this).identifierString(), "flux-density", "field-intensity" } );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );
    this->setPostProcessSaveAllFieldsAvailable( { FieldTag::vectorPotential(this).identifierString() } );
    super_type::initPostProcess();

    if ( !this->postProcessExportsFields().empty() )
    {
#if 0
        std::string geoExportType="static";//change_coords_only, change, static
#else
        bool useStaticExporter = boption(_name="exporter.use-static-mesh",_prefix=this->prefix(),_vm=this->clovm());
        std::string geoExportType = useStaticExporter? "static":"change";
#endif
        M_exporter = exporter( _mesh=this->mesh(),
                               _name="Export",
                               _geo=geoExportType,
                               _path=this->exporterPath() );

        // restart exporter
        if ( M_exporter->doExport() && this->doRestart() && this->restartPath().empty() )
            M_exporter->restart(this->timeInitial());
    }

#if 0
    if ( this->modelProperties().postProcess().hasJsonProperties( this->keyword() ) )
    {
        auto const& j_pp = this->modelProperties().postProcess().jsonProperties( this->keyword() );
        std::string ppTypeMeasures = "Measures";
        if ( j_pp.contains( ppTypeMeasures ) )
        {
            auto j_pp_measures = j_pp.at( ppTypeMeasures );
            for ( auto const& [j_pp_measureskey,j_pp_measuresval] : j_pp_measures.items() )
            {
                if ( j_pp_measureskey == "Normal-Heat-Flux" )
                {
                    for ( auto const& [j_pp_measures_nhfkey,j_pp_measures_nhfval] : j_pp_measuresval.items() )
                    {
                        auto indexesAllCases = ModelIndexes::generateAllCases( j_pp_measures_nhfval );
                        for ( auto const& indexes : indexesAllCases )
                        {
                            ModelMeasuresNormalFluxGeneric ppFlux;
                            ppFlux.setup( j_pp_measures_nhfval, indexes.replace( j_pp_measures_nhfkey ), indexes );
                            if ( !ppFlux.markers().empty() )
                                M_postProcessMeasuresNormalHeatFlux[ppFlux.name()] = ppFlux;
                        }
                    }
                }
            }
        }
    }
#endif
    auto se = this->symbolsExpr();
    this->template initPostProcessMeshes<mesh_type>( se );

    // start or restart the export of measures
    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasures().restart( this->timeInitial() );
    }

    double tElpased = this->timerTool("Constructor").stop("initPostProcess");
    this->log("Magnetic","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );

}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    auto algebraicFactory = std::make_shared<model_algebraic_factory_type>( this->shared_from_this(),this->backend() );
    this->setAlgebraicFactory( algebraicFactory );

    if ( M_preconditionerAttachAms )
    {
        this->algebraicFactory()->attachAuxiliarySparseMatrix( "G", M_preconditionerAmsMatrixG );
        this->algebraicFactory()->attachAuxiliaryVector( "Px", M_preconditionerAmsVectorOnes.at(0) );
        if ( M_preconditionerAmsVectorOnes.size() > 1 )
        this->algebraicFactory()->attachAuxiliaryVector( "Py", M_preconditionerAmsVectorOnes.at(1) );
        if ( M_preconditionerAmsVectorOnes.size() > 2 )
            this->algebraicFactory()->attachAuxiliaryVector( "Pz", M_preconditionerAmsVectorOnes.at(2) );

        //this->algebraicFactory()->attachAuxiliarySparseMatrix("a_beta",NULL);
    }

#if 0
    if ( M_timeStepping == "Theta" )
    {
        M_timeStepThetaSchemePreviousContrib = this->backend()->newVector( this->algebraicBlockVectorSolution()->vectorMonolithic()->mapPtr() );
        algebraicFactory->addVectorResidualAssembly( M_timeStepThetaSchemePreviousContrib, 1.0, "Theta-Time-Stepping-Previous-Contrib", true );
        algebraicFactory->addVectorLinearRhsAssembly( M_timeStepThetaSchemePreviousContrib, -1.0, "Theta-Time-Stepping-Previous-Contrib", false );
        if ( M_stabilizationGLS )
            algebraicFactory->dataInfos().addVectorInfo( "time-stepping.previous-solution", this->backend()->newVector( this->algebraicBlockVectorSolution()->vectorMonolithic()->mapPtr() ) );
    }
#endif
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

    super_physics_type::updateInformationObjectFromCurrentType( p["Physics"] );

#if 0
    // Physics
    subPt.emplace( "time mode", std::string( (this->isStationary())?"Stationary":"Transient") );
    p["Physics2"] = subPt;
#endif
    // Boundary Conditions
    M_boundaryConditions->updateInformationObject( p["Boundary Conditions"] );

    // Materials properties
    if ( this->materialsProperties() )
        this->materialsProperties()->updateInformationObject( p["Materials Properties"] );
    // FunctionSpace
    nl::json subPt;
    subPt.clear();
    subPt["VectorPotential"] = M_spaceVectorPotential->journalSection().to_string();
    p.emplace( "Function Spaces",  subPt );

    this->modelFields().updateInformationObject( p["Fields"] );

#if 0
    if ( !this->isStationary() )
    {
        subPt.clear();
        subPt.emplace( "initial time", this->timeStepBase()->timeInitial() );
        subPt.emplace( "final time", this->timeStepBase()->timeFinal() );
        subPt.emplace( "time step", this->timeStepBase()->timeStep() );
        subPt.emplace( "type", M_timeStepping );
        p["Time Discretization"] = subPt;
    }
#endif

    // Algebraic Solver
    if ( this->algebraicFactory() )
    {
        this->algebraicFactory()->updateInformationObject( p["Algebraic Solver"] );
    }
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
MAGNETIC_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics") )
        tabInfo->add( "Physics", super_physics_type::tabulateInformations( jsonInfo.at("Physics"), tabInfoProp ) );
#if 0
    if ( jsonInfo.contains("Physics2") )
    {
        Feel::Table tabInfoPhysics;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoPhysics, jsonInfo.at("Physics2"), tabInfoProp );
        tabInfo->add( "Physics2", TabulateInformations::New( tabInfoPhysics, tabInfoProp ) );
    }
#endif
    if ( this->materialsProperties() && jsonInfo.contains("Materials Properties") )
        tabInfo->add( "Materials Properties", this->materialsProperties()->tabulateInformations(jsonInfo.at("Materials Properties"), tabInfoProp ) );

    if ( jsonInfo.contains("Boundary Conditions") )
        tabInfo->add( "Boundary Conditions", boundary_conditions_type::tabulateInformations( jsonInfo.at("Boundary Conditions"), tabInfoProp ) );

    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );

    if ( jsonInfo.contains("Function Spaces") )
    {
        auto const& jsonInfoFunctionSpaces = jsonInfo.at("Function Spaces");
        auto tabInfoFunctionSpaces = TabulateInformationsSections::New( tabInfoProp );

        nl::json::json_pointer jsonPointerSpaceVectorPotential( jsonInfoFunctionSpaces.at( "VectorPotential" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerSpaceVectorPotential ) )
            tabInfoFunctionSpaces->add( "VectorPotential", TabulateInformationTools::FromJSON::tabulateInformationsFunctionSpace( JournalManager::journalData().at( jsonPointerSpaceVectorPotential ), tabInfoProp ) );

        tabInfo->add( "Function Spaces", tabInfoFunctionSpaces );
    }

    // fields
    if ( jsonInfo.contains("Fields") )
        tabInfo->add( "Fields", TabulateInformationTools::FromJSON::tabulateInformationsModelFields( jsonInfo.at("Fields"), tabInfoProp ) );
#if 0
    if ( jsonInfo.contains("Time Discretization") )
    {
        Feel::Table tabInfoTimeDiscr;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoTimeDiscr, jsonInfo.at("Time Discretization"), tabInfoProp );
        tabInfo->add( "Time Discretization", TabulateInformations::New( tabInfoTimeDiscr, tabInfoProp ) );
    }

    if ( jsonInfo.contains("Finite element stabilization") )
    {
        Feel::Table tabInfoStab;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoStab, jsonInfo.at("Finite element stabilization"), tabInfoProp );
        tabInfo->add( "Finite element stabilization", TabulateInformations::New( tabInfoStab, tabInfoProp ) );
    }
#endif
    if ( jsonInfo.contains( "Algebraic Solver" ) )
        tabInfo->add( "Algebraic Solver", model_algebraic_factory_type::tabulateInformations( jsonInfo.at("Algebraic Solver"), tabInfoProp ) );

    return tabInfo;
}


MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    int previousParam = 0;
    while ( true )
    {
        this->modelProperties().parameters().updateParameterValues();
        auto paramValues = this->modelProperties().parameters().toParameterValues();
        this->materialsProperties()->updateParameterValues( paramValues );
        for ( auto [physicName,physicData] : this->physics/*FromCurrentType*/() )
            physicData->updateParameterValues( paramValues );

        this->updateParameterValues_postProcess( paramValues, prefixvm("postprocess",this->keyword(),"_" ) );

        if ( paramValues.size() == previousParam )
            break;
        previousParam = paramValues.size();

        this->setParameterValues( paramValues );
    }
}
MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    this->log("Magnetic","setParameterValues", "start");

    for ( auto const& [param,val] : paramValues )
        M_currentParameterValues[param] = val;

    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->modelProperties().initialConditions().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        physicData->setParameterValues( paramValues );

    super_type::super_model_meshes_type::setParameterValues( paramValues );

    M_boundaryConditions->setParameterValues( paramValues );

    this->log("Magnetic","setParameterValues", "finish");
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    M_boundaryConditions = std::make_shared<boundary_conditions_type>( this->shared_from_this() );
    if ( !this->modelProperties().boundaryConditions().hasSection( this->keyword() ) )
        return;
    M_boundaryConditions->setup( this->modelProperties().boundaryConditions().section( this->keyword() ) );
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::updateAlgebraicDofEliminationIds()
{
    for ( auto const& [bcName,bcData] : M_boundaryConditions->magneticPotentialImposed() )
        bcData->updateDofEliminationIds( *this, FieldTag::vectorPotential(this).identifierString(), this->spaceVectorPotential() );
    for ( auto const& [bcName,bcData] : M_boundaryConditions->magneticInsulation() )
        bcData->updateDofEliminationIds( *this, FieldTag::vectorPotential(this).identifierString(), this->spaceVectorPotential() );

    if ( M_nullSpaceMethod == "saddle-point" )
    {
        this->updateDofEliminationIds( FieldTag::lagrangeMultiplierCoulombGauge(this).identifierString(),
                                       this->spaceLagrangeMultiplierCoulombGauge(),
                                       boundaryfaces( support( this->spaceLagrangeMultiplierCoulombGauge() ) )
                                       );
    }
}

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("Magnetic","solve", "start");
    this->timerTool("Solve").start();

    this->setStartBlockSpaceIndex( 0 );

    this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();

    this->algebraicFactory()->solve( M_solverName, this->algebraicBlockVectorSolution()->vectorMonolithic() );

    this->algebraicBlockVectorSolution()->localize();

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("Magnetic","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}



MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    this->exportResults( time, mfields, se, this->exprPostProcessExports( se ) );
}
#if 0
MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
void
MAGNETIC_CLASS_TEMPLATE_TYPE::executePostProcessMeasures( double time )
{
    auto mfields = this->modelFields();
    this->executePostProcessMeasures( time, mfields, this->symbolsExpr( mfields ) );
}
#endif

MAGNETIC_CLASS_TEMPLATE_DECLARATIONS
bool
MAGNETIC_CLASS_TEMPLATE_TYPE::checkResults() const
{
    const_cast<self_type*>(this)->updateParameterValues();
    auto se = this->symbolsExpr();
    return super_type::checkResults( se );
}


} // end namespace FeelModels
} // end namespace Feel
