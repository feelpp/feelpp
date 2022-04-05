/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>

//#include <feel/feelmodels/modelmesh/createmesh.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameter.hpp>

namespace Feel
{
namespace FeelModels
{

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::CoefficientFormPDE( typename super_type::super2_type::infos_ptrtype const& infosPDE,
                                                            std::string const& prefix,
                                                            std::string const& keyword,
                                                            worldcomm_ptr_t const& worldComm,
                                                            std::string const& subPrefix,
                                                            ModelBaseRepository const& modelRep )
    :
    super_type( infosPDE, prefix, keyword, worldComm, subPrefix, modelRep ),
    ModelBase( prefix, keyword, worldComm, subPrefix, modelRep, ModelBaseCommandLineOptions( coefficientformpde_options( prefix ) ) )
{}


COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    this->log("CoefficientFormPDE","init", "start" );
    this->timerTool("Constructor").start();

    this->initModelProperties();

    this->initPhysics( this->shared_from_this(), this->modelProperties().models() );

    this->initMaterialProperties();

    this->initMesh();

    this->materialsProperties()->addMesh( this->mesh() );

    this->initFunctionSpaces();

    this->initBoundaryConditions();

    // start or restart time step scheme
    if ( !this->isStationary() )
        this->initTimeStep();

    // stabilization gls
    if ( this->M_applyStabilization && !this->M_stabilizationGLSParameter )
    {
        typedef StabilizationGLSParameter<mesh_type, nOrderUnknown> stab_gls_parameter_impl_type;
        this->M_stabilizationGLSParameter.reset( new stab_gls_parameter_impl_type( this->mesh(), prefixvm(this->prefix(),"stabilization.gls.parameter"), this->clovm() ) );
        this->M_stabilizationGLSParameter->init();
    }

    // update constant parameters into
    this->updateParameterValues();

    // update initial conditions
    this->updateInitialConditions( this->symbolsExpr() );

    // post-process
    this->initPostProcess();

    // backend
    this->initAlgebraicBackend();

    // subspaces index
    size_type currentStartIndex = 0;
    this->setStartSubBlockSpaceIndex( this->unknownName(), currentStartIndex++ );

     // vector solution
    auto bvs = this->initAlgebraicBlockVectorSolution( 1 );
    bvs->operator()(0) = this->fieldUnknownPtr();
    // init petsc vector associated to the block
    bvs->buildVector( this->backend() );

    // solver type
    M_solverName = soption(_prefix=this->prefix(),_name="solver",_vm=this->clovm());
    if ( M_solverName == "automatic" )
        this->updateAutomaticSolverSelection();

    // algebraic solver
    if ( buildModelAlgebraicFactory )
        this->initAlgebraicFactory();

    this->setIsUpdatedForUse( true );

    double tElapsedInit = this->timerTool("Constructor").stop("init");
    if ( this->scalabilitySave() ) this->timerTool("Constructor").save();
    this->log("CoefficientFormPDE","init",(boost::format("finish in %1% s")%tElapsedInit).str() );
}




COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initFunctionSpaces()
{
    this->log("CoefficientFormPDE","initFunctionSpaces", "start" );
    this->timerTool("Constructor").start();

    if( !M_Xh )
    {
        bool useExtendedDoftable = is_lagrange_polynomialset_v<typename space_unknown_type::fe_type> && !space_unknown_type::fe_type::isContinuous;

        auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
        // functionspace
        if ( mom->isDefinedOnWholeMesh( this->physicsAvailableFromCurrentType() ) )
        {
            this->M_rangeMeshElements = elements(this->mesh());
            M_Xh = space_unknown_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),_extended_doftable=useExtendedDoftable );
        }
        else
        {
            this->M_rangeMeshElements = markedelements(this->mesh(), mom->markers( this->physicsAvailableFromCurrentType() ));
            M_Xh = space_unknown_type::New( _mesh=this->mesh(), _worldscomm=this->worldsComm(),_range=this->M_rangeMeshElements,_extended_doftable=useExtendedDoftable );
        }
    }
    else
    {
        this->M_rangeMeshElements = elements( support( M_Xh ) );
    }

    M_fieldUnknown.reset( new element_unknown_type( M_Xh,this->unknownName() ) );

    double tElpased = this->timerTool("Constructor").stop("initFunctionSpaces");
    this->log("CoefficientFormPDE","initFunctionSpaces",(boost::format("finish in %1% s")%tElpased).str() );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initBoundaryConditions()
{
    M_boundaryConditions = std::make_shared<boundary_conditions_type>( this->shared_from_this() );
    if ( !this->modelProperties().boundaryConditions().hasSection( this->keyword() ) )
        return;

    M_boundaryConditions->setup( this->modelProperties().boundaryConditions().section( this->keyword() ) );

    for ( auto const& [bcName,bcData] : M_boundaryConditions->dirichlet() )
        bcData->updateDofEliminationIds( *this, this->unknownName(), this->spaceUnknown() );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    this->log("CoefficientFormPDE","initTimeStep", "start" );
    this->timerTool("Constructor").start();

    CHECK( this->timeStepping() == "BDF" || this->timeStepping() == "Theta" ) << "invalid time-stepping";

    std::string myFileFormat = soption(_name="ts.file-format");// without prefix

    int bdfOrder = 1;
    if ( this->timeStepping() == "BDF" )
        bdfOrder = ioption(_prefix=this->prefix(),_name="bdf.order",_vm=this->clovm());
    int nConsecutiveSave = std::max( 3, bdfOrder ); // at least 3 is required when restart with theta scheme

    M_bdfUnknown = this->createBdf( this->spaceUnknown(),this->unknownName(), bdfOrder, nConsecutiveSave, myFileFormat );

    if (!this->doRestart())
    {
        // up current time
        this->updateTime( M_bdfUnknown->timeInitial() );
    }
    else
    {
        // start time step
        double tir = M_bdfUnknown->restart();
        // load a previous solution as current solution
        *this->fieldUnknownPtr() = M_bdfUnknown->unknown(0);
        // up initial time
        this->setTimeInitial( tir );
        // up current time
        this->updateTime( tir );
    }

    double tElapsed = this->timerTool("Constructor").stop("initTimeStep");
    this->log("CoefficientFormPDE","initTimeStep", (boost::format("finish in %1% s") %tElapsed).str() );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initPostProcess()
{
    this->log("CoefficientFormPDE","initPostProcess", "start");
    this->timerTool("Constructor").start();

    this->initBasePostProcess();

    auto se = this->symbolsExpr();
    this->template initPostProcessMeshes<mesh_type>( se );

    double tElpased = this->timerTool("Constructor").stop("initPostProcess");
    this->log("CoefficientFormPDE","initPostProcess",(boost::format("finish in %1% s")%tElpased).str() );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::initAlgebraicFactory()
{
    auto algebraicFactory = std::make_shared<typename ModelAlgebraic::model_algebraic_factory_type>( this->shared_from_this(),this->backend() );
    this->setAlgebraicFactory( algebraicFactory );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateAutomaticSolverSelection()
{
    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( 0/*rowStartInVector*/ ) );
    auto trialSymbolNames = tse.names();
    bool isLinear = true;

    if ( this->hasSymbolDependencyInCoefficients( trialSymbolNames, se ) 
            || this->hasSymbolDependencyInBoundaryConditions( trialSymbolNames, se ) 
            || ( this->applyStabilization() && this->stabilizationGLS_applyShockCapturing() ) )
    {
        isLinear = false;
    }

    M_solverName = ( isLinear )? "Linear" : "Newton";
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = 1;//this->nBlockMatrixGraph();

    size_type thePat = size_type(Pattern::COUPLED);
    if ( is_lagrange_polynomialset_v<typename space_unknown_type::fe_type> && !space_unknown_type::fe_type::isContinuous )
        thePat = size_type(Pattern::EXTENDED);
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    myblockGraph(0,0) = stencil(_test=this->spaceUnknown(),
                                _trial=this->spaceUnknown(),
                                _pattern=thePat )->graph();
    return myblockGraph;
}


COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

    super_type::super_physics_type::updateInformationObjectFromCurrentType( p["Physics"] );

    // Boundary Conditions
    M_boundaryConditions->updateInformationObject( p["Boundary Conditions"] );

    // FunctionSpace
    nl::json subPt;
    subPt[this->unknownName()] = M_Xh->journalSection().to_string();
    p.emplace( "Function Spaces",  subPt );

    this->modelFields().updateInformationObject( p["Fields"] );

#if 0
    if ( M_algebraicFactory )
        M_algebraicFactory->updateInformationObject( p["Algebraic Solver"] );
#endif
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics") )
        tabInfo->add( "Physics", super_type::super_physics_type::tabulateInformations( jsonInfo.at("Physics"), tabInfoProp ) );

    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );

    if ( jsonInfo.contains("Boundary Conditions") )
        tabInfo->add( "Boundary Conditions", boundary_conditions_type::tabulateInformations( jsonInfo.at("Boundary Conditions"), tabInfoProp ) );

    if ( jsonInfo.contains("Function Spaces") )
    {
        auto const& jsonInfoFunctionSpaces = jsonInfo.at("Function Spaces");
        auto tabInfoFunctionSpaces = TabulateInformationsSections::New( tabInfoProp );
        nl::json::json_pointer jsonPointerSpaceUnknown( jsonInfoFunctionSpaces.at( this->unknownName() ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerSpaceUnknown ) )
            tabInfoFunctionSpaces->add( this->unknownName(), TabulateInformationTools::FromJSON::tabulateInformationsFunctionSpace( JournalManager::journalData().at( jsonPointerSpaceUnknown ), tabInfoProp ) );
        tabInfo->add( "Function Spaces", tabInfoFunctionSpaces );
    }

    if ( jsonInfo.contains("Fields") )
        tabInfo->add( "Fields", TabulateInformationTools::FromJSON::tabulateInformationsModelFields( jsonInfo.at("Fields"), tabInfoProp ) );

    return tabInfo;
}

#if 0
COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
#if 0
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||-----------Info : CoefficientFormPDE----------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository();
    *_ostr << "\n   Physical Model"
           << "\n     -- name      : " << this->physic()
           << "\n     -- time mode : " << std::string( (this->isStationary())?"Stationary":"Transient");
    *_ostr << "\n   Boundary conditions"
           << M_bcDirichletMarkerManagement.getInfoDirichletBC()
           << M_bcNeumannMarkerManagement.getInfoNeumannBC()
           << M_bcRobinMarkerManagement.getInfoRobinBC();
    *_ostr << this->materialsProperties()->getInfoMaterialParameters()->str();
#if 0
    *_ostr << "\n   Mesh Discretization"
           << "\n     -- mesh filename      : " << this->meshFile()
           << "\n     -- number of element : " << this->mesh()->numGlobalElements()
           << "\n     -- order             : " << nOrderGeo;
#endif
    *_ostr << "\n   Space Unknown Discretization"
           << "\n     -- name of unknown         : " << this->unknownName()
           << "\n     -- symbol of unknown : " << this->unknownSymbol()
           << "\n     -- basis : " << this->unknownBasis()
           << "\n     -- number of dof : " << this->spaceUnknown()->nDof() << " (" << this->spaceUnknown()->nLocalDof() << ")";
    if ( !this->isStationary() )
    {
        *_ostr << "\n   Time Discretization"
               << "\n     -- initial time : " << this->timeStepBase()->timeInitial()
               << "\n     -- final time   : " << this->timeStepBase()->timeFinal()
               << "\n     -- time step    : " << this->timeStepBase()->timeStep()
               << "\n     -- type : " << this->timeStepping();
        if ( this->timeStepping() == "BDF" )
            *_ostr << " ( order=" << this->timeStepBdfUnknown()->timeOrder() << " )";
    }
    if ( this->algebraicFactory() )
        *_ostr << this->algebraicFactory()->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";
#endif
    return _ostr;
}
#endif

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    this->log("CoefficientFormPDE","startTimeStep", "start");

    // start time step
    if (!this->doRestart())
        M_bdfUnknown->start( M_bdfUnknown->unknowns() );
     // up current time
    this->updateTime( M_bdfUnknown->time() );

    // update all expressions in bc or in house prec
    this->updateParameterValues();

    this->log("CoefficientFormPDE","startTimeStep", "finish");
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->log("CoefficientFormPDE","updateTimeStep", "start");
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    bool rebuildCstAssembly = false;
    if ( this->timeStepping() == "BDF" )
    {
        int previousTimeOrder = this->timeStepBdfUnknown()->timeOrder();
        M_bdfUnknown->next( this->fieldUnknown() );
        int currentTimeOrder = this->timeStepBdfUnknown()->timeOrder();
        rebuildCstAssembly = previousTimeOrder != currentTimeOrder && this->timeStepBase()->strategy() == TS_STRATEGY_DT_CONSTANT;
        this->updateTime( this->timeStepBdfUnknown()->time() );
    }
    else if ( this->timeStepping() == "Theta" )
    {
        M_bdfUnknown->next( this->fieldUnknown() );
        this->updateTime( this->timeStepBdfUnknown()->time() );
    }

    if ( rebuildCstAssembly )
        this->setNeedToRebuildCstPart(true);

    this->updateParameterValues();

    this->timerTool("TimeStepping").stop("updateTimeStep");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("CoefficientFormPDE","updateTimeStep", "finish");
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );

    this->setParameterValues( paramValues );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    this->log("CoefficientFormPDE","setParameterValues", "start");

    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->modelProperties().initialConditions().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        physicData->setParameterValues( paramValues );

    M_boundaryConditions->setParameterValues( paramValues );

    this->log("CoefficientFormPDE","setParameterValues", "finish");
}


COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("CoefficientFormPDE","solve", "start");
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
    this->log("CoefficientFormPDE","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto mctx = this->modelContext( XVec, this->rowStartInVector() );
    if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    {
        auto previousSol = data.vectorInfo( "time-stepping.previous-solution");
        auto mctxPrevious = this->modelContextNoTrialSymbolsExpr( previousSol, this->rowStartInVector() );
        mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    }
    this->updateLinearPDE( data, mctx );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto mctx = this->modelContext( XVec, this->rowStartInVector() );
    this->updateLinearPDEDofElimination( data, mctx );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    vector_ptrtype& U = data.initialGuess();
    auto mctx = this->modelContext( U, this->rowStartInVector() );
    this->updateNewtonInitialGuess( data, mctx );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto mctx = this->modelContext( XVec, this->rowStartInVector() );
    this->updateJacobian( data, mctx );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( ModelAlgebraic::DataUpdateJacobian & data ) const
{
    if ( !M_boundaryConditions->hasTypeDofElimination() )
        return;

    this->log("CoefficientFormPDE","updateJacobianDofElimination","start" );

    this->updateDofEliminationIds( this->unknownName(), data );

    this->log("CoefficientFormPDE","updateJacobianDofElimination","finish" );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    auto mctx = this->modelContext( XVec, this->rowStartInVector() );
    if ( data.hasVectorInfo( "time-stepping.previous-solution" ) )
    {
        auto previousSol = data.vectorInfo( "time-stepping.previous-solution");
        auto mctxPrevious = this->modelContextNoTrialSymbolsExpr( previousSol, this->rowStartInVector() );
        mctx.setAdditionalContext( "time-stepping.previous-model-context", std::move( mctxPrevious ) );
    }
    this->updateResidual( data, mctx );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( ModelAlgebraic::DataUpdateResidual & data ) const
{
    if ( !M_boundaryConditions->hasTypeDofElimination() )
        return;
    this->log("CoefficientFormPDE","updateResidualDofElimination","start" );

    this->updateDofEliminationIds( this->unknownName(), data );

    this->log("CoefficientFormPDE","updateResidualDofElimination","finish" );
}

COEFFICIENTFORMPDE_CLASS_TEMPLATE_DECLARATIONS
void
COEFFICIENTFORMPDE_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
#if 0
    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    this->exportResults( time, mfields, se, this->exprPostProcessExports( se ) );
#endif
}

} // namespace Feel
} // namespace FeelModels
